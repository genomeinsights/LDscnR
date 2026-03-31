#' Create a GDS file from a genotype matrix
#'
#' Converts a genotype matrix and accompanying SNP map into a GDS file
#' compatible with \pkg{SNPRelate}, and returns an open GDS handle.
#'
#' The genotype matrix is expected in \strong{individuals x SNPs} format.
#' SNP metadata are taken from \code{map}.
#'
#' @param geno Numeric genotype matrix of allele counts (0,1,2) with individuals in rows and SNPs in columns. NA's are allowed.
#' @param map Data frame or \code{data.table} with one row per SNP. Must contain
#'   the columns \code{marker}, \code{Chr}, and \code{Pos}.
#' @param gds_path File path where the GDS file will be written.
#'
#' @return An open GDS object.
#'
#' @details
#' Internally, the genotype matrix is transposed because
#' \code{SNPRelate::snpgdsCreateGeno()} expects SNPs in rows when
#' \code{snpfirstdim = TRUE}.
#'
#' Sample IDs are generated automatically as \code{ind_1}, \code{ind_2}, etc.
#'
#' @export
create_gds_from_geno <- function(geno, map, gds_path) {
  stopifnot(ncol(geno) == nrow(map))
  SNPRelate::snpgdsCreateGeno(
    gds_path,
    genmat         = t(round(geno)),     # SNPRelate expects SNP × sample if snpfirstdim=TRUE
    sample.id      = paste0("ind_", seq_len(nrow(geno))), ## sample name not important
    snp.id         = map$marker,
    snp.chromosome = map$Chr,
    snp.position   = map$Pos,
    snpfirstdim    = TRUE
  )
  SNPRelate::snpgdsOpen(gds_path)
}


.read_gds_ids <- function(gds) {
  list(
    snp_id  = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.id")),
    snp_chr = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.chromosome")),
    snp_pos = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.position"))
  )
}
#' Build an LD edge list from a subset of SNPs
#'
#' Computes pairwise linkage disequilibrium (LD; \eqn{r^2}) for a selected
#' subset of SNPs in a GDS file and returns the result as an edge list with
#' genomic coordinates and physical distances.
#'
#' Depending on the chosen options, the function can return:
#' \itemize{
#'   \item a standard pairwise edge list, with each SNP pair represented once,
#'   \item a SNP-centered symmetric representation, where each edge contributes
#'   one record to each SNP,
#'   \item an optionally symmetry-adjusted SNP-centered edge list for local
#'   neighborhood summaries.
#' }
#'
#' @param gds An open GDS object.
#' @param idx Optional integer vector of SNP indices (1-based, referring to the
#'   SNP order in the GDS file). If \code{NULL}, all SNPs are used.
#' @param slide_win_ld Integer SNP window size passed to
#'   \code{SNPRelate::snpgdsLDMat()}. If \code{> 0}, LD is computed within a
#'   sliding window of this size. If \code{<= 0}, all pairwise LD values are
#'   computed.
#' @param method LD-method paste on to \code{SNPRelate::snpgdsLDMat}, default="r"
#'   (EM algorithm assuming HWE). If strong deviations from HWE are expected use "corr".
#' @param cores Number of CPU threads used in LD computation.
#' @param by_chr Logical; if \code{TRUE}, LD is computed separately within each
#'   chromosome. If \code{FALSE}, LD is computed across the full SNP subset.
#'
#' @return
#' A \code{data.table} with columns:
#' \describe{
#'   \item{SNP1, SNP2}{SNP identifiers.}
#'   \item{Chr1, Chr2}{Chromosomes of the two SNPs.}
#'   \item{pos1, pos2}{Physical positions of the two SNPs.}
#'   \item{r2}{Pairwise LD (\eqn{r^2}).}
#'   \item{d}{Absolute physical distance between the SNPs.}
#' }
#'
#' @details
#' When \code{by_chr = TRUE}, SNPs are split by chromosome and LD is computed
#' independently within each chromosome subset before results are combined.
#'
#' The SNP-centered representation is useful for downstream summaries that
#' require all LD partners of each focal SNP.
#'
#' @export
get_el <- function (gds, idx = NULL, SNP_id = NULL, slide_win_ld = 1000,
                    method = "r", cores = 1, by_chr = FALSE)
{
  ids <- .read_gds_ids(gds)
  if (is.null(SNP_id)) {
    if (missing(idx) || is.null(idx)) {
      snp_ids <- ids$snp_id
    } else {
      snp_ids <- ids$snp_id[as.integer(idx)]
    }
  } else {
    snp_ids <- SNP_id
  }
  slide <- if (slide_win_ld > 0) as.integer(slide_win_ld) else -1L

  compute_one <- function(local_idx) {
    ldmat <- SNPRelate::snpgdsLDMat(gds, snp.id = ids$snp_id[local_idx],
                                    method = method, slide = slide, verbose = FALSE,
                                    num.thread = as.integer(cores))
    el <- data.table::as.data.table(reshape2::melt(ldmat$LD^2,

                                                   value.name = "r2"))
    if (slide_win_ld > 0) {
      el[, `:=`(Var1, Var1 + Var2)]
    }
    el <- el[Var1 > Var2]
    el <- el[is.finite(r2)]
    el[, `:=`(Chr1, ids$snp_chr[local_idx][Var1])]
    el[, `:=`(Chr2, ids$snp_chr[local_idx][Var2])]
    el[, `:=`(pos1, ids$snp_pos[local_idx][Var1])]
    el[, `:=`(pos2, ids$snp_pos[local_idx][Var2])]
    el[, `:=`(SNP1, ids$snp_id[local_idx][Var1])]
    el[, `:=`(SNP2, ids$snp_id[local_idx][Var2])]
    el[, `:=`(d, abs(pos1 - pos2))]
    el[, .(SNP1, SNP2, Chr1, Chr2, pos1, pos2, r2, d)]
  }
  if (!by_chr) {
    return(compute_one(which(ids$snp_id %in% snp_ids)))
  }
  chr_vec <- ids$snp_chr
  chr_levels <- unique(chr_vec)
  el_list <- vector("list", length(chr_levels))
  # i <- 1
  for (i in seq_along(chr_levels)) {
    ch <- chr_levels[i]

    local_idx <- which(chr_vec == ch & ids$snp_id %in% snp_ids)

    if (length(local_idx) < 2L)
      next
    el_list[[i]] <- compute_one(local_idx)
  }

  data.table::rbindlist(el_list, use.names = TRUE)
}


.get_n_inds <- function(gds) {
  length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "sample.id")))
}
