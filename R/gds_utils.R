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
#' @param cores Number of CPU threads used in LD computation.
#' @param by_chr Logical; if \code{TRUE}, LD is computed separately within each
#'   chromosome. If \code{FALSE}, LD is computed across the full SNP subset.
#' @param symmetric Logical; if \code{TRUE}, returns a SNP-centered edge list in
#'   which each pair contributes one row for each of the two SNPs.
#' @param edge_symmetry Logical; if \code{TRUE}, applies an additional balancing
#'   step to the SNP-centered representation to reduce asymmetry in left/right
#'   genomic neighborhoods in chromosome edges (removes edge effects).
#'
#' @return
#' If \code{symmetric = FALSE} and \code{edge_symmetry = FALSE}, a
#' \code{data.table} with columns:
#' \describe{
#'   \item{Var1, Var2}{Indices of the paired SNPs within the analyzed subset.}
#'   \item{r2}{Pairwise LD (\eqn{r^2}).}
#'   \item{Chr1, Chr2}{Chromosomes of the two SNPs.}
#'   \item{pos1, pos2}{Physical positions of the two SNPs.}
#'   \item{SNP1, SNP2}{SNP identifiers.}
#'   \item{d}{Absolute physical distance between the SNPs.}
#' }
#'
#' If \code{symmetric = TRUE} or \code{edge_symmetry = TRUE}, a SNP-centered
#' \code{data.table} with columns:
#' \describe{
#'   \item{SNP}{Reference SNP.}
#'   \item{pos1}{Position of the reference SNP.}
#'   \item{pos2}{Position of the paired SNP.}
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
get_el <- function(gds,
                   idx = NULL,
                   slide_win_ld = 1000,
                   cores = 1,
                   by_chr = FALSE,
                   symmetric = FALSE,
                   edge_symmetry = FALSE) {

  ids <- .read_gds_ids(gds)

  if (missing(idx) || is.null(idx)) {
    idx <- seq_along(ids$snp_id)
  } else {
    idx <- as.integer(idx)
  }

  if (anyNA(idx) || length(idx) == 0L)
    stop("`idx` must be a non-empty integer vector.")
  if (min(idx) < 1L || max(idx) > length(ids$snp_id))
    stop("`idx` is out of range for this GDS.")

  slide <- if (slide_win_ld > 0) as.integer(slide_win_ld) else -1L

  # ------------------------------------------------------------
  # Helper to compute LD for one index subset
  # ------------------------------------------------------------
  compute_one <- function(local_idx) {

    snp_ids <- ids$snp_id[local_idx]

    ldmat <- SNPRelate::snpgdsLDMat(
      gds,
      snp.id = snp_ids,
      method = "r",
      slide = slide,
      verbose = FALSE,
      num.thread = as.integer(cores)
    )

    el <- data.table::as.data.table(
      reshape2::melt(ldmat$LD^2, value.name = "r2")
    )

    if (slide_win_ld > 0) {
      el[, Var1 := Var1 + Var2]
    }

    el <- el[Var1 > Var2]
    el <- el[is.finite(r2)]

    el[, Chr1 := ids$snp_chr[local_idx][Var1]]
    el[, Chr2 := ids$snp_chr[local_idx][Var2]]
    el[, pos1 := ids$snp_pos[local_idx][Var1]]
    el[, pos2 := ids$snp_pos[local_idx][Var2]]
    el[, SNP1 := ids$snp_id[local_idx][Var1]]
    el[, SNP2 := ids$snp_id[local_idx][Var2]]
    el[, d := abs(pos1 - pos2)]

    el
  }

  # ------------------------------------------------------------
  # by_chr = FALSE → original behavior
  # ------------------------------------------------------------
  if (!by_chr) {
    return(compute_one(idx))
  }

  # ------------------------------------------------------------
  # by_chr = TRUE → split per chromosome
  # ------------------------------------------------------------
  chr_vec <- ids$snp_chr[idx]
  chr_levels <- unique(chr_vec)

  el_list <- vector("list", length(chr_levels))

  for (i in seq_along(chr_levels)) {
    ch <- chr_levels[i]
    local_idx <- idx[chr_vec == ch]

    if (length(local_idx) < 2L)
      next

    el <-  compute_one(local_idx)

    if(symmetric || edge_symmetry){

      el <- data.table::rbindlist(list(
        el[, .(SNP = SNP1, pos = pos1, pos_other = pos2, r2, d)],
        el[, .(SNP = SNP2, pos = pos2, pos_other = pos1, r2, d)]
      ))

      # if(edge_symmetry){
      #
      #   el[, side := data.table::fifelse(pos_other > pos, "R", "L")]
      #
      #   maxL <- if (any(el$side == "L")) max(el$d[el$side == "L"]) else 0
      #   maxR <- if (any(el$side == "R")) max(el$d[el$side == "R"]) else 0
      #   S <- min(maxL, maxR)
      #
      #   if (S > 0 && maxL != maxR) {
      #     long_side <- if (maxL > maxR) "L" else "R"
      #     tail_el <- el[side == long_side & d > S]
      #     if (nrow(tail_el) > 0)
      #       el <- data.table::rbindlist(list(el, tail_el))
      #   }
      #
      # }

     el <- el[,.(SNP,pos1=pos,pos2=pos_other,r2,d)]

    }

    el_list[[i]] <- el
  }



  data.table::rbindlist(el_list, use.names = TRUE)
}


.get_n_inds <- function(gds) {
  length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "sample.id")))
}
