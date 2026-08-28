#' Create a GDS file from a genotype matrix
#'
#' Converts a genotype matrix and accompanying SNP map into a GDS file
#' compatible with \pkg{SNPRelate}, and returns an open GDS handle.
#'
#' The genotype matrix is expected in \strong{individuals x SNPs} format.
#' SNP metadata are taken from \code{map}.
#'
#' @param geno Numeric genotype matrix of allele counts (0,1,2) with individuals in rows and SNPs in columns. NA's are allowed. Rownames must contain individual ID names.
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
#' @param SNP_id Optional vector of SNP ids to restrict the edge list to
#'   (alternative to \code{idx}).
#' @param el_floor Minimum \eqn{r^2} for a pair to be retained, applied
#'   \strong{during} construction rather than afterwards. Default \code{0},
#'   which returns every finite pair and so reproduces the previous behaviour
#'   exactly.
#'
#'   The point is peak memory, not the returned object: a sliding window of
#'   \code{slide_win_ld} over \eqn{n} markers materialises on the order of
#'   \eqn{n \times} \code{slide_win_ld} rows before any downstream filter sees
#'   them, which for a 15,000-marker chromosome at the default window is ~15M
#'   rows and on the order of a gigabyte. A caller that is going to discard
#'   low-\eqn{r^2} pairs anyway should say so here and never pay for them.
#'
#'   \strong{Only safe when the caller genuinely does not need the discarded
#'   pairs.} Anything that summarises the LD \emph{distribution} does need them:
#'   [compute_LD_decay()] fits decay across the whole range, and
#'   [compute_ld_w()] takes a median within a window, which a floor would bias
#'   upward. Both therefore leave this at \code{0}. Threshold-based consumers
#'   such as [ld_complexity_reduction()], which joins at
#'   \code{ld_from_rho(b, c, rho)}, can safely floor anywhere below their own
#'   threshold.
#' @export
get_el <- function (gds, idx = NULL, SNP_id = NULL, slide_win_ld = 1000,
                    method = "r", cores = 1, by_chr = FALSE, el_floor = 0)
{
  if (!is.numeric(el_floor) || length(el_floor) != 1L || is.na(el_floor) ||
      el_floor < 0 || el_floor > 1)
    stop("`el_floor` must be a single number in [0, 1].")
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

  ## The LD matrix is turned into an edge list a BLOCK OF COLUMNS AT A TIME.
  ## Melting it whole (the previous approach) cost several full-size copies --
  ## `LD^2` duplicated the matrix, reshape2::melt built every cell as three
  ## columns including the ones about to be discarded, and each subsequent filter
  ## copied the table again -- which on a 15M-cell chromosome dominated peak
  ## memory. Working in blocks keeps the transient cost proportional to
  ## `cells_per_block` instead of to the whole matrix, and r2 is squared only for
  ## the pairs that survive filtering. Output is unchanged: same columns, same
  ## rows, same (column-major) order.
  cells_per_block <- 2e6

  compute_one <- function(local_idx) {
    ld <- SNPRelate::snpgdsLDMat(gds, snp.id = ids$snp_id[local_idx],
                                 method = method, slide = slide, verbose = FALSE,
                                 num.thread = as.integer(cores))$LD
    chr <- ids$snp_chr[local_idx]
    pos <- ids$snp_pos[local_idx]
    sid <- ids$snp_id[local_idx]
    n_snp <- length(sid)
    n_row <- nrow(ld); n_col <- ncol(ld)

    block <- max(1L, as.integer(cells_per_block %/% max(n_row, 1L)))
    parts <- vector("list", ceiling(n_col / block))
    k <- 0L

    for (start in seq.int(1L, n_col, by = block)) {
      cols <- start:min(start + block - 1L, n_col)
      v <- as.vector(ld[, cols, drop = FALSE])
      i <- rep.int(seq_len(n_row), length(cols))     # row index within the matrix
      j <- rep(cols, each = n_row)                   # column index

      ## slide mode: LD[i, j] is the pair (j, j + i); full mode: the pair (i, j),
      ## of which only the lower triangle is kept. SNP1 is the higher-indexed SNP
      ## in both cases, as before.
      a <- if (slide_win_ld > 0) i + j else i
      b <- j

      ## The floor is applied HERE, before the block is materialised as a
      ## data.table, so discarded pairs are never allocated. v is the LD value on
      ## the correlation scale, and r2 = v^2 is what the column reports.
      keep <- is.finite(v) & a > b & a <= n_snp
      if (el_floor > 0) keep <- keep & (v * v >= el_floor)
      if (!any(keep)) next
      a <- a[keep]; b <- b[keep]
      k <- k + 1L
      parts[[k]] <- data.table::data.table(
        SNP1 = sid[a], SNP2 = sid[b], Chr1 = chr[a], Chr2 = chr[b],
        pos1 = pos[a], pos2 = pos[b],
        r2   = v[keep]^2,                            # squared only for survivors
        d    = abs(pos[a] - pos[b])
      )
    }

    if (k == 0L) {
      return(data.table::data.table(
        SNP1 = sid[0], SNP2 = sid[0], Chr1 = chr[0], Chr2 = chr[0],
        pos1 = pos[0], pos2 = pos[0], r2 = numeric(0), d = pos[0]))
    }
    data.table::rbindlist(parts[seq_len(k)], use.names = FALSE)
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
