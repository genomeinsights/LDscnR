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

#' Build an LD edge list for a subset of SNPs (unidirectional)
#'
#' Computes pairwise LD (\eqn{r^2}) for a specified subset of SNPs in a GDS file
#' and returns an edge list with genomic coordinates and distances. Each SNP pair
#' appears once (unidirectional), which is sufficient for undirected clustering.
#'
#' @param gds An open GDS object.
#' @param idx Integer vector of SNP indices (1-based, referring to the GDS SNP index).
#' @param slide_win_ld Integer. If > 0, LD is estimated using a sliding window of
#'   this size; if \code{<= 0}, all pairwise LD values are computed.
#' @param n_cores Integer. Number of threads to use in LD computation.
#'
#' @return A \code{data.table} with columns \code{Var1}, \code{Var2}, \code{r2},
#'   \code{Chr1}, \code{Chr2}, \code{pos1}, \code{pos2}, \code{SNP1}, \code{SNP2}, \code{d}.
#'
#' @export
get_el <- function(gds,
                   idx = NULL,
                   slide_win_ld = 1000,
                   cores = 1,
                   by_chr = FALSE) {

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

    el_list[[i]] <- compute_one(local_idx)
  }

  data.table::rbindlist(el_list, use.names = TRUE)
}


.get_n_inds <- function(gds) {
  length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "sample.id")))
}
