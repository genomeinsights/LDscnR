
#' @export
create_gds_from_geno <- function(geno, map, gds_path) {
  stopifnot(ncol(geno) == nrow(map))
  snpgdsCreateGeno(
    gds_path,
    genmat         = t(round(geno)),     # SNPRelate expects SNP × sample if snpfirstdim=TRUE
    sample.id      = paste0("ind_", seq_len(nrow(geno))), ## sample name not important
    snp.id         = map$marker,
    snp.chromosome = map$Chr,
    snp.position   = map$Pos,
    snpfirstdim    = TRUE
  )
  snpgdsOpen(gds_path)
}


.read_gds_ids <- function(gds) {
  list(
    snp_id  = read.gdsn(index.gdsn(gds, "snp.id")),
    snp_chr = read.gdsn(index.gdsn(gds, "snp.chromosome")),
    snp_pos = read.gdsn(index.gdsn(gds, "snp.position"))
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
get_el <- function(gds, idx, slide_win_ld = 1000, n_cores = 1) {

  ids <- .read_gds_ids(gds)

  if(missing(idx)){
    idx <- seq_along(ids$snp_id)
  }else{
    idx <- as.integer(idx)
  }

  if (anyNA(idx) || length(idx) == 0L) stop("`idx` must be a non-empty integer vector.")
  if (min(idx) < 1L || max(idx) > length(ids$snp_id)) stop("`idx` is out of range for this GDS.")

  snp_ids <- ids$snp_id[idx]
  slide <- if (slide_win_ld > 0) as.integer(slide_win_ld) else -1L

  ldmat <- SNPRelate::snpgdsLDMat(
    gds,
    snp.id = snp_ids,
    method = "r",
    slide = slide,
    verbose = FALSE,
    num.thread = as.integer(n_cores)
  )

  el <- data.table::as.data.table(
    reshape2::melt(ldmat$LD^2, value.name = "r2")
  )

  ## Sliding-window LDMat uses a banded/offset representation: Var2 is an offset.
  ## Convert to absolute indices within the SNP subset.
  if (slide_win_ld > 0) {
    el[, Var1 := Var1 + Var2]
  }

  ## Keep only one copy per pair (drop diagonal + duplicates)
  ## Choose one triangle consistently; here: Var1 > Var2
  el <- el[Var1 > Var2]

  ## Drop missing/non-finite r2
  el <- el[is.finite(r2)]

  ## Attach metadata
  el[, Chr1 := ids$snp_chr[idx][Var1]]
  el[, Chr2 := ids$snp_chr[idx][Var2]]
  el[, pos1 := ids$snp_pos[idx][Var1]]
  el[, pos2 := ids$snp_pos[idx][Var2]]
  el[, SNP1 := ids$snp_id[idx][Var1]]
  el[, SNP2 := ids$snp_id[idx][Var2]]

  ## Physical distance
  el[, d := abs(pos1 - pos2)]

  el
}
