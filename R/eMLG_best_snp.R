#' Best single representative SNP per eMLG block
#'
#' Post-processing step for [ld_prune_and_eMLG()]. For every block that has a
#' stored eMLG consensus (`has_eMLG = TRUE`, i.e. clusters larger than
#' `min_n_loci_eMLG`), this finds the member SNP whose genotype is most
#' correlated with the block consensus, and (optionally) returns that SNP's
#' genotype with its missing calls filled from the consensus.
#'
#' The eMLG consensus averages a block's markers, which is ideal when you want
#' one value per LD block but can dilute signal that is genuinely SNP-specific
#' (differing even between strongly-linked markers). In that situation a single
#' representative SNP is preferable. The clustering already returns a
#' representative marker per block (`result$groups$representative`), but it is
#' chosen for cluster centrality (highest median \eqn{r^2} to the rest of the
#' block), not to reproduce the consensus. This function instead selects, per
#' block, the member that best matches the consensus, reports how much better it
#' does than the centrality representative, and can hand back a ready-to-use
#' genotype: the observed SNP calls (preserving SNP-specific signal) with only
#' the missing calls repaired from the consensus.
#'
#' Correlation is evaluated on the absolute value (`|r|`): a member SNP can be in
#' repulsion with the arbitrarily-polarized consensus, and its coding would
#' simply be flipped to serve as a proxy. When filling, the consensus is oriented
#' to the chosen SNP first (`2 - consensus` when their correlation is negative);
#' such blocks are flagged in the `flipped` column.
#'
#' @param result A list returned by [ld_prune_and_eMLG()]. Must contain a
#'   non-empty `eMLG` matrix (built with `compute_unflagged_eMLG = TRUE` and/or
#'   flagged merges) and its `groups` table.
#' @param GTs The individuals-by-markers genotype matrix (0/1/2 dosage, `NA` for
#'   missing) that was passed to [ld_prune_and_eMLG()]. Its rows are matched to
#'   `result$eMLG` by row name, so both need individual IDs as row names.
#' @param fill Logical; if `TRUE` (default) also return a genotype matrix of the
#'   best SNP per block with missing calls filled from the (oriented) consensus.
#'   Observed SNP calls are never overwritten.
#' @param round_fill Logical; if `TRUE` (default) filled values are rounded to the
#'   nearest 0/1/2 so the returned matrix is an integer genotype. If `FALSE` the
#'   fractional consensus dosage is inserted (a numeric matrix).
#'
#' @return A list with:
#' \describe{
#'   \item{stats}{A `data.table`, one row per eMLG block, with `group_id`,
#'     `representative` (centrality choice), `best_marker` (consensus-optimal
#'     member), `n_loci`, `score` (eMLG fidelity), `best_r` (signed) and
#'     `best_abs_r`, `rep_abs_r` (the representative's `|r|` to the consensus),
#'     `rep_is_best`, `flipped`, and -- when `fill = TRUE` -- `n_obs` (observed
#'     SNP calls), `n_filled` (calls repaired from consensus) and `n_resid_na`
#'     (still missing because the consensus was also `NA`).}
#'   \item{geno}{When `fill = TRUE`, an individuals-by-blocks matrix of the
#'     best SNP per block, consensus-filled. Columns are named by `group_id`
#'     (aligning with `result$eMLG`); the underlying SNP is in
#'     `stats$best_marker`. `NULL` when `fill = FALSE`.}
#' }
#'
#' Blocks without a stored consensus (`has_eMLG = FALSE`, i.e. 1--2-marker
#' clusters) are already represented by a single SNP in `result$pruned` and are
#' not included here.
#'
#' @seealso [ld_prune_and_eMLG()], [ld_complexity_reduction()]
#'
#' @examples
#' \dontrun{
#' res  <- ld_prune_and_eMLG(GTs, stage1, LD_decay = ld_decay, rho = 0.95,
#'                           compute_unflagged_eMLG = TRUE, min_n_loci_eMLG = 3)
#' best <- eMLG_best_snp(res, GTs)
#' head(best$stats)
#' # drop-in single-SNP alternative to the consensus matrix res$eMLG:
#' dim(best$geno)
#' }
#' @export
eMLG_best_snp <- function(result, GTs, fill = TRUE, round_fill = TRUE) {
  if (is.null(result$eMLG) || !NCOL(result$eMLG))
    stop("`result$eMLG` is empty; run ld_prune_and_eMLG() with eMLG generation enabled.")
  if (is.null(result$groups))
    stop("`result$groups` is missing.")

  emlg   <- result$eMLG
  groups <- data.table::as.data.table(result$groups)

  ## Align individuals: by row name when both are named, otherwise positionally
  ## (result$eMLG rows are built directly from GTs rows, in the same order).
  if (!is.null(rownames(emlg)) && !is.null(rownames(GTs))) {
    if (!all(rownames(emlg) %in% rownames(GTs)))
      stop("Some individuals in `result$eMLG` are absent from `GTs`.")
    G <- GTs[rownames(emlg), , drop = FALSE]
  } else {
    if (nrow(GTs) != nrow(emlg))
      stop("`GTs` and `result$eMLG` differ in number of individuals and have no ",
           "row names to align by.")
    G <- GTs
  }

  gi <- groups[match(colnames(emlg), groups$group_id)]
  if (anyNA(gi$group_id))
    stop("Some `result$eMLG` columns have no matching row in `result$groups`.")
  members <- gi$members

  n_ind <- nrow(emlg); n_blk <- ncol(emlg)
  best_marker <- character(n_blk)
  best_r  <- rep(NA_real_, n_blk)
  rep_r   <- rep(NA_real_, n_blk)
  flipped <- logical(n_blk)
  n_obs <- n_filled <- n_resid_na <- integer(n_blk)
  geno <- if (fill) matrix(NA_real_, n_ind, n_blk,
                           dimnames = list(rownames(emlg), colnames(emlg))) else NULL

  for (j in seq_len(n_blk)) {
    cons <- emlg[, j]
    mk   <- members[[j]]
    mk   <- mk[mk %in% colnames(G)]

    if (!length(mk)) {                       # no usable member -> keep representative
      best_marker[j] <- gi$representative[j]
    } else {
      rr <- suppressWarnings(stats::cor(G[, mk, drop = FALSE], cons,
                                        use = "pairwise.complete.obs"))[, 1]
      if (all(is.na(rr))) {                  # correlation undefined -> keep representative
        best_marker[j] <- gi$representative[j]
      } else {
        bi <- which.max(abs(rr))
        best_marker[j] <- mk[bi]
        best_r[j]      <- rr[bi]
        flipped[j]     <- isTRUE(rr[bi] < 0)
      }
      if (gi$representative[j] %in% mk)
        rep_r[j] <- rr[match(gi$representative[j], mk)]
    }

    if (fill) {
      snp  <- G[, best_marker[j]]
      co   <- if (flipped[j]) 2 - cons else cons
      miss <- is.na(snp)
      out  <- snp
      fv   <- co[miss]
      if (round_fill) fv <- pmin(pmax(round(fv), 0), 2)
      out[miss] <- fv
      geno[, j]     <- out
      n_obs[j]      <- sum(!miss)
      n_filled[j]   <- sum(miss & !is.na(co))
      n_resid_na[j] <- sum(is.na(out))
    }
  }
  if (fill && round_fill) storage.mode(geno) <- "integer"

  stats <- data.table::data.table(
    group_id       = colnames(emlg),
    representative = gi$representative,
    best_marker    = best_marker,
    n_loci         = gi$n_loci,
    score          = gi$score,
    best_r         = best_r,
    best_abs_r     = abs(best_r),
    rep_abs_r      = abs(rep_r),
    rep_is_best    = best_marker == gi$representative,
    flipped        = flipped
  )
  if (fill)
    stats[, c("n_obs", "n_filled", "n_resid_na") := list(n_obs, n_filled, n_resid_na)]

  list(stats = stats[], geno = geno)
}
