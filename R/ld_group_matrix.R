#' One consensus-optimal SNP per `ld_prune_and_eMLG()` group -- a different
#' scale from `ld_unit_matrix()`
#'
#' [eMLG_best_snp()] picks, per group, the member SNP most correlated with
#' the group's eMLG consensus, and fills its missing calls from that
#' consensus. This function turns that into a ready-to-use genotype matrix,
#' one column per eligible group.
#'
#' **This is not a fourth `ld_unit_matrix()` representation, and its output
#' cannot be handed to `ld_outlier_test(statistic = "unit")` as-is.**
#' [ld_prune_and_eMLG()]'s groups are not guaranteed to be the same
#' partition as `stage1`'s Stage-1 units: an unflagged Stage-1 cluster maps
#' 1:1 onto one group, but a flagged cluster can be MERGED with
#' physically-nearby, correlated neighbours into one group spanning several
#' Stage-1 clusters (see [ld_prune_and_eMLG()]'s distance-restricted dynamic
#' cut). `ld_outlier_test(statistic = "unit")` expects one p-value per
#' Stage-1 unit, in the order [ld_unit_matrix()] would return them at the
#' same `size_floor` -- a different count, order and identity in general
#' from this function's columns. Testing this matrix's columns and feeding
#' the result to `ld_outlier_test()` will therefore either error (a length
#' mismatch, if you are fortunate) or silently score a p-value against the
#' wrong unit. If you test these columns, you are testing at the
#' consolidated (Stage-2) scale `ld_prune_and_eMLG()` produced, and any
#' significance/multiple-testing workflow built on the result is your own,
#' not `ld_outlier_test()`'s.
#'
#' @param GTs Genotype dosage matrix (individuals x SNPs, column names = markers).
#' @param prune_result An [ld_prune_and_eMLG()] result (must contain a
#'   non-empty `eMLG` matrix and its `groups` table).
#' @param map data.frame/data.table with `marker`, `Chr`, `Pos`, covering
#'   every marker `prune_result`'s groups reference.
#' @param size_floor Minimum markers per group to be included (default
#'   8L), applied to `prune_result$groups$n_loci` -- independent of
#'   whether the Stage-1 clusters that fed `prune_result` themselves
#'   cleared any size floor.
#' @param best_snp_args Extra arguments to [eMLG_best_snp()]. Must not set
#'   `fill = FALSE` -- this function needs the filled genotype matrix, not
#'   just [eMLG_best_snp()]'s stats table.
#'
#' @return An individuals x groups matrix (each column is an observed SNP
#'   with gaps filled from its group's consensus). Carries
#'   `attr(, "units")`: a data.table with `unit_id` (`prune_result`'s
#'   `group_id`, a character id, e.g. `"U12"` for an unmerged Stage-1
#'   cluster or `"F3"` for a merged run), `Chr`, `from`, `to`, `n_markers`
#'   (`= n_loci`) and `best_marker` (the underlying SNP each column
#'   actually is), aligned to columns.
#'
#' @seealso [eMLG_best_snp()], [ld_prune_and_eMLG()], [ld_unit_matrix()],
#'   [ld_outlier_test()]
#' @export
ld_group_matrix <- function(GTs, prune_result, map, size_floor = 8L, best_snp_args = list()) {
  args <- utils::modifyList(list(result = prune_result, GTs = GTs), best_snp_args)
  res <- do.call(eMLG_best_snp, args)
  if (is.null(res$geno))
    stop("ld_group_matrix() needs a genotype matrix -- pass `fill = TRUE` (the ",
         "default) via `best_snp_args`, not `fill = FALSE`.")

  groups <- data.table::as.data.table(prune_result$groups)
  eligible <- groups[has_eMLG == TRUE & n_loci >= size_floor]
  if (!nrow(eligible))
    stop("No `prune_result` group with a stored eMLG clears size_floor = ", size_floor, ".")
  keep_id <- colnames(res$geno)[colnames(res$geno) %in% eligible$group_id]
  m <- res$geno[, keep_id, drop = FALSE]

  g <- eligible[match(keep_id, group_id)]
  mp <- data.table::as.data.table(map)
  span <- data.table::rbindlist(lapply(g$members, function(mk) {
    p <- mp[marker %chin% mk, Pos]
    data.table::data.table(from = min(p), to = max(p))
  }))
  best <- res$stats[match(keep_id, group_id), best_marker]
  attr(m, "units") <- data.table::data.table(
    unit_id = g$group_id, Chr = g$Chr, from = span$from, to = span$to,
    n_markers = g$n_loci, best_marker = best
  )
  m
}
