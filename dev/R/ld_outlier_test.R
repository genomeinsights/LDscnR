#' Test LD-complexity-reduction clusters and assemble significant ones into regions
#'
#' The core detection step of the stage-1-cluster outlier-region method: one of
#' two genuinely different algorithms in this package for turning per-marker
#' p-values plus LD structure into significant regions (the other is the
#' consistency C-score family behind [ld_scan()]/[ld_outlier_regions()]). No
#' null here -- see [ld_outlier_perm()] for the permutation null and
#' [ld_region_rotation()] for annotation-overlap validation, both of which
#' reuse this function's result rather than duplicate its logic.
#'
#' `statistic = "simes"` aggregates marker-level `p_obs` (`min(n * p_(i) / i)`
#' per cluster, so the penalty scales with cluster size); `statistic = "unit"`
#' takes `p_obs` already computed per unit -- built with [ld_unit_matrix()] and
#' scanned with your own engine -- and does no aggregation, only BH and assembly.
#'
#' `min_r2_rho` for `assembly = "stage2_discovered"` is read from
#' `stage1$params$rho`, NOT a separate argument: the whole point of this
#' assembly rule is that it introduces no parameter beyond what built `stage1`,
#' and a caller-supplied value could silently disagree with that. `stage1` must
#' therefore come from an `ld_complexity_reduction()` new enough to record
#' `params` (added alongside this function); an older `stage1` object errors
#' rather than guessing.
#'
#' @param stage1 An `ld_complexity_reduction` object (must have `$params$rho`).
#' @param map data.frame/data.table with `marker`, `Chr`, `Pos`, aligned to
#'   `stage1`'s marker universe.
#' @param p_obs Per-marker p-values (if `statistic = "simes"`) or per-unit
#'   p-values in unit order (if `statistic = "unit"`).
#' @param statistic `"simes"` or `"unit"`.
#' @param size_floor Minimum markers per tested unit (default 8L).
#' @param alpha BH level (default 0.05).
#' @param assembly `"stage2_discovered"` (default) -- re-run
#'   [ld_prune_and_eMLG()] over the significant clusters only, so a region can
#'   only span discovered signal -- or `"physical"`, a plain merge of
#'   significant clusters within `gap` on a chromosome. The two are not
#'   equivalent in general; `"physical"` is retained as a near-free check
#'   (agreement between an arbitrary rule and a motivated one is evidence
#'   neither is an artefact of the other), not as an equally-preferred
#'   alternative.
#' @param GTs,LD_decay Required only when `assembly = "stage2_discovered"`.
#' @param score_threshold,distance_threshold [ld_prune_and_eMLG()] parameters
#'   for `"stage2_discovered"`. `ld_w_threshold` is fixed at 0 internally, not
#'   exposed: used only to merge already-significant clusters it is not a
#'   gate, only a speed filter, and 0 filters nothing.
#' @param gap Physical merge distance in bp, used only when
#'   `assembly = "physical"`.
#'
#' @return An `ld_outlier_test` object: `units` (data.table, one row per tested
#'   unit -- `unit_id`, `Chr`, `from`, `to`, `n_markers`, `p`, `q`,
#'   `significant`), `regions` (data.table, one row per assembled region --
#'   `Chr`, `from`, `to`, `n_units`, `n_markers`, `occupancy`), `params`
#'   (everything above, resolved).
#'
#' @seealso [ld_unit_matrix()], [ld_outlier_perm()], [ld_region_rotation()],
#'   [ld_scan()] (the C-score alternative), [ld_complexity_reduction()]
#' @export
ld_outlier_test <- function(stage1, map, p_obs,
                            statistic = c("simes", "unit"),
                            size_floor = 8L, alpha = 0.05,
                            assembly = c("stage2_discovered", "physical"),
                            GTs = NULL, LD_decay = NULL,
                            score_threshold = 0.80, distance_threshold = 1e5,
                            gap = 3e5) {
  statistic <- match.arg(statistic)
  assembly <- match.arg(assembly)
  if (assembly == "stage2_discovered" && is.null(stage1$params$rho))
    stop("stage1$params$rho is missing -- this stage1 object predates ld_complexity_reduction() ",
         "recording its own rho. Rebuild it, or pass assembly = \"physical\".")
  stop("ld_outlier_test(): not yet implemented (signature under review)")
}

#' @export
print.ld_outlier_test <- function(x, ...) {
  p <- x$params
  cat(sprintf("<ld_outlier_test> %s / %s | floor = %d, alpha = %.2f | %d/%d units -> %d regions\n",
              p$statistic, p$assembly, p$size_floor, p$alpha,
              sum(x$units$significant), nrow(x$units), nrow(x$regions)))
  invisible(x)
}
