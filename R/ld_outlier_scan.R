#' LD-aware outlier regions from stage-1 clusters, end to end
#'
#' The sibling entry point to [ld_scan()], composing [ld_outlier_test()],
#' [ld_outlier_perm()] and [ld_region_rotation()] the same way [ld_scan()]
#' composes its own C-score sub-functions -- for one-call convenience, while
#' each piece stays independently usable for anyone who wants to reuse region
#' assembly on a different region set, or run the rotation null on regions from
#' elsewhere.
#'
#' This is a genuinely different method from [ld_scan()], not an alternative
#' parameterisation of it: that one computes a consistency C-score (Fang et al.
#' 2021) over an `ld_w` quantile grid and thresholds it; this one tests
#' [ld_complexity_reduction()] clusters directly (Simes over members, or your
#' own test on an [ld_unit_matrix()] variable), BH-corrects, and assembles
#' significant units into regions. Both turn p-values plus LD structure into
#' significant regions; pick one rather than have it picked for you.
#'
#' @inheritParams ld_outlier_test
#' @param p_perm Optional surrogate p-values (same shape `p_obs` had) --
#'   supplying this runs [ld_outlier_perm()] and adds `$null` to the result.
#'   `NULL` (default) skips it.
#' @param annotation Optional bed-like data.table -- supplying this (with
#'   `chrom_lengths`) runs [ld_region_rotation()] on the assembled regions and
#'   adds `$rotation` to the result. `NULL` (default) skips it.
#' @param chrom_lengths Required if `annotation` is supplied.
#' @param n_rotations,rotation_scheme Passed to [ld_region_rotation()].
#' @param B,perm_level,cores Passed to [ld_outlier_perm()].
#' @param verbose Print progress (default `TRUE`).
#'
#' @return An `ld_outlier_scan` object: `test` (the `ld_outlier_test` result --
#'   `units`, `regions`), `null` (the `ld_outlier_perm` result, or `NULL`),
#'   `rotation` (the `ld_region_rotation` result, or `NULL`), `params`.
#'
#' @seealso [ld_scan()], [ld_unit_matrix()], [ld_outlier_test()],
#'   [ld_outlier_perm()], [ld_region_rotation()]
#' @export
ld_outlier_scan <- function(stage1, map, p_obs, p_perm = NULL,
                            statistic = c("simes", "unit"),
                            size_floor = 8L, alpha = 0.05,
                            assembly = c("stage2_discovered", "physical"),
                            GTs = NULL, LD_decay = NULL,
                            score_threshold = 0.80, distance_threshold = 1e5,
                            gap = 3e5,
                            annotation = NULL, chrom_lengths = NULL,
                            n_rotations = 10000L,
                            rotation_scheme = c("within", "genome"),
                            B = NULL, perm_level = c("units", "regions"),
                            cores = 1L, verbose = TRUE) {
  statistic <- match.arg(statistic)
  assembly <- match.arg(assembly)
  rotation_scheme <- match.arg(rotation_scheme)
  perm_level <- match.arg(perm_level)
  if (!is.null(annotation) && is.null(chrom_lengths))
    stop("`chrom_lengths` is required when `annotation` is supplied.")

  test <- ld_outlier_test(stage1, map, p_obs, statistic = statistic,
                          size_floor = size_floor, alpha = alpha, assembly = assembly,
                          GTs = GTs, LD_decay = LD_decay,
                          score_threshold = score_threshold,
                          distance_threshold = distance_threshold, gap = gap)
  if (verbose) print(test)

  null <- if (is.null(p_perm)) NULL else
    ld_outlier_perm(test, stage1, map, p_perm, GTs = GTs, LD_decay = LD_decay,
                    B = B, level = perm_level, cores = cores, verbose = verbose)
  if (verbose && !is.null(null)) print(null)

  rotation <- if (is.null(annotation)) NULL else
    ld_region_rotation(test$regions, annotation, chrom_lengths, scheme = rotation_scheme,
                       n_rotations = n_rotations, seed = 1L)
  if (verbose && !is.null(rotation)) print(rotation)

  structure(list(test = test, null = null, rotation = rotation,
                 params = c(test$params, list(n_rotations = n_rotations,
                                               rotation_scheme = rotation_scheme))),
            class = "ld_outlier_scan")
}

#' @export
print.ld_outlier_scan <- function(x, ...) {
  print(x$test)
  if (!is.null(x$null)) print(x$null)
  if (!is.null(x$rotation)) print(x$rotation)
  invisible(x)
}
