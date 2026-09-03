#' Span-preserving rotation null: do these regions overlap an annotation more
#' than chance?
#'
#' General-purpose, not tied to any one outlier-detection method: `regions` can
#' come from [ld_outlier_test()], from [ld_scan()]/[ld_outlier_regions()]'s
#' C-score family, or from anywhere else -- this only needs a table of genomic
#' intervals.
#'
#' A raw overlap RATE is not interpretable on its own: a wider region overlaps a
#' fixed annotation more readily regardless of whether it is better localised,
#' so a set of wide, poorly-localised regions can show a higher raw rate than a
#' set of narrow, well-localised ones. The rotation preserves each region's
#' OBSERVED SPAN and only randomises its position, so that advantage is present
#' in the null exactly as in the observation and cannot inflate the result --
#' read the fold and the rotation p-value, not the raw rate.
#'
#' @param regions data.table with `Chr`/`chr_num`, `from`, `to` (one row per
#'   region).
#' @param annotation data.table with the same chromosome column and `start`,
#'   `end` (the intervals being tested against -- an EcoPeak set, a QTL panel,
#'   a gene list, anything).
#' @param chrom_lengths data.table with the same chromosome column and `len`.
#' @param scheme `"within"` (default) preserves each region's chromosome
#'   assignment and rotates only its position on that chromosome -- the right
#'   default whenever the annotation is non-uniformly distributed among
#'   chromosomes, since `"genome"` would then credit a region merely for
#'   landing on an annotation-rich chromosome. `"genome"` rotates across the
#'   whole concatenated genome, reassigning chromosome as well as position.
#' @param n_rotations Rotation draws (default 10000L).
#' @param seed Seed for the rotation draws (default 1L).
#'
#' @return An `ld_region_rotation` object: `observed` (regions overlapping the
#'   annotation), `null_mean`, `fold` (`observed / null_mean`), `p` (one-sided),
#'   `per_region` (data.table, one row per region: `on_peak` logical),
#'   `params`.
#'
#' @seealso [ld_outlier_test()], [ld_outlier_perm()]
#' @export
ld_region_rotation <- function(regions, annotation, chrom_lengths,
                               scheme = c("within", "genome"),
                               n_rotations = 10000L, seed = 1L) {
  scheme <- match.arg(scheme)
  stop("ld_region_rotation(): not yet implemented (signature under review)")
}

#' @export
print.ld_region_rotation <- function(x, ...) {
  cat(sprintf("<ld_region_rotation> observed %d/%d | null mean %.2f | fold %.2fx | p = %.4f\n",
              x$observed, nrow(x$per_region), x$null_mean, x$fold, x$p))
  invisible(x)
}
