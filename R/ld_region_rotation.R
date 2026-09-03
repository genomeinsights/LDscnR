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
  R <- data.table::as.data.table(regions)
  A <- data.table::as.data.table(annotation)
  L <- data.table::as.data.table(chrom_lengths)
  ## Accept "chr", "Chr" or "chr_num" -- whichever the caller's table already uses.
  ## First bug found here: checking only for "Chr"/"chr_num" missed the (very common)
  ## case where a column is already named "chr", which then errored trying to rename
  ## "chr_num" onto a table that never had that name.
  .chr_col <- function(d) intersect(c("chr", "Chr", "chr_num"), names(d))[1]
  cR <- .chr_col(R); cA <- .chr_col(A); cL <- .chr_col(L)
  if (is.na(cR) || is.na(cA) || is.na(cL))
    stop("regions/annotation/chrom_lengths must each have a chr/Chr/chr_num column.")
  if (cR != "chr") data.table::setnames(R, cR, "chr")
  if (cA != "chr") data.table::setnames(A, cA, "chr")
  if (cL != "chr") data.table::setnames(L, cL, "chr")
  data.table::setkey(A, chr, start, end)

  ## number of R's rows overlapping A, by any amount
  n_overlap <- function(d) {
    if (!nrow(d)) return(0L)
    ov <- data.table::foverlaps(d[, .(chr, from, to)], A,
                                by.x = c("chr", "from", "to"), by.y = c("chr", "start", "end"),
                                type = "any", mult = "first", nomatch = NA)
    sum(!is.na(ov$start))
  }
  observed <- n_overlap(R)

  len_of <- stats::setNames(L$len, L$chr)
  rot_once <- function() {
    d <- data.table::copy(R)[, w := to - from]
    if (scheme == "within") {
      Lc <- len_of[d$chr]
      off <- stats::runif(nrow(d), 0, Lc)
      d[, `:=`(from = off, to = off + w)]
    } else {
      d[, chr := sample(names(len_of), .N, replace = TRUE)]
      Lc <- len_of[d$chr]
      off <- stats::runif(nrow(d), 0, Lc)
      d[, `:=`(from = off, to = pmin(off + w, Lc))]
    }
    n_overlap(d)
  }
  set.seed(seed)
  null <- vapply(seq_len(n_rotations), function(i) rot_once(), 0L)

  ov <- data.table::foverlaps(R[, .(chr, from, to)], A, by.x = c("chr", "from", "to"),
                              by.y = c("chr", "start", "end"), type = "any",
                              mult = "first", nomatch = NA)
  per_region <- data.table::data.table(chr = R$chr, from = R$from, to = R$to,
                                       on_peak = !is.na(ov$start))

  structure(list(
    observed = observed, null_mean = mean(null), fold = observed / max(mean(null), 1e-9),
    p = (1 + sum(null >= observed)) / (n_rotations + 1),
    per_region = per_region,
    params = list(scheme = scheme, n_rotations = n_rotations, seed = seed)
  ), class = "ld_region_rotation")
}

#' @export
print.ld_region_rotation <- function(x, ...) {
  cat(sprintf("<ld_region_rotation> observed %d/%d | null mean %.2f | fold %.2fx | p = %.4f\n",
              x$observed, nrow(x$per_region), x$null_mean, x$fold, x$p))
  invisible(x)
}
