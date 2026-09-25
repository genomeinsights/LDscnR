## Builds a zero-argument closure that returns one span-preserving,
## bounds-respecting relocation of `R` (a data.table with `chr`, `from`, `to`)
## for the given scheme. `len_of` is a named numeric vector, chromosome ->
## length. Not exported; factored out of ld_region_rotation() so its
## placement logic can be unit-tested directly (via `LDscnR:::.build_relocator`)
## on the raw relocated coordinates, not only on downstream overlap counts.
.build_relocator <- function(R, len_of, scheme) {
  n_reg <- nrow(R)
  w <- R$to - R$from

  if (n_reg == 0L) {
    ## Guard explicitly rather than relying on apply()/max.col() on a
    ## zero-row matrix in the "genome" branch below, which does not degrade
    ## to a zero-length result the way the rest of this function's vectorised
    ## operations do.
    return(function() R)
  }

  if (scheme == "within") {
    Lc <- len_of[R$chr]
    if (anyNA(Lc))
      stop("region chromosome(s) not found in `chrom_lengths`: ",
           paste(unique(R$chr[is.na(Lc)]), collapse = ", "))
    space <- Lc - w
    if (any(bad <- space < 0))
      stop(sprintf(
        "region %d spans %.0f bp but its chromosome (%s) is only %.0f bp long -- ",
        which(bad)[1], w[which(bad)[1]], R$chr[which(bad)[1]], Lc[which(bad)[1]]),
        "cannot place under scheme = \"within\".")
    function() {
      off <- stats::runif(n_reg, 0, space)
      data.table::copy(R)[, `:=`(from = off, to = off + w)]
    }
  } else {
    chr_names <- names(len_of)
    if (anyNA(len_of))
      stop("`chrom_lengths` contains NA length(s) for: ",
           paste(chr_names[is.na(len_of)], collapse = ", "))
    ## Placement space for region i on chromosome j is max(len_j - w_i, 0):
    ## the number of offsets at which the region's full span still fits.
    ## Sampling chromosomes with probability proportional to this (not
    ## uniformly) is what "reflects the available placement space" means --
    ## a region is more likely to land on a chromosome that actually offers
    ## more room for it, and chromosomes it cannot fit on get probability 0
    ## rather than truncating the region there.
    ##
    ## Selection WEIGHT is space + 1 on every chromosome the region actually
    ## fits on (len_j >= w_i), not space itself: space alone is exactly 0
    ## when a region's span exactly equals a chromosome's length -- the
    ## tightest possible fit, with exactly one valid placement (off = 0),
    ## not zero. Treating that as zero probability wrongly excluded the
    ## chromosome from selection, and could make row_tot = 0 for a region
    ## that DOES fit somewhere, raising the "does not fit anywhere" error on
    ## a region that fits everywhere it was checked -- reproduced directly:
    ## a single chromosome exactly the length of the region errored before
    ## this fix. The "+1" is the same add-one convention this package
    ## already uses for its p-values (never let a real possibility be
    ## weighted at exactly zero); for any chromosome with real headroom
    ## (space >> 1) it is negligible. The offset actually drawn still uses
    ## the real `space_mat` (0 at an exact fit, forcing off = 0), never the
    ## weight -- only chromosome SELECTION uses the smoothed weight.
    fits_mat <- outer(w, len_of, function(wi, Lj) Lj >= wi)
    space_mat <- outer(w, len_of, function(wi, Lj) pmax(Lj - wi, 0))
    weight_mat <- ifelse(fits_mat, space_mat + 1, 0)
    row_tot <- rowSums(weight_mat)
    if (any(bad <- row_tot <= 0))
      stop(sprintf(
        "region %d (span %.0f bp) does not fit on any chromosome in `chrom_lengths` ",
        which(bad)[1], w[which(bad)[1]]),
        "under scheme = \"genome\".")
    cum_mat <- t(apply(weight_mat, 1L, cumsum))
    function() {
      u <- stats::runif(n_reg) * row_tot
      col_idx <- max.col(cum_mat >= u, ties.method = "first")
      chosen_chr <- chr_names[col_idx]
      chosen_space <- space_mat[cbind(seq_len(n_reg), col_idx)]
      off <- stats::runif(n_reg, 0, chosen_space)
      data.table::copy(R)[, `:=`(chr = chosen_chr, from = off, to = off + w)]
    }
  }
}

#' Span-preserving random-relocation null: do these regions overlap an
#' annotation more than chance?
#'
#' General-purpose, not tied to any one outlier-detection method: `regions` can
#' come from [ld_outlier_test()], from [ld_scan()]/[ld_outlier_regions()]'s
#' C-score family, or from anywhere else -- this only needs a table of genomic
#' intervals.
#'
#' A raw overlap RATE is not interpretable on its own: a wider region overlaps a
#' fixed annotation more readily regardless of whether it is better localised,
#' so a set of wide, poorly-localised regions can show a higher raw rate than a
#' set of narrow, well-localised ones. Each null draw relocates every region to
#' an independently and uniformly chosen new position that still preserves its
#' OBSERVED SPAN and keeps it entirely on a valid chromosome, so that advantage
#' is present in the null exactly as in the observation and cannot inflate the
#' result -- read the fold and the null p-value, not the raw rate. This is a
#' random-relocation null (as in `bedtools shuffle` or GAT), not a circular
#' rotation: each region is placed independently, not shifted together as one
#' configuration, so it does not preserve the spacing between regions. That
#' independence matches how these regions are generated -- each is its own
#' discovery from an independently significant Stage-1/Stage-2 unit, not one
#' jointly-patterned point process whose internal geometry should be held
#' fixed.
#'
#' @param regions data.table with `Chr`/`chr_num`, `from`, `to` (one row per
#'   region).
#' @param annotation data.table with the same chromosome column and `start`,
#'   `end` (the intervals being tested against -- an EcoPeak set, a QTL panel,
#'   a gene list, anything).
#' @param chrom_lengths data.table with the same chromosome column and `len`.
#' @param scheme `"within"` (default) preserves each region's chromosome
#'   assignment and draws its new position uniformly from every offset at
#'   which its full span still fits on that chromosome -- the right default
#'   whenever the annotation is non-uniformly distributed among chromosomes,
#'   since `"genome"` would then credit a region merely for landing on an
#'   annotation-rich chromosome. `"genome"` also reassigns chromosome: for
#'   each region, a chromosome is drawn with probability proportional to its
#'   OWN valid placement space for that region's span, `max(chrom_len -
#'   span, 0)` -- not uniformly among chromosomes, which would over-place
#'   regions on short chromosomes relative to the room they actually offer,
#'   and under-place them on long ones. A region whose span exceeds every
#'   available chromosome (or, under `"within"`, its own chromosome) cannot be
#'   placed at all and raises an error identifying the offending region rather
#'   than truncating its span or silently dropping it.
#' @param n_rotations Relocation draws (default 10000L).
#' @param seed Seed for the relocation draws (default 1L).
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
  relocate_once <- .build_relocator(R, len_of, scheme)
  set.seed(seed)
  null <- vapply(seq_len(n_rotations), function(i) n_overlap(relocate_once()), 0L)

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
