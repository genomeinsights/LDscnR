#' Pre-compute the LD r^2-edge graph over a candidate marker set
#'
#' Builds the within-chromosome r^2 linkage graph among a set of candidate
#' markers so that clustering any gated subset later ([ld_regions()]) is a pure
#' graph lookup, with no further genotype access. Two markers are linked when
#' their r^2 exceeds a (by default LD-decay-informed) threshold `r2_link`.
#'
#' The r^2 threshold is decay-relative and self-calibrates per chromosome:
#' `r2_link = ld_from_rho(b, c, rho_ld)`. `rho_ld = 0.75` is the canonical value,
#' used by both the stickleback and the simulation analyses. Note the direction:
#' a HIGHER `rho_ld` gives a LOWER r^2 link and therefore MORE merging.
#'
#' Insensitivity is bounded, not general. On the stickleback panel 0.60 and 0.75
#' give the same 17 regions, the same `q_R`, and the same null-gate verdicts --
#' the whole difference is one marker joining one region -- even though the r^2
#' link moves from 0.405 to 0.273. It does not hold indefinitely: by `rho_ld = 0.9`
#' (r^2 ~ 0.14) the linkage is loose enough to overmerge, and a chromosome carrying
#' several distinct peaks collapses into one region. Sweep it if your decay fits
#' differ materially from these.
#' The physical extent of an
#' outlier region is bounded by a hard distance cap `dcap` (default 500 kb)
#' rather than by LD decay, because a region's footprint is set by
#' selection/recombination, not by background LD decay -- long-range r^2 links are
#' the rule inside inversions, so a decay-derived distance would fragment them.
#' Supply a scalar `r2_link` to override the decay-relative r^2, or `rho_d` to
#' derive the distance from LD decay (capped at `dcap`) instead of the hard cap.
#'
#' @param markers Character vector of candidate markers (the "universe" that any
#'   downstream `tau_C`/`l_min` gate is a subset of).
#' @param GTs Genotype dosage matrix (individuals x SNPs; column names = markers).
#' @param map data.frame/data.table with columns `marker`, `Chr`, `Pos`.
#' @param decay_sum Per-chromosome LD-decay summary with `Chr`, `b`, `c`,
#'   `a_pred` (see [compute_LD_decay()]).
#' @param rho_ld Decay-relative r^2 link threshold (default 0.9).
#' @param dcap Hard distance cap in bp for the gap-split (default 5e5).
#' @param r2_link Optional scalar r^2 overriding the decay-relative value.
#' @param rho_d Optional decay-relative distance; overrides the hard cap but is
#'   itself capped at `dcap`.
#'
#' @return An object of class `ld_edges`: a per-chromosome list of member
#'   markers, positions, the chromosome's r^2 / distance thresholds, and the r^2
#'   edge list (a two-column character matrix of linked marker pairs). Consumed
#'   by [ld_regions()].
#' @seealso [ld_regions()]
#' @export
ld_edges <- function(markers, GTs, map, decay_sum, rho_ld = 0.75, dcap = 5e5,
                     r2_link = NULL, rho_d = NULL) {
  map <- data.table::as.data.table(map)
  m <- map[map$marker %in% markers, ]
  ds <- data.table::as.data.table(decay_sum)
  ac <- if ("a_pred" %in% names(ds)) ds$a_pred else ds$a
  cc <- if ("c" %in% names(ds)) ds$c else rep(1, nrow(ds))
  r2_by <- if (!is.null(r2_link)) stats::setNames(rep(r2_link, nrow(ds)), ds$Chr)
           else stats::setNames(ld_from_rho(ds$b, cc, rho_ld), ds$Chr)
  d_by  <- if (is.null(rho_d)) stats::setNames(rep(dcap, nrow(ds)), ds$Chr)
           else stats::setNames(pmin(d_from_rho(ac, rho_d), dcap), ds$Chr)
  chrs <- unique(m$Chr)
  structure(stats::setNames(lapply(chrs, function(ch) {
    mk <- m$marker[m$Chr == ch]; pos <- m$Pos[m$Chr == ch]
    r2c <- if (!is.na(r2_by[ch])) r2_by[[ch]] else 0.5
    dc  <- if (!is.na(d_by[ch]))  d_by[[ch]]  else Inf
    ed <- NULL
    if (length(mk) > 1L) {
      R <- suppressWarnings(stats::cor(GTs[, mk], use = "pairwise.complete.obs")^2)
      R[!is.finite(R)] <- 0; R[lower.tri(R, diag = TRUE)] <- 0
      w <- which(R >= r2c, arr.ind = TRUE)
      if (nrow(w)) ed <- cbind(mk[w[, 1]], mk[w[, 2]])
    }
    list(marker = mk, pos = pos, r2c = r2c, dc = dc, edges = ed)
  }), chrs), class = "ld_edges")
}

#' Cluster gated markers into outlier regions
#'
#' Partitions a gated marker set into outlier regions by single-linkage on the
#' cached r^2-edge graph ([ld_edges()]) -- i.e. connected components of the
#' induced r^2 subgraph -- followed by a physical gap-split wherever consecutive
#' members are more than the chromosome's distance cap apart. Needs no genotype
#' access, so any `tau_C` / `l_min` can be scored ad hoc; reproduces the
#' clustering step of [ld_outlier_clusters()] exactly.
#'
#' @param markers Character vector of gated markers, a subset of the [ld_edges()]
#'   universe, e.g. `names(C)[C >= tau_C]` for a C-score `C` ([ld_cscore()]).
#' @param edges An `ld_edges` object from [ld_edges()].
#'
#' @return A list of character vectors, one per outlier region (its member
#'   markers). Apply the region-size filter with
#'   `regions[lengths(regions) >= l_min]`.
#' @seealso [ld_edges()], [ld_cscore()]
#' @export
ld_regions <- function(markers, edges) {
  markers <- unique(markers)
  out <- list()
  for (E in edges) {
    inmk <- E$marker %in% markers
    mkc <- E$marker[inmk]
    if (!length(mkc)) next
    posc <- E$pos[inmk]
    if (length(mkc) == 1L) { out <- c(out, list(mkc)); next }
    comp <- stats::setNames(seq_along(mkc), mkc)          # union-find on the induced subgraph
    if (!is.null(E$edges) && nrow(E$edges)) {
      keep <- E$edges[, 1] %in% mkc & E$edges[, 2] %in% mkc
      ee <- E$edges[keep, , drop = FALSE]
      for (r in seq_len(nrow(ee))) {
        ca <- comp[[ee[r, 1]]]; cb <- comp[[ee[r, 2]]]
        if (ca != cb) comp[comp == cb] <- ca
      }
    }
    cl <- .split_gap(as.integer(factor(comp)), posc, E$dc)
    out <- c(out, split(mkc, cl))
  }
  unname(out)
}

## internal: split a single-linkage cluster wherever consecutive members (by
## position) are more than `dcap` bp apart. Not exported.
.split_gap <- function(cl, pos, dcap) {
  if (!is.finite(dcap)) return(cl)
  out <- integer(length(cl)); nxt <- 0L
  for (g in unique(cl)) {
    gi <- which(cl == g); o <- order(pos[gi])
    run <- integer(length(gi)); run[o] <- cumsum(c(TRUE, diff(pos[gi][o]) > dcap))
    out[gi] <- nxt + run; nxt <- nxt + max(run)
  }
  out
}
