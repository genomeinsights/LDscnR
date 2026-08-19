#' Build the C-score and LD-edge graph for a region scan
#'
#' The data-prep step for [ld_region_stability()] (and for ad-hoc clustering at
#' any `tau_C` / `l_min`): from a phenotype (or a precomputed association
#' p-vector) it returns the full per-SNP consistency C-score ([ld_cscore()]) and
#' the cached LD r^2-edge graph ([ld_edges()]). Unlike [ld_outlier_regions()] it
#' draws no structure-aware null and commits to no threshold -- it just
#' assembles the two objects a threshold-free region ranking needs.
#'
#' The edge graph is deliberately built over the **entire** lit-up universe
#' (`names(C)[C > 0]`), not a gated subset, so that a downward `tau_C` sweep
#' never references a marker missing from the graph.
#'
#' Two entry points:
#' * supply `y` (and `K` or `prep`) to run a fast EMMAX scan for the p-values;
#' * or supply `pvals` directly (e.g. LFMM p-values) to skip EMMAX -- then `K`
#'   is not needed, though `GTs` is still required for the LD graph.
#'
#' `pvals` / the EMMAX output must be row-aligned to `ld_ws` (and to the columns
#' of `GTs` and rows of `map`), exactly as elsewhere in the package.
#'
#' @param y Phenotype vector (length = individuals). Ignored if `pvals` is given.
#' @param GTs Genotype dosage matrix (individuals x SNPs; column names = markers).
#' @param K Kinship / GRM for EMMAX. Not needed when `pvals` or `prep` is supplied.
#' @param ld_ws Local-LD support matrix (SNPs x `ld_w` windows; see [ld_cscore()]).
#' @param map data.frame/data.table with `marker`, `Chr`, `Pos` (aligned to `GTs`).
#' @param decay_sum Per-chromosome LD-decay summary (see [ld_edges()]).
#' @param alpha,rho,qstar C-score parameters passed to [ld_cscore()].
#' @param rho_ld,dcap,rho_d Clustering parameters passed to [ld_edges()].
#' @param prep Optional `emmax_prep` from [emmax_setup()] to reuse.
#' @param pvals Optional precomputed per-SNP p-values (skips EMMAX).
#'
#' @return A list with `C` (full per-SNP C-score), `edges` (an [ld_edges()]
#'   object over the `C > 0` universe), and `prep` (the `emmax_prep`, or `NULL`
#'   when `pvals` was supplied). Feed `C` and `edges` to [ld_region_stability()].
#' @seealso [ld_region_stability()], [ld_cscore()], [ld_edges()],
#'   [ld_outlier_regions()]
#' @export
ld_cscore_scan <- function(y = NULL, GTs, K = NULL, ld_ws, map, decay_sum,
                           alpha = 0.05, rho = colnames(ld_ws),
                           qstar = seq(0, 0.95, by = 0.05),
                           rho_ld = 0.9, dcap = 5e5, rho_d = NULL,
                           prep = NULL, pvals = NULL) {
  if (is.null(pvals)) {
    if (is.null(prep)) {
      if (is.null(K)) stop("supply one of `pvals`, `prep`, or `K` (with `y`).")
      prep <- emmax_setup(GTs, K)
    }
    if (is.null(y)) stop("`y` is required unless `pvals` is supplied.")
    pvals <- emmax_fast(prep, y)
  }
  C <- ld_cscore(pvals, ld_ws, alpha = alpha, rho = rho, qstar = qstar)
  edges <- ld_edges(names(C)[C > 0], GTs, map, decay_sum,
                    rho_ld = rho_ld, dcap = dcap, rho_d = rho_d)
  list(C = C, edges = edges, prep = prep)
}

#' Rank outlier regions by cross-parameter stability (threshold-free)
#'
#' Ranks LD-aware outlier regions without committing to any single `tau_C` or
#' `l_min`, by asking *in how many* of a `tau_C` x `l_min` grid of cells each
#' region survives -- the C-score idea lifted from the SNP to the region level.
#' It sets no threshold; it only orders regions, so it sidesteps the arbitrary
#' cutoff a `maxC`-style filter would impose.
#'
#' For every `tau_C` on the grid the lit-up SNPs are re-clustered ([ld_regions()])
#' and each base region is mapped to the largest re-clustered region covering any
#' of its members (its `maxsz` at that `tau_C`); the region *survives* a cell
#' `(tau_C, l_min)` when `maxsz >= l_min`. Two complementary statistics fall out:
#' * **`stability`** = the fraction of the grid's cells in which the region
#'   survives (`n_cells / (|tau_grid| * |l_min_grid|)`) -- a threshold-free
#'   region-level consistency score; rank regions by it (or by the raw
#'   `n_cells`). Large, consistently-elevated blocks survive nearly every cell;
#'   single high-C spikes survive only the low-`l_min` cells.
#' * **`persist_tau`** = the largest `tau_C` on the grid at which the region
#'   still holds `>= base_lmin` consistent SNPs (a continuous strength).
#'
#' **Adaptive grid (default `tau_grid = "auto"`, `l_min_grid = "auto"`).** The
#' grid covers exactly where regions are produced rather than a fixed box: the
#' `tau_C` axis is the distinct observed C values (data-driven breakpoints, from
#' the most permissive up to `maxC` where the last region vanishes, thinned to
#' `max_tau`), and the `l_min` axis runs from 1 to the largest region at the most
#' permissive `tau_C` (capped at `lmin_cap`). Supply explicit vectors to override.
#'
#' **Null-aware value (supply `null`).** Given a structure-aware null bundle
#' ([structured_null()]), a cell `(tau_C, l_min)` is *clean* when its region-level
#' FDR -- mean surrogate regions divided by observed regions -- is at or below
#' `fdr`. (A strict "no surrogate region at all" rule is degenerate at `B >= 100`,
#' where the mean count is essentially never exactly zero; the FDR gate is its
#' robust form.) The region then scores **`stability_null`** = the number of cells
#' where it is found *and* the cell is clean, divided by *all* grid cells. Keeping
#' the denominator at all cells means that when the null lights up much of the
#' parameter space the value drops for every region -- automatically discounting
#' datasets (or regions) whose apparent consistency lives where the null is also
#' active. Ranking uses `stability_null` when a null is supplied, and the plain
#' grid-survival `stability` otherwise.
#'
#' @param x Either a named per-SNP C-score vector ([ld_cscore()] /
#'   [ld_cscore_scan()]), or an `ld_outlier_regions` object built with
#'   `keep_cache = TRUE` (its `C`, `edges`, `regions`, `l_min`, and `null` are used).
#' @param edges An [ld_edges()] object over the full `C > 0` universe (from
#'   [ld_cscore_scan()] or [ld_edges()]). Required for the vector method.
#' @param regions Optional list of base regions (member-marker vectors) to rank.
#'   If `NULL`, regions are clustered at `base_tau` and filtered to `>= base_lmin`.
#' @param null Optional [structured_null()] bundle (needs `C_surr` over the same
#'   `edges` universe); enables the null-aware `stability_null` value.
#' @param tau_grid,l_min_grid Grid axes, or `"auto"` (default) for the data-driven
#'   grid described above.
#' @param base_tau,base_lmin The `tau_C` / `l_min` defining the base region set
#'   when `regions` is `NULL`, and the `l_min` reference for `persist_tau`
#'   (defaults: `min(tau_grid)`, `2`).
#' @param fdr Cell is null-clean when its region-level FDR (mean surrogate regions
#'   / observed regions) is at or below this (default `0.05`).
#' @param max_tau,lmin_cap Caps on the auto grid's axis lengths (defaults 40, 20).
#' @param map Optional `marker`/`Chr`/`Pos` table to annotate each region with
#'   `Chr`, `start`, `end`.
#' @param ... Passed on from the generic to the default method.
#'
#' @return An `ld_region_stability` data.table, ordered most- to least-stable,
#'   with `region`, (`Chr`, `start`, `end` if `map`), `size`, `maxC`,
#'   `persist_tau`, `n_cells`, `stability`, `n_clean`, `stability_null`,
#'   `rank_stability`, `rank_persist`. The grid, cell total, and whether a null was
#'   used are stored in the `"grid"` attribute.
#' @seealso [ld_cscore_scan()], [ld_regions()], [ld_outlier_regions()], [structured_null()]
#' @export
ld_region_stability <- function(x, ...) UseMethod("ld_region_stability")

#' @rdname ld_region_stability
#' @export
ld_region_stability.ld_outlier_regions <- function(x, tau_grid = "auto",
                                                    l_min_grid = "auto", map = NULL, ...) {
  if (is.null(x$edges))
    stop("no edge cache on this result -- re-run ld_outlier_regions(..., keep_cache = TRUE).")
  ld_region_stability.default(x$C, edges = x$edges, regions = x$regions, null = x$null,
                              tau_grid = tau_grid, l_min_grid = l_min_grid,
                              base_lmin = x$params$l_min, map = map, ...)
}

#' @rdname ld_region_stability
#' @export
ld_region_stability.default <- function(x, edges, regions = NULL,
                                        null = NULL, tau_grid = "auto", l_min_grid = "auto",
                                        base_tau = NULL, base_lmin = 2L, fdr = 0.05,
                                        max_tau = 40L, lmin_cap = 20L, map = NULL, ...) {
  C <- x
  if (!is.numeric(C) || is.null(names(C)))
    stop("`x` must be a named numeric C-score vector (or an ld_outlier_regions object).")

  ## --- adaptive grid: cover exactly where regions are produced -----------
  ## tau axis = the distinct observed C values (data-driven breakpoints, from the
  ## most permissive up to maxC where the last region vanishes), thinned to max_tau.
  if (identical(tau_grid, "auto")) {
    pos <- sort(unique(C[C > 0])); if (!length(pos)) pos <- 0
    if (length(pos) > max_tau) pos <- as.numeric(stats::quantile(pos, seq(0, 1, length.out = max_tau), type = 1))
    tau_grid <- unique(pos)
  } else tau_grid <- sort(unique(tau_grid))
  if (is.null(base_tau)) base_tau <- min(tau_grid)
  if (is.null(regions)) {
    regions <- ld_regions(names(C)[C >= base_tau], edges)
    regions <- regions[lengths(regions) >= base_lmin]
  }
  ## l_min axis = 1 up to the largest region at the most permissive tau (capped).
  if (identical(l_min_grid, "auto")) {
    allsz <- lengths(ld_regions(names(C)[C >= base_tau], edges))
    l_min_grid <- seq_len(min(lmin_cap, max(1L, if (length(allsz)) max(allsz) else 1L)))
  } else l_min_grid <- sort(unique(as.integer(l_min_grid)))
  ntot <- length(tau_grid) * length(l_min_grid)
  use_null <- !is.null(null) && !is.null(null$C_surr)

  if (!length(regions)) {
    empty <- data.table::data.table(region = integer(), Chr = character(),
      start = numeric(), end = numeric(), size = integer(), maxC = numeric(),
      persist_tau = numeric(), n_cells = integer(), stability = numeric(),
      n_clean = integer(), stability_null = numeric(),
      rank_stability = integer(), rank_persist = integer())
    return(structure(empty, class = c("ld_region_stability", class(empty)),
                     grid = list(tau = tau_grid, l_min = l_min_grid, n_cells_total = ntot, null = use_null)))
  }

  ncells <- ncells_clean <- integer(length(regions))
  persist <- rep(NA_real_, length(regions))
  for (t in tau_grid) {                                  # ascending -> last hit = largest persist_tau
    cr <- ld_regions(names(C)[C >= t], edges)
    sz <- if (length(cr)) lengths(cr) else integer(0)
    cm <- if (length(cr)) stats::setNames(rep(sz, sz), unlist(cr)) else stats::setNames(integer(0), character(0))
    ## which l_min cells are null-CLEAN at this tau: region-level FDR (mean surrogate
    ## regions / observed regions) at or below `fdr`. Robust to the stochastic null,
    ## where "any surrogate region at all" is almost never exactly zero at B>=100.
    clean_L <- rep(TRUE, length(l_min_grid))
    if (use_null) {
      ssz <- lapply(null$C_surr, function(Csp) { mk <- names(Csp)[Csp >= t]
        if (length(mk)) lengths(ld_regions(mk, edges)) else integer(0) })
      clean_L <- vapply(seq_along(l_min_grid), function(j) { L <- l_min_grid[j]
        nreg <- mean(vapply(ssz, function(s) sum(s >= L), numeric(1)))
        oreg <- sum(sz >= L)
        nreg <= fdr * max(oreg, 1)
      }, logical(1))
    }
    for (i in seq_along(regions)) {
      ms <- cm[regions[[i]]]; ms <- ms[!is.na(ms)]
      mx <- if (length(ms)) max(ms) else 0L
      hit <- l_min_grid <= mx                            # region found at (t, L)
      ncells[i]       <- ncells[i]       + sum(hit)
      ncells_clean[i] <- ncells_clean[i] + sum(hit & clean_L)   # ... and the null is silent there
      if (mx >= base_lmin) persist[i] <- t
    }
  }

  POS <- CH <- NULL
  if (!is.null(map)) {
    mp <- data.table::as.data.table(map)
    POS <- stats::setNames(mp$Pos, mp$marker)
    CH  <- stats::setNames(as.character(mp$Chr), mp$marker)
  }
  out <- data.table::rbindlist(lapply(seq_along(regions), function(i) {
    m <- regions[[i]]
    data.table::data.table(
      region = i,
      Chr   = if (!is.null(CH))  names(sort(table(CH[m]), decreasing = TRUE))[1] else NA_character_,
      start = if (!is.null(POS)) min(POS[m]) else NA_real_,
      end   = if (!is.null(POS)) max(POS[m]) else NA_real_,
      size  = length(m), maxC = max(C[m]), persist_tau = persist[i],
      n_cells = ncells[i], stability = ncells[i] / ntot,
      n_clean = ncells_clean[i], stability_null = ncells_clean[i] / ntot)
  }))
  ## rank by the null-aware value when a null is supplied, else by plain grid-survival
  key <- if (use_null) "n_clean" else "n_cells"
  out[, "rank_stability" := data.table::frank(-out[[key]], ties.method = "min")]
  out[, "rank_persist"   := data.table::frank(-out$persist_tau, ties.method = "min")]
  data.table::setorderv(out, c(key, "persist_tau"), order = -1L)
  structure(out[], class = c("ld_region_stability", class(out)),
            grid = list(tau = tau_grid, l_min = l_min_grid, n_cells_total = ntot, null = use_null))
}

#' @export
print.ld_region_stability <- function(x, ...) {
  g <- attr(x, "grid")
  cat(sprintf("<ld_region_stability> %d region(s) over a %d x %d tau_C x l_min grid (%d cells)%s\n",
              nrow(x), length(g$tau), length(g$l_min), g$n_cells_total,
              if (isTRUE(g$null)) "; null-aware (ranked by stability_null)" else ""))
  if (nrow(x)) {
    cols <- intersect(c("region", "Chr", "start", "size", "maxC", "persist_tau",
                        "n_cells", "stability", "n_clean", "stability_null"), names(x))
    print(utils::head(as.data.frame(x)[, cols, drop = FALSE], 15L), row.names = FALSE)
  }
  invisible(x)
}
