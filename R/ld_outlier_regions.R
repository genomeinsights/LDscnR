#' Call LD-aware outlier regions with a null-calibrated threshold
#'
#' The end-to-end LD-scan pipeline. From a phenotype, genotypes and a kinship it
#' (1) runs a fast EMMAX scan, (2) computes the per-SNP consistency C-score
#' ([ld_cscore()]) integrating the `ld_w` window and stringency, (3) draws a
#' structure-aware null ([structured_null()]) and calibrates the C-score threshold
#' `tau_C` to a target FDR at the region unit ([calibrate_tauc()]), and (4)
#' clusters the SNPs above `tau_C` into outlier regions ([ld_regions()]) keeping
#' those with at least `l_min` SNPs.
#'
#' Observed statistics and the null share one engine (fast EMMAX with `K`), which
#' is essential for a valid calibration: the null must be able to reproduce the
#' observed structure-driven false positives. For a method without a fast null
#' (e.g. LFMM), compute its C-score separately and transfer this run's `tau_C`
#' with [gc_map_tauc()].
#'
#' The two swept parameters (`rho`, `q*`) are integrated into the C-score; `alpha`
#' is fixed; clustering is LD-decay-informed (`rho_ld`) with a hard distance cap
#' (`dcap`). The two operational knobs are `tau_C` (calibrated here, not set) and
#' `l_min` (region size; the region-level null makes a larger `l_min` a strong,
#' clean false-positive control).
#'
#' @param y Phenotype vector (length = individuals).
#' @param GTs Genotype dosage matrix (individuals x SNPs; column names = markers).
#' @param K Kinship / GRM (see [structured_null()] on the recommended GRM).
#' @param ld_ws Local-LD support matrix (SNPs x `ld_w` windows; see [ld_cscore()]).
#' @param map data.frame/data.table with `marker`, `Chr`, `Pos` (aligned to `GTs`).
#' @param decay_sum Per-chromosome LD-decay summary (see [ld_edges()]).
#' @param alpha Fixed within-candidate FDR for the C-score (default 0.05).
#' @param rho_ld,dcap,rho_d Clustering parameters (see [ld_edges()]; defaults 0.9,
#'   500 kb hard cap, `rho_d = NULL`).
#' @param l_min Minimum region size in SNPs (default 2), or `"auto"` to set it
#'   data-drivenly from the structure null via [calibrate_lmin()].
#' @param lmin_q,lmin_tau Passed to [calibrate_lmin()] when `l_min = "auto"`: the
#'   upper quantile of the null max-region-size distribution (default 0.99) and the
#'   reference C-threshold at which it is read (default 0.05).
#' @param fdr Target region-level FDR for `tau_C` (default 0.05).
#' @param basis,coords Structure basis for the null (see [structured_null()]).
#' @param B Number of null surrogates (default 100).
#' @param rho,qstar C-score integration grids (see [ld_cscore()]).
#' @param tau_grid `tau_C` grid searched for the FDR target.
#' @param keep_cache If `TRUE`, also return the null bundle and edge cache for
#'   ad-hoc re-estimation at other `tau_C` / `l_min`.
#' @param seed RNG seed for the null.
#'
#' @return An `ld_outlier_regions` object: `regions` (list of member-marker
#'   vectors), `summary` (per-region `Chr`, `start`, `end`, `size`, `maxC`),
#'   `tau_C`, `C` (full per-SNP C-score), `fdr_curve`, `params`, and -- if
#'   `keep_cache` -- `null` and `edges`.
#' @seealso [ld_cscore()], [structured_null()], [calibrate_tauc()], [ld_regions()],
#'   [gc_map_tauc()], [evaluate_ors()]
#' @export
ld_outlier_regions <- function(y, GTs, K, ld_ws, map, decay_sum,
                               alpha = 0.05, rho_ld = 0.9, dcap = 5e5, rho_d = NULL,
                               l_min = 2L, lmin_q = 0.99, lmin_tau = 0.05, fdr = 0.05,
                               basis = c("genetic", "spatial"), coords = NULL,
                               B = 100L, rho = colnames(ld_ws), qstar = seq(0, 0.95, by = 0.05),
                               tau_grid = seq(0.05, 1, by = 0.05), keep_cache = FALSE, seed = 1L) {
  basis <- match.arg(basis)
  prep <- emmax_setup(GTs, K)
  null <- structured_null(y, GTs, K, ld_ws, basis = basis, coords = coords, B = B,
                          alpha = alpha, rho = rho, qstar = qstar, prep = prep, seed = seed)
  C <- null$C_obs
  edges <- ld_edges(null$universe, GTs, map, decay_sum, rho_ld = rho_ld, dcap = dcap, rho_d = rho_d)
  if (identical(l_min, "auto"))                       # data-driven region-size floor from the null
    l_min <- calibrate_lmin(null, edges, tau = lmin_tau, q = lmin_q)
  fdr_curve <- null_fdr(null, edges, tau_grid, l_min)
  tau_C <- calibrate_tauc(null, edges, l_min, fdr, tau_grid)
  regions <- if (!is.na(tau_C)) {
    reg <- ld_regions(names(C)[C >= tau_C], edges); reg[lengths(reg) >= l_min]
  } else list()

  mp <- data.table::as.data.table(map)
  POS <- stats::setNames(mp$Pos, mp$marker); CH <- stats::setNames(as.character(mp$Chr), mp$marker)
  summary <- if (length(regions)) data.table::rbindlist(lapply(seq_along(regions), function(i) {
    m <- regions[[i]]
    data.table::data.table(region = i, Chr = names(sort(table(CH[m]), decreasing = TRUE))[1],
      start = min(POS[m]), end = max(POS[m]), size = length(m), maxC = max(C[m])) }))
    else data.table::data.table(region = integer(), Chr = character(), start = numeric(),
      end = numeric(), size = integer(), maxC = numeric())

  out <- list(regions = regions, summary = summary, tau_C = tau_C, C = C, fdr_curve = fdr_curve,
              params = list(alpha = alpha, rho_ld = rho_ld, dcap = dcap, l_min = l_min,
                            fdr = fdr, basis = basis, B = B))
  if (keep_cache) { out$null <- null; out$edges <- edges }
  structure(out, class = "ld_outlier_regions")
}

#' @export
print.ld_outlier_regions <- function(x, ...) {
  cat(sprintf("<ld_outlier_regions> %d region(s) at tau_C = %s (l_min = %d, FDR <= %.2f, %d null surrogates)\n",
              length(x$regions), format(x$tau_C), x$params$l_min, x$params$fdr, x$params$B))
  if (nrow(x$summary)) print(x$summary[order(-x$summary$size)][seq_len(min(nrow(x$summary), 10L))])
  invisible(x)
}
