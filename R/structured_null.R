#' Structure-aware null for tau_C calibration
#'
#' Draws `B` surrogate phenotypes that keep the observed phenotype's *structure*
#' (its spatial or genetic autocorrelation) but are *orthogonal* to it, runs a
#' fast EMMAX scan for each ([emmax_fast()]), and stores the resulting per-SNP
#' C-scores ([ld_cscore()]). Because a surrogate reproduces the structure-driven
#' false-positive rate while carrying none of the true signal, the empirical FDR
#' of the C-score against these surrogates ([null_fdr()], [calibrate_tauc()])
#' gives a truth-free `tau_C` that prices in structure confounding -- unlike a
#' naive phenotype permutation, which breaks the structure<->phenotype link and
#' is too permissive.
#'
#' Surrogates keep the structure but are made orthogonal to `y` by Gram-Schmidt
#' (`resid(lm(surrogate ~ y))`):
#' * `basis = "genetic"`: `surrogate ~ MVN(0, K)` -- same genetic-relatedness
#'   covariance as the observed structure (use when confounding is genetic, e.g.
#'   ecotype divergence).
#' * `basis = "spatial"`: `surrogate ~ MVN(0, Gaussian kernel over coords)` --
#'   same spatial autocorrelation (use when confounding is spatial, e.g. an
#'   environmental gradient); requires `coords`.
#'
#' @param y Observed phenotype vector (length = number of individuals).
#' @param GTs Genotype dosage matrix (individuals x SNPs; column names = markers).
#' @param K Kinship / GRM used both to correct EMMAX and as the genetic structure
#'   basis. For empirical data prefer a GRM built from LD-independent markers so it
#'   captures structure without absorbing the signal.
#' @param ld_ws Local-LD support matrix for the C-score (see [ld_cscore()]).
#' @param basis Structure basis, `"genetic"` (default) or `"spatial"`.
#' @param coords Two-column matrix/data.frame of coordinates (required for
#'   `basis = "spatial"`).
#' @param B Number of surrogates (default 100).
#' @param alpha,rho,qstar C-score parameters passed to [ld_cscore()].
#' @param seed RNG seed for reproducible surrogates.
#'
#' @return An `ld_null` object: `C_obs` (observed sparse C, only `C > 0`),
#'   `C_surr` (list of `B` sparse surrogate C-vectors), `universe` (every marker
#'   any surrogate lights up), and the settings. Build an [ld_edges()] over
#'   `universe`, then call [null_fdr()] / [calibrate_tauc()].
#' @seealso [null_fdr()], [calibrate_tauc()], [emmax_setup()], [ld_cscore()]
#' @export
structured_null <- function(y, GTs, K, ld_ws, basis = c("genetic", "spatial"),
                            coords = NULL, B = 100L, alpha = 0.05,
                            rho = colnames(ld_ws), qstar = seq(0, 0.95, by = 0.05),
                            seed = 1L) {
  basis <- match.arg(basis)
  n <- length(y)
  prep <- emmax_setup(GTs, K)
  if (basis == "genetic") {
    eK <- eigen(K, symmetric = TRUE)
  } else {
    if (is.null(coords)) stop("`coords` is required for basis = \"spatial\".")
    Dm <- as.matrix(stats::dist(coords)); l <- stats::median(Dm[lower.tri(Dm)])
    eK <- eigen(exp(-0.5 * (Dm / l)^2), symmetric = TRUE)
  }
  Lv <- pmax(eK$values, 0); Vk <- eK$vectors
  gen <- function() {
    s <- as.numeric(Vk %*% (sqrt(Lv) * stats::rnorm(n)))
    as.numeric(stats::resid(stats::lm(s ~ y)))
  }
  sparseC <- function(pv) { C <- ld_cscore(pv, ld_ws, alpha, rho, qstar); C[C > 0] }
  C_obs <- sparseC(emmax_fast(prep, y))
  set.seed(seed)
  C_surr <- vector("list", B)
  for (b in seq_len(B)) C_surr[[b]] <- sparseC(emmax_fast(prep, gen()))
  universe <- unique(c(names(C_obs), unlist(lapply(C_surr, names))))
  structure(list(C_obs = C_obs, C_surr = C_surr, universe = universe,
                 basis = basis, B = B, params = list(alpha = alpha, rho = rho, qstar = qstar)),
            class = "ld_null")
}

#' Region-level FDR of the C-score from a structure-aware null
#'
#' Turns an [structured_null()] bundle into an empirical false-discovery-rate
#' curve over `tau_C`, at the region unit: for each `tau_C`, cluster the gated
#' SNPs ([ld_regions()]) and count regions of at least `l_min` SNPs, both for the
#' observed C and (averaged) for the surrogates. `FDR = mean surrogate regions /
#' observed regions`. Because structure-only surrogates almost never produce large
#' clusters, a large `l_min` makes the null very clean -- the region-size filter
#' is itself a strong false-positive control.
#'
#' @param bundle An `ld_null` object from [structured_null()].
#' @param edges An [ld_edges()] object built over `bundle$universe`.
#' @param tau_grid C-score thresholds to evaluate (default `seq(0.05, 1, 0.05)`).
#' @param l_min Minimum region size in SNPs (default 2).
#'
#' @return A data.frame with `tau_C`, `n_obs`, `n_null`, `fdr`.
#' @seealso [structured_null()], [calibrate_tauc()], [ld_edges()]
#' @export
null_fdr <- function(bundle, edges, tau_grid = seq(0.05, 1, by = 0.05), l_min = 2L) {
  n_reg <- function(Csp, tau) { mk <- names(Csp)[Csp >= tau]
    if (!length(mk)) return(0L); sum(lengths(ld_regions(mk, edges)) >= l_min) }
  n_obs  <- vapply(tau_grid, function(t) n_reg(bundle$C_obs, t), numeric(1))
  n_null <- vapply(tau_grid, function(t) mean(vapply(bundle$C_surr, n_reg, numeric(1), tau = t)), numeric(1))
  data.frame(tau_C = tau_grid, n_obs = n_obs, n_null = n_null, fdr = pmin(1, n_null / pmax(n_obs, 1)))
}

#' Calibrate tau_C to a target FDR
#'
#' The smallest `tau_C` whose region-level empirical FDR ([null_fdr()]) is at or
#' below `fdr` (with observed regions present). This is the truth-free operating
#' threshold for the C-score.
#'
#' @param bundle An `ld_null` object from [structured_null()].
#' @param edges An [ld_edges()] object over `bundle$universe`.
#' @param l_min Minimum region size (default 2).
#' @param fdr Target FDR (default 0.05).
#' @param tau_grid C-score grid (default `seq(0.05, 1, 0.05)`).
#'
#' @return The calibrated `tau_C` (numeric scalar), or `NA` if none reaches the
#'   target.
#' @seealso [null_fdr()], [gc_map_tauc()]
#' @export
calibrate_tauc <- function(bundle, edges, l_min = 2L, fdr = 0.05, tau_grid = seq(0.05, 1, by = 0.05)) {
  f <- null_fdr(bundle, edges, tau_grid, l_min)
  ok <- which(f$fdr <= fdr & f$n_obs > 0)
  if (length(ok)) f$tau_C[min(ok)] else NA_real_
}

#' Map a calibrated tau_C onto another method's C-scale (genomic control)
#'
#' Transfers a `tau_C` calibrated on a reference method (e.g. EMMAX, which has a
#' fast structure-aware null) to a target method that has no cheap null (e.g.
#' LFMM), by quantile-matching their C-score distributions: find the upper-tail
#' position of `tau_C` in the reference C, then read off the target C at the same
#' quantile. This is genomic control on the C-score -- it operates the target
#' method at the reference method's FDR anchor, correcting for a shifted
#' (e.g. inflated) C distribution.
#'
#' @param tau_C Calibrated reference threshold ([calibrate_tauc()]).
#' @param C_ref Full per-SNP C-score of the reference method ([ld_cscore()]).
#' @param C_target Full per-SNP C-score of the target method.
#'
#' @return The target method's `tau_C` on its own C-scale (numeric scalar).
#' @seealso [calibrate_tauc()]
#' @export
gc_map_tauc <- function(tau_C, C_ref, C_target) {
  q <- mean(C_ref < tau_C)
  as.numeric(stats::quantile(C_target, q))
}
