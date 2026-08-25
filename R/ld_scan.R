#' LD-aware outlier regions from your own scan, end to end
#'
#' The one-call entry point. Give it per-SNP p-values for your observed data and
#' for `B` permuted datasets, and it runs the whole pipeline: reduce both to
#' C-scores ([ld_null_from_p()]), build the LD-edge graph ([ld_edges()]), check
#' the null against the observed background ([ld_gate()]), call regions and test
#' each at its own locus ([ld_region_scan()]), and -- unless you turn it off --
#' rank them by the second-tier C-score over the `(tau, l_min)` grid
#' ([ld_region_c2()]).
#'
#' It is engine-agnostic by construction. The pipeline never sees your phenotype
#' or your association model, only p-values, so EMMAX, LFMM, BayPass, a plain GLM
#' or anything else all work, and two engines can be compared on identical
#' footing. The one requirement is that **every surrogate is scanned with the same
#' engine and settings as the observed data** -- otherwise the comparison is not
#' like-for-like and nothing downstream can detect it.
#'
#' How you permute is the scientific choice, and it is yours. The surrogates
#' define what "no signal" means: shuffling labels within localities holds
#' between-locality composition fixed and therefore treats a clinal signal as
#' null, whereas shuffling globally breaks the cline and treats it as real. Decide
#' which hypothesis you are testing, permute accordingly, and state it -- the
#' `basis` argument is carried through every printed output for exactly that
#' reason.
#'
#' The result prints the **gate before the regions**, deliberately. A p-value from
#' a null whose surrogates reproduce the observed background is not interpretable,
#' and no amount of downstream care recovers it.
#'
#' @param p_obs Observed per-SNP p-values, aligned to the rows of `ld_ws`.
#' @param p_perm Surrogate p-values: a SNPs-by-`B` matrix, a list of vectors, or a
#'   `function(b)` streaming them from disk (see [ld_null_from_p()]).
#' @param ld_ws Local-LD support matrix (SNPs x `ld_w` windows; see [compute_ld_w()]).
#' @param map data.frame/data.table with `marker`, `Chr`, `Pos`.
#' @param GTs Genotype dosage matrix (individuals x SNPs, column names = markers),
#'   used only to compute the r^2 between markers for the edge graph.
#' @param decay_sum Per-chromosome LD-decay summary (see [compute_LD_decay()]).
#' @param tau C-score threshold for forming regions (default 0.05).
#' @param l_min Minimum region size in markers (default 3).
#' @param fdr BH level for `q_R` (default 0.05).
#' @param alpha,rho,qstar C-score parameters (see [ld_cscore()]).
#' @param rho_ld,dcap,rho_d Clustering parameters (see [ld_edges()]).
#' @param B Number of surrogates; required only if `p_perm` is a function.
#' @param c2 Compute the second-tier C-score ranking (default `TRUE`). This
#'   re-runs the region test over the whole grid and is by far the most expensive
#'   step; set `FALSE` for a quick look.
#' @param tau_grid,lmin_grid The `C^(2)` grid (see [ld_region_c2()]).
#' @param basis,engine Labels describing how surrogates were built and which
#'   engine produced the p-values. Documentation, carried through the output.
#' @param cores Cores for reducing surrogates to C-scores (Unix only).
#' @param verbose Print progress.
#'
#' @return An `ld_scan` object with `gate`, `regions` (the `ld_region_scan`),
#'   `c2` (or `NULL`), `null`, `edges` and `params`.
#'
#' @examples
#' \dontrun{
#' ## p_obs: your observed scan. p_list: one p-value vector per permuted dataset,
#' ## each produced by re-running THAT SAME scan on permuted labels.
#' fit <- ld_scan(p_obs, p_list, ld_ws, map, GTs, decay_sum,
#'                tau = 0.05, l_min = 3, rho_ld = 0.6,
#'                basis = "within-locality permutation", engine = "EMMAX")
#' fit                          # gate first, then the significant regions
#' fit$gate                     # is this null usable at all?
#' fit$regions$regions          # per-region s_R, p_R, q_R
#' fit$c2$regions               # the same regions ranked by C^(2)
#' }
#' @seealso [ld_null_from_p()], [ld_gate()], [ld_region_scan()], [ld_region_c2()]
#' @export
ld_scan <- function(p_obs, p_perm, ld_ws, map, GTs, decay_sum,
                    tau = 0.05, l_min = 3L, fdr = 0.05,
                    alpha = 0.05, rho = colnames(ld_ws), qstar = seq(0, 0.95, by = 0.05),
                    rho_ld = 0.9, dcap = 5e5, rho_d = NULL, B = NULL,
                    c2 = TRUE, tau_grid = seq(0.02, 0.5, by = 0.02),
                    lmin_grid = c(1L, 2L, 3L, 5L, 10L, 15L, 20L),
                    basis = "user-supplied", engine = "user-supplied",
                    cores = 1L, verbose = TRUE) {
  mp <- data.table::as.data.table(map)
  if (!all(c("marker", "Chr", "Pos") %in% names(mp)))
    stop("`map` needs columns `marker`, `Chr` and `Pos`.")

  null <- ld_null_from_p(p_obs, p_perm, ld_ws, B = B, alpha = alpha, rho = rho,
                         qstar = qstar, basis = basis, engine = engine,
                         cores = cores, verbose = verbose)
  if (!length(null$universe))
    stop("No marker reaches C > 0 in the observed data or any surrogate -- nothing to cluster.")

  if (verbose) { cat(sprintf("[ld_scan] building the LD-edge graph over %d marker(s)...\n",
                             length(null$universe))); utils::flush.console() }
  edges <- ld_edges(null$universe, GTs, mp[, c("marker", "Chr", "Pos"), with = FALSE],
                    decay_sum, rho_ld = rho_ld, dcap = dcap, rho_d = rho_d)

  gate <- ld_gate(null, edges, tau = tau, l_min = l_min)   # warns if the null is mis-specified
  scan <- ld_region_scan(null, edges, tau = tau, l_min = l_min, fdr = fdr)
  c2r <- if (isTRUE(c2)) ld_region_c2(null, edges, tau_grid = tau_grid, lmin_grid = lmin_grid,
                                      fdr = fdr, anchor_tau = tau, anchor_lmin = l_min,
                                      verbose = verbose) else NULL

  structure(list(gate = gate, regions = scan, c2 = c2r, null = null, edges = edges,
                 params = list(tau = tau, l_min = l_min, fdr = fdr, rho_ld = rho_ld,
                               dcap = dcap, alpha = alpha, B = null$B,
                               basis = basis, engine = engine)),
            class = "ld_scan")
}

#' @export
print.ld_scan <- function(x, ...) {
  p <- x$params
  cat(sprintf("<ld_scan> %s / %s | tau = %.2f, l_min = %d, rho_ld = %.2f, B = %d\n",
              p$engine, p$basis, p$tau, p$l_min, p$rho_ld, p$B))
  cat("\n-- gate (read this first) --\n"); print(x$gate)
  if (!all(x$gate$pass))
    cat("\n!! This null failed the gate. The regions below are NOT interpretable.\n")
  cat("\n-- regions --\n"); print(x$regions)
  if (!is.null(x$c2)) { cat("\n-- second-tier ranking --\n"); print(x$c2) }
  invisible(x)
}
