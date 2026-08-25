#' Build a null bundle from your own observed and permuted p-values
#'
#' The engine-agnostic entry point to the LD-scan pipeline. Where
#' [structured_null()] generates its own surrogates by running a fast EMMAX scan
#' on MVN-resampled phenotypes, this constructor takes p-values you have already
#' computed -- from EMMAX, LFMM, BayPass, a GLM, anything -- for the observed data
#' and for `B` permuted (or otherwise surrogate) datasets, and reduces them to the
#' same `ld_null` bundle the rest of the pipeline consumes.
#'
#' This is what makes the null *yours*: the pipeline never sees your phenotype,
#' your covariates or your association model, so the surrogate construction is
#' entirely under your control. Permute in whatever way makes your null hypothesis
#' the one you mean -- within-locality label shuffling, a spatial kernel draw, a
#' case-control swap -- run your scan on each, and hand the p-values here. The
#' only requirement is that every surrogate is scanned with the *same* engine and
#' settings as the observed data, or the comparison is not like-for-like.
#'
#' Surrogate C-scores are stored *sparsely* (`C > 0` only), so `B = 200`
#' whole-genome surrogates cost a few MB, not gigabytes. `p_perm` is consumed one
#' surrogate at a time and reduced immediately, so a `function(b)` accessor lets
#' you stream surrogates off disk without ever holding the full
#' markers-by-`B` p-value matrix in memory.
#'
#' @param p_obs Numeric vector of observed per-SNP p-values, aligned to the rows
#'   of `ld_ws` (i.e. `length(p_obs) == nrow(ld_ws)`).
#' @param p_perm The surrogate p-values, in any of three forms:
#'   * a numeric **matrix** with SNPs in rows and surrogates in columns;
#'   * a **list** of numeric vectors, one per surrogate;
#'   * a **function** `f(b)` returning the p-value vector of surrogate `b`, in
#'     which case `B` must be given. Use this to stream from disk.
#' @param ld_ws Local-LD support matrix (SNPs x `ld_w` windows; see [ld_cscore()]).
#' @param B Number of surrogates. Required only when `p_perm` is a function;
#'   otherwise inferred, and checked if supplied.
#' @param alpha,rho,qstar C-score parameters, passed to [ld_cscore()]. These must
#'   match between observed and surrogates -- they are applied identically here.
#' @param basis A label for how the surrogates were constructed (e.g.
#'   `"within-locality permutation"`). Carried through to the printed output and
#'   the gate table; it is documentation, not behaviour.
#' @param engine A label for the association engine that produced the p-values
#'   (e.g. `"EMMAX"`, `"LFMM"`). Also documentation.
#' @param cores Number of cores for reducing surrogates to C-scores
#'   ([parallel::mclapply()], so >1 is Unix-only). The reduction is the expensive
#'   step; the scan itself is cheap.
#' @param verbose Print progress.
#'
#' @return An `ld_null` object -- the same class [structured_null()] returns, so
#'   everything downstream applies: `C_obs` (full observed C-vector), `C_surr`
#'   (list of `B` sparse surrogate C-vectors), `universe` (every marker lit up by
#'   the observed or any surrogate), and the settings.
#' @seealso [ld_gate()], [ld_region_scan()], [ld_scan()], [structured_null()]
#' @export
ld_null_from_p <- function(p_obs, p_perm, ld_ws, B = NULL, alpha = 0.05,
                           rho = colnames(ld_ws), qstar = seq(0, 0.95, by = 0.05),
                           basis = "user-supplied", engine = "user-supplied",
                           cores = 1L, verbose = TRUE) {
  if (length(p_obs) != nrow(ld_ws))
    stop(sprintf("`p_obs` has %d values but `ld_ws` has %d rows -- they must be aligned.",
                 length(p_obs), nrow(ld_ws)))
  if (is.null(rownames(ld_ws))) stop("`ld_ws` must have rownames (the marker IDs).")
  .check_p_names(p_obs, ld_ws, "`p_obs`")

  ## ---- normalise the three accepted forms of `p_perm` to a get(b) accessor ----
  if (is.function(p_perm)) {
    if (is.null(B)) stop("`B` is required when `p_perm` is a function.")
    get_b <- p_perm
  } else if (is.matrix(p_perm)) {
    if (nrow(p_perm) != nrow(ld_ws))
      stop(sprintf("`p_perm` has %d rows but `ld_ws` has %d -- surrogates must be SNPs-by-B.",
                   nrow(p_perm), nrow(ld_ws)))
    if (!is.null(B) && B != ncol(p_perm))
      stop(sprintf("`B` = %d disagrees with ncol(p_perm) = %d.", B, ncol(p_perm)))
    B <- ncol(p_perm); get_b <- function(b) p_perm[, b]
  } else if (is.list(p_perm)) {
    if (!is.null(B) && B != length(p_perm))
      stop(sprintf("`B` = %d disagrees with length(p_perm) = %d.", B, length(p_perm)))
    B <- length(p_perm); get_b <- function(b) p_perm[[b]]
  } else {
    stop("`p_perm` must be a matrix (SNPs x B), a list of vectors, or a function(b).")
  }
  if (B < 1L) stop("Need at least one surrogate.")
  if (B < 20L)
    warning(sprintf(paste("Only %d surrogate(s): the smallest attainable region p-value is",
                          "1/(1+B) = %.3f, and the gate's median is read off %d value(s).",
                          "B >= 100 is recommended."), B, 1 / (1 + B), B), call. = FALSE)

  if (verbose) { cat(sprintf("[ld_null_from_p] reducing %d surrogate(s) to sparse C-scores on %d core(s)...\n",
                             B, cores)); utils::flush.console() }
  C_obs <- ld_cscore(p_obs, ld_ws, alpha, rho, qstar)
  one <- function(b) {
    pv <- get_b(b)
    if (length(pv) != nrow(ld_ws))
      stop(sprintf("Surrogate %d has %d p-values but `ld_ws` has %d rows.", b, length(pv), nrow(ld_ws)))
    .check_p_names(pv, ld_ws, sprintf("Surrogate %d", b))
    C <- ld_cscore(pv, ld_ws, alpha, rho, qstar)
    C[C > 0]                                   # sparse: only the markers that light up
  }
  C_surr <- if (cores > 1L) parallel::mclapply(seq_len(B), one, mc.cores = cores)
            else lapply(seq_len(B), one)
  bad <- which(!vapply(C_surr, is.numeric, logical(1)))
  if (length(bad))                             # mclapply returns try-errors rather than stopping
    stop(sprintf("Surrogate(s) %s failed to reduce: %s",
                 paste(bad, collapse = ", "), as.character(C_surr[[bad[1]]])))

  universe <- unique(c(names(C_obs)[C_obs > 0], unlist(lapply(C_surr, names), use.names = FALSE)))
  if (verbose) cat(sprintf("[ld_null_from_p] universe = %d marker(s) with C > 0 in the observed or any surrogate\n",
                           length(universe)))
  structure(list(C_obs = C_obs, C_surr = C_surr, universe = universe,
                 basis = basis, engine = engine, B = B,
                 params = list(alpha = alpha, rho = rho, qstar = qstar)),
            class = "ld_null")
}

#' @export
print.ld_null <- function(x, ...) {
  cat(sprintf("<ld_null> B = %d surrogate(s) | basis: %s | engine: %s\n",
              x$B, x$basis %||% "?", x$engine %||% "?"))
  cat(sprintf("  observed: %d marker(s) with C > 0 of %d | universe: %d | p-floor 1/(1+B) = %.4f\n",
              sum(x$C_obs > 0), length(x$C_obs), length(x$universe), 1 / (1 + x$B)))
  invisible(x)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
