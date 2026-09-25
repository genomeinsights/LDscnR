#' Permutation null for `ld_outlier_test()`: how many discoveries under no signal
#'
#' Not separate machinery. This runs the SAME [ld_outlier_test()] pipeline once
#' per column of `p_perm` and reports the surrogate distribution of discovery
#' counts against the observed one -- the null falls out of code that already
#' exists rather than being reimplemented for validation.
#'
#' What "no signal" means is entirely determined by how `p_perm`'s surrogates
#' were generated, and that is the caller's design, not this function's: which
#' unit is permuted, what structure is held fixed. A scheme that breaks the
#' structure your association model corrects for produces an anticonservative
#' null that no amount of downstream care repairs -- this function cannot detect
#' that and does not try to.
#'
#' @param obs An `ld_outlier_test` result -- the observed call, whose
#'   `$params` are reused for every surrogate so the null is calibrated under
#'   exactly the same `statistic`/`size_floor`/`alpha`/`assembly` as the
#'   observed count.
#' @param stage1,map,GTs,LD_decay Same as [ld_outlier_test()]; only `GTs`/
#'   `LD_decay` are actually needed here for `assembly = "stage2_discovered"`,
#'   passed through unchanged for every surrogate.
#' @param p_perm Surrogate p-values, same shape `p_obs` had in `obs`: a
#'   markers/units-by-`B` matrix, a list of vectors, or a `function(b)`.
#' @param B Number of surrogates; required only if `p_perm` is a function.
#' @param level `"units"` (default; count significant units per surrogate) or
#'   `"regions"` (count assembled regions per surrogate) -- the two need not
#'   agree, since assembly can merge several significant units into one region.
#' @param cores Cores for the surrogate loop (Unix only).
#' @param verbose Print progress (default `TRUE`).
#'
#' @return An `ld_outlier_perm` object: `observed` (the count from `obs`,
#'   at `level`), `surrogates` (integer vector, length `B`), `p` (one-sided,
#'   `(1 + #surrogates >= observed) / (B + 1)`, always well-defined, including
#'   when `observed = 0`), `null_obs_ratio`
#'   (`mean(surrogates) / observed` -- a null-to-observed DISCOVERY-COUNT
#'   ratio, not an observed false-discovery rate: it says how many
#'   discoveries a structure-preserving null generates relative to how many
#'   were actually found, nothing about which of the observed discoveries,
#'   if any, are real. A value near or above 1 warns that the observed count
#'   is unremarkable next to this null; a value well below 1 does not by
#'   itself establish that individual reported regions are genuine, or that
#'   the nominal BH false-discovery rate is calibrated -- both `p` and
#'   `null_obs_ratio` describe the AGGREGATE discovery count, at `level`,
#'   under whatever structure `p_perm`'s surrogates were built to preserve;
#'   neither is a per-region statement. `NA` when `observed = 0`: the ratio
#'   is undefined there, not zero -- there is no observed count for the
#'   null to be a fraction of, regardless of the surrogate mean. `p`
#'   remains valid and interpretable in that case), `realised_fdr`
#'   (deprecated alias for `null_obs_ratio`, identical value including the
#'   `NA` case, kept for scripts written against the old name -- new code
#'   should read `null_obs_ratio`), `params`.
#'
#' @seealso [ld_outlier_test()], [ld_region_rotation()]
#' @export
ld_outlier_perm <- function(obs, stage1, map, p_perm,
                            GTs = NULL, LD_decay = NULL,
                            B = NULL, level = c("units", "regions"),
                            cores = 1L, verbose = TRUE) {
  level <- match.arg(level)
  if (!inherits(obs, "ld_outlier_test"))
    stop("`obs` must be the result of ld_outlier_test().")
  p <- obs$params

  ## normalise p_perm to a get(b) accessor, same three forms ld_null_from_p() accepts
  if (is.function(p_perm)) {
    if (is.null(B)) stop("`B` is required when `p_perm` is a function.")
    get_b <- p_perm
  } else if (is.matrix(p_perm)) {
    B <- ncol(p_perm)
    get_b <- function(b) p_perm[, b]
  } else if (is.list(p_perm)) {
    B <- length(p_perm)
    get_b <- function(b) p_perm[[b]]
  } else stop("`p_perm` must be a matrix, a list of vectors, or a function(b).")

  ## THE FAST PATH: when level = "units", region assembly is never touched at all --
  ## not run and thrown away, not run with a shortcut, simply never called. Assembly
  ## cannot change which units are significant (it runs strictly after BH), so for a
  ## unit-level null it is pure waste, and for assembly = "stage2_discovered" it was the
  ## DOMINANT cost per surrogate (a full ld_prune_and_eMLG() call, eMLG computation
  ## included). Only level = "regions" needs the real assembly, once per surrogate.
  ##
  ## `units_base` is built ONCE, outside the surrogate loop, and reused for every draw --
  ## it depends only on (stage1, map, size_floor), never on the surrogate phenotype, so
  ## rebuilding it per surrogate was pure waste even after .ld_outlier_units() itself was
  ## fixed from 27.7s to 0.5s (a separate, more fundamental fix: a named-vector lookup by
  ## marker name over ~790k names is a linear scan per lookup in base R, not a hash).
  units_base <- .ld_outlier_units(stage1, map, p$size_floor)
  one <- function(b) {
    if (level == "units") {
      u <- .ld_outlier_tested_units(stage1, map, get_b(b), p$statistic, p$size_floor, p$alpha,
                                    units = units_base)
      return(sum(u$significant))
    }
    r <- ld_outlier_test(stage1, map, get_b(b), statistic = p$statistic,
                         size_floor = p$size_floor, alpha = p$alpha,
                         assembly = p$assembly, GTs = GTs, LD_decay = LD_decay,
                         score_threshold = p$score_threshold,
                         distance_threshold = p$distance_threshold, gap = p$gap)
    nrow(r$regions)
  }
  if (verbose) { cat(sprintf("[ld_outlier_perm] %d surrogates, level = %s\n", B, level))
    utils::flush.console() }
  surrogates <- if (cores > 1L) unlist(parallel::mclapply(seq_len(B), one, mc.cores = cores))
                else vapply(seq_len(B), one, 0L)

  observed <- if (level == "units") sum(obs$units$significant) else nrow(obs$regions)
  ## Undefined, not zero, when nothing was observed: 0/0 is not "the null
  ## found nothing relative to the observed count", it is "there was no
  ## observed count to compare the null to" -- a different, degenerate
  ## case. The count-tail p-value stays well-defined either way (surrogates
  ## >= 0 is true whenever observed = 0, correctly giving p close to 1) and
  ## is not touched by this.
  null_obs_ratio <- if (observed > 0) mean(surrogates) / observed else NA_real_
  structure(list(
    observed = observed, surrogates = surrogates,
    p = (1 + sum(surrogates >= observed)) / (B + 1),
    null_obs_ratio = null_obs_ratio,
    realised_fdr = null_obs_ratio,  ## deprecated alias -- see @return
    params = list(level = level, B = B, test_params = p)
  ), class = "ld_outlier_perm")
}

#' @export
print.ld_outlier_perm <- function(x, ...) {
  ratio_str <- if (is.na(x$null_obs_ratio)) "NA (0 observed)" else sprintf("%.1f%%", 100 * x$null_obs_ratio)
  cat(sprintf("<ld_outlier_perm> observed %d | surrogate mean %.2f | p = %.4f | null/obs ratio %s\n",
              x$observed, mean(x$surrogates), x$p, ratio_str))
  invisible(x)
}
