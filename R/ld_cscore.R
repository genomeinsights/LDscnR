#' Per-SNP consistency C-score
#'
#' The LD-scan C-score: for each SNP, the fraction of `(rho, q*)` cells in which
#' it is BOTH a local-LD candidate (its `ld_w` at window `rho` exceeds the `q*`
#' quantile of that window) AND significant (Benjamini-Hochberg FDR `< alpha`
#' among the candidates of that cell). It integrates over the two nuisance
#' parameters that have a natural `[0, 1]` domain -- the `ld_w` window `rho` and
#' the local-LD quantile `q*` -- at a FIXED significance level. A high C means a
#' SNP is consistently prioritised-and-significant however the LD window and
#' stringency are set, concentrating on true, LD-supported signal.
#'
#' `alpha` is fixed (default 0.05) rather than swept: `rho` and `q*` span `[0, 1]`
#' so integrating over them is assumption-free, but `alpha` has no natural range,
#' and the structure-aware null ([structured_null()]) prices any
#' anticonservativeness into the `tau_C` threshold -- so sweeping `alpha` is
#' unnecessary (and empirically neutral). Passing a *vector* `alpha` folds it into
#' the integration too (an extra sweep axis), which is useful only as a
#' robustness check that the fixed choice does not matter; keep the default for
#' analysis.
#'
#' The C-score is a per-SNP statistic: clustering into outlier regions
#' ([ld_regions()]) and the region-size filter `l_min` are applied downstream,
#' not folded in.
#'
#' @section Attribution:
#' **The C-score is not original to this package.** It was introduced by
#' Fang, Kemppainen, Momigliano & Merila (2021), *Population structure limits
#' parallel evolution in sticklebacks*, Molecular Biology and Evolution
#' 38(10):4205-4221, \doi{10.1093/molbev/msab144}, which defines the C-score by
#' name together with `tau_C = 0.05`, `l_min = 10`, a 500 kb single-linkage
#' distance cap, and the definition of a region as a cluster of significant
#' loci. There the score integrates over a **menu of 144 tests** -- three
#' clustering parameterisations crossed with four multiple-testing corrections,
#' one of them the LD-cluster permutation of Li, Kemppainen, Rastas & Merila
#' (2018), \doi{10.1111/1755-0998.12893}, made affordable by complexity
#' reduction. Cite Fang et al. (2021) when you use the C-score.
#'
#' What this implementation changes is the shape of the integration, not the
#' idea. Integration here is **continuous** over `rho` and `q*` -- two nuisance
#' parameters with a natural `[0, 1]` domain -- rather than over a discrete menu
#' of analyses; and the score stays a **per-SNP** statistic, whereas Fang et al.
#' test cluster-level units (SMLAs) directly. Keeping it per-SNP is what allows
#' `tau_C` and `l_min` to be swept as a second, nested scale downstream
#' ([ld_region_c2()]).
#'
#'
#' @param p Numeric vector of per-SNP association p-values (from any method),
#'   aligned to the rows of `ld_ws`.
#' @param ld_ws Numeric matrix of local-LD support: SNPs in rows, `ld_w` windows
#'   (`rho`) in columns, `rownames` = markers. See [compute_ld_w()].
#' @param alpha Within-candidate BH-FDR level (fixed scalar; default 0.05). A
#'   vector adds an extra integration axis (robustness check only).
#' @param rho `ld_w` window columns of `ld_ws` to integrate over (default: all).
#' @param qstar Local-LD quantile thresholds to integrate over
#'   (default `seq(0, 0.95, by = 0.05)`).
#'
#' @return Named numeric vector of per-SNP C-scores in `[0, 1]`
#'   (`names` = markers). Gate with a threshold `tau_C`
#'   (`names(C)[C >= tau_C]`) calibrated by [calibrate_tauc()].
#' @references
#' Fang B, Kemppainen P, Momigliano P, Merila J (2021). Population structure
#' limits parallel evolution in sticklebacks. *Molecular Biology and Evolution*
#' 38(10):4205-4221. \doi{10.1093/molbev/msab144}
#'
#' Li Z, Kemppainen P, Rastas P, Merila J (2018). Linkage disequilibrium
#' clustering-based approach for association mapping with tightly linked
#' genomewide data. *Molecular Ecology Resources*. \doi{10.1111/1755-0998.12893}
#' @seealso [ld_regions()], [ld_region_c2()], [structured_null()]
#' @export
ld_cscore <- function(p, ld_ws, alpha = 0.05,
                      rho = colnames(ld_ws), qstar = seq(0, 0.95, by = 0.05)) {
  stopifnot(length(p) == nrow(ld_ws))
  ncell <- length(rho) * length(qstar) * length(alpha)
  cnt <- integer(nrow(ld_ws))
  for (rc in rho) {
    lw <- ld_ws[, rc]
    for (q in qstar) {
      thr <- stats::quantile(lw, q, na.rm = TRUE)
      cand <- which(lw >= thr)
      if (!length(cand)) next
      qv <- stats::p.adjust(p[cand], "BH")
      for (al in alpha) {
        hit <- cand[which(qv < al)]          # which() drops NA q-values
        if (length(hit)) cnt[hit] <- cnt[hit] + 1L
      }
    }
  }
  stats::setNames(cnt / ncell, rownames(ld_ws))
}
