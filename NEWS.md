# LDscnR 0.9.3

## Exact fixed-effect covariates in EMMAX

* `emmax_setup()` gains `covar`: a vector, factor, matrix or data.frame with
  one row per individual. It is turned into a design matrix once
  (`model.matrix()`, so factors become indicator columns; zero-variance
  columns are dropped), and the fixed effects become `X = [1, covar]`. REML
  is fitted with that `X`, and `emmax_fast()` tests each SNP by generalised
  least squares conditional on it (Frisch-Waugh-Lovell in the whitened
  space), on `n - ncol(X) - 1` residual degrees of freedom. The prep keeps
  its covariates for every phenotype, which is what permutation loops need.
  Covariates that are rank-deficient together with the intercept are
  rejected. Without `covar`, p-values are bit-for-bit those of 0.9.2.
* Regressing the phenotype on the covariates first, the previous approach,
  leaves the genotypes unadjusted. That shrinks the test statistic of a SNP
  by about `sqrt(1 - R^2(SNP ~ covar))` and loses power exactly for SNPs
  correlated with a covariate, and with `n - 2` degrees of freedom it is
  slightly anti-conservative.
* `emmax()` gains `covar` and `covar_method = c("gls", "residualise")`.
  `"gls"` (default) is the exact model above; `"residualise"` is the old
  method. `Covar` is deprecated: it still works, selects `"residualise"` and
  gives the same numbers as before, with a once-per-session warning.
* `emmax_fast()` accepts a matrix of phenotypes (one per column, e.g.
  permutations) and returns a SNPs x phenotypes matrix. Its result carries
  the attributes `df2` and `covar_names`. A SNP with nothing left to test
  (monomorphic, or explained completely by the covariates) gets `NA`, not
  `NaN`. The prep adds `X`, `Xt_X`, `covar_names` and `df2`; `Xo` and `xot`
  are kept, and preps saved by earlier versions still work.
* `ld_outlier_regions()`, `structured_null()` and `ld_cscore_scan()` pass a
  new `covar` argument to `emmax_setup()`. `structured_null()` makes its
  surrogates orthogonal to `[1, covar, y]` instead of `y` only when the model
  has covariates. `ld_outlier_clusters()` gains `covar` for
  `null = "permutation"`, which permutes `Y` with the covariates fixed; its
  permuted scans now call `emmax_fast()` instead of an inline copy of it
  (results without `covar` are unchanged).

## Fixes

* `ld_unit_matrix()`, `ld_prune_and_eMLG()` and `make_eMLGs()` look up each
  cluster's genotype columns by position, matched once for all clusters,
  instead of subsetting `GTs` by marker name once per cluster. Subsetting by
  name re-hashes every column name on each call (about 40 ms at 1M markers),
  so building tens of thousands of units took about half an hour; it now
  takes seconds. Results are identical.
* `ld_complexity_reduction()` again stores `params$rho` (equal to `r2_rho`),
  so callers written for 0.9.1 that pass `min_r2_rho = stage1$params$rho` to
  `ld_prune_and_eMLG()` keep working. In 0.9.2 they passed `NULL` and failed
  on the first chromosome with flagged clusters with `'names' attribute [n]
  must be the same length as the vector [0]`. `ld_prune_and_eMLG()` now
  rejects a `NULL` or out-of-range `r2_rho` or `dist_rho` with a clear
  message.

# LDscnR 0.9.2

* `ld_complexity_reduction()` gains an optional `r2_th` argument, a fixed r2
  clustering threshold applied to every chromosome instead of the
  decay-relative threshold derived from `r2_rho`. The default (`NULL`) keeps the
  decay-relative threshold, so existing results are unchanged.

* `ld_outlier_test()` gains the same optional `r2_th` argument, a fixed r2
  merge threshold for Stage-2 assembly (`assembly = "stage2_discovered"`).
  If it is `NULL` and `stage1` was built with a fixed `r2_th`, Stage-2
  assembly uses that same value. Otherwise it keeps the decay-relative
  threshold. `ld_outlier_perm()` passes the threshold on to every surrogate.

* `ld_outlier_scan()` passes `r2_th` on to `ld_outlier_test()`, and
  `ld_prune_and_eMLG()` uses a fixed `r2_th` stored on `stage1` when none is
  given.

## Consistent argument names

One parameter now has one name across the package. The old names still work
but give a deprecation warning.

* Fixed r2 threshold `r2_th`: was `min_r2` in `ld_prune_and_eMLG()` and
  `dynamic_cut_eMLG()`, and `r2_link` in `ld_outlier_clusters()`.
* Decay-relative level for the r2 threshold `r2_rho`: was `rho` in
  `ld_complexity_reduction()`, `min_r2_rho` in `ld_prune_and_eMLG()`,
  `rho_ld` in `ld_outlier_clusters()` and `rho_r2` in `score_thresholds()`.
* Decay-relative level for a distance `dist_rho`: was `rho` in
  `ld_prune_and_eMLG()` and `compute_ld_w()`, and `rho_d` in
  `ld_outlier_clusters()` and `score_thresholds()`.
* `alpha` and `cluster_alpha`: were `fdr` and `cluster_fdr` in
  `ld_outlier_clusters()`, and `fdr` in `rmsc_threshold()`.
* `size_floor` and `cores`: were `l_min` and `ncores` in `make_eMLGs()`.
* `LD_decay`: was `ld_decay` in `compute_ld_w()`.
* `stage1`: was `pruned_stage1` in `plot_pruning_comparison()`.
* `GTs`: was `geno` in `create_gds_from_geno()` and `run_ohtadstats()`.
* `Y` (phenotype): was `y` in `emmax_fast()`.
* `slide`: was `slide_win_ld` in `compute_ld_w()` and `get_el()`.
* `ld_method`: was `method` in `get_el()`.
* `maf_min_decay`: was `min_maf_decay` in `compute_LD_decay()`.
* `params` of `ld_complexity_reduction()`, `ld_prune_and_eMLG()` and
  `ld_outlier_test()` use the new names (`r2_rho`, `r2_th`, `dist_rho`,
  `r2_th_resolved`). Stage-1 objects saved by earlier versions, which store
  `params$rho`, are still read correctly.

## Changed default

* `get_el()` now defaults to `ld_method = "corr"`, the same as
  `compute_LD_decay()` and `compute_ld_w()`. It previously defaulted to
  `"r"`. Every internal caller already passed the method explicitly, so
  results from the main workflow do not change.

# LDscnR 0.9.1

* Truth-scoring helpers for simulated data are now exported: `flag_true_qtns()`
  (detectable QTN by MAF and share of additive variance), `qtn_ld_table()`
  (r^2 and distance between markers and QTN) and `score_thresholds()`
  (decay-relative match limits). They are method-agnostic and can score any
  scan of simulated data, including marker-wise and window-based scans.

* `emmax_setup()` / `emmax_fast()` (fast EMMAX reusing one relationship-matrix
  transformation across phenotypes) are documented as part of the public API.

* `coef_ld_dec()`: removed a stray no-op line. The `sim_ex` documentation now
  describes its `map` component.

# LDscnR 0.9.0

## Current outlier workflow

* `ld_unit_matrix()` builds one testable variable per phenotype-blind Stage-1
  unit. `ld_outlier_test()` applies multiple-testing correction to unit-level
  p-values and assembles significant units into reported regions.
  `ld_outlier_perm()` repeats the same analysis under user-supplied surrogate
  p-values to compare discovery counts with a study-appropriate null.

* `ld_region_rotation()` separately tests overlap with external annotations
  by independently relocating regions while preserving their spans. This is
  not a phenotype null or a test of individual-region significance.

* The C-score tutorials have been moved outside the installed vignette set.
  The earlier implementation remains available for reproducibility, but the
  Stage-1 workflow is the method described in the current manuscript.

* Documentation now distinguishes LD reduction from independence, clarifies
  the dataset-specific size floor and the LFMM/Simes route, and avoids
  treating independent relocation draws as independent discoveries.

## Historical changes to retained C-score code

The items below describe the earlier implementation, not the recommended
Stage-1 outlier workflow.

### Breaking changes

* The default `rho_ld` is now **0.75** everywhere (`ld_edges()`, `ld_scan()`,
  `ld_outlier_regions()`, `ld_region_stability()`), replacing 0.9. The package
  previously disagreed with itself: four formal defaults at 0.9, `ld_scan()`'s
  own example at 0.6, and every analysis built on the package at 0.6 or 0.75.
  0.75 is the value used by both the stickleback and the simulation analyses.

  On the stickleback panel the choice turns out not to matter across this range:
  0.60 and 0.75 give the same 17 regions, the same `q_R` values and the same
  null-gate verdicts, differing only by a single marker joining a single region,
  even though the r^2 link moves from 0.405 to 0.273. Insensitivity is bounded
  though -- by 0.9 the link is loose enough to overmerge distinct peaks -- so the
  documentation for `ld_edges()` now states the direction (higher `rho_ld` means
  a lower r^2 link and more merging) and the range over which it was checked,
  in place of the previous unqualified "shown insensitive" claim.

### Bug fixes

* `ld_cscore()` no longer errors on NA p-values. `p.adjust()` propagates NA, so
  the hit index contained NA and assignment failed with "NAs are not allowed in
  subscripted assignments". NA p-values are routine (monomorphic markers, EMMAX
  failures, LFMM missingness) and the failure surfaced mid-way through surrogate
  loops. C-scores for all non-NA markers are unchanged.

* `ld_cscore()` and `ld_null_from_p()` now reject named p-value vectors whose
  names disagree with `rownames(ld_ws)`, instead of silently scoring markers
  against other markers' p-values. Unnamed vectors remain legal.

* `ld_gate()` no longer reports a basis as passing when the observed data
  produced no regions. "Observed found nothing, surrogates found plenty" is a
  gate failure, and was previously reported as a pass because the ratio was `NA`.

### New features

* `ld_scan()` and its four stages -- `ld_null_from_p()`, `ld_gate()`,
  `ld_region_scan()`, `ld_region_c2()` -- take p-values from the observed data
  and from user-supplied permuted datasets, making the pipeline independent of
  the association engine and of how the null is constructed.

* `null_fdr()`, `calibrate_tauc()` and `calibrate_lmin()` are documented as
  diagnostics rather than as the way to choose an operating point.

* `ld_cscore()` now credits the C-score to Fang et al. (2021),
  <doi:10.1093/molbev/msab144>, which introduced it.
