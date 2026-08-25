# LDscnR (development version)

## Breaking changes

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

## Bug fixes

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

## New features

* `ld_scan()` and its four stages -- `ld_null_from_p()`, `ld_gate()`,
  `ld_region_scan()`, `ld_region_c2()` -- take p-values from the observed data
  and from user-supplied permuted datasets, making the pipeline independent of
  the association engine and of how the null is constructed.

* `null_fdr()`, `calibrate_tauc()` and `calibrate_lmin()` are documented as
  diagnostics rather than as the way to choose an operating point.

* `ld_cscore()` now credits the C-score to Fang et al. (2021),
  <doi:10.1093/molbev/msab144>, which introduced it.
