# Archived vignette source

`LDscnR_outlier_regions_from_pvalues.Rmd` documents `ld_scan()` and the
consistency C-score family it wraps -- an earlier, genuinely different
approach to the outlier-region problem than the package's current primary
method, and not the method behind any current manuscript result.

It is kept here, outside `vignettes/`, so its source remains readable and
part of the repository's history without being built or installed as one
of the package's tutorial vignettes (this directory is excluded from the
built package via `.Rbuildignore`). Its content was reviewed and corrected
for known factual issues before archiving (an inaccurate claim that
BayPass supplies p-values directly, a leftover "prefer `ld_scan()`"
recommendation that contradicted its own now-outdated status, and an
undercaveated presentation of `q_R`) rather than left as-is, since even an
archived document should not state things that are simply wrong.

For the current, primary outlier-analysis workflow, see
`vignette("LDscnR_stage1_outlier_regions")` (in `vignettes/`) or
`vignette("LDscnR_outlier_analysis")` (also in `vignettes/`, itself
labelled as the older, non-primary C-score method but still installed,
since it does not share this file's issues).
