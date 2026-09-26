# LDscnR

**Chromosome-wise LD-decay estimation, LD-based marker pruning, eMLG generation, and LD-aware outlier-region detection**

`LDscnR` estimates linkage disequilibrium (LD) decay and groups correlated markers without using the phenotype. These Stage-1 groups can supply a less redundant marker set for a relationship matrix, a consensus genotype for each group, or the units for association testing. Significant Stage-1 units can then be assembled into outlier regions for reporting. Grouping reduces redundancy; it does not prove that different groups are statistically independent. See the vignettes for worked examples.

------------------------------------------------------------------------

## Key features

- **Chromosome-specific LD-decay estimation**

- **Background LD estimation** from inter-chromosomal SNP pairs

- **Decay-rate vs. chromosome-size model** to stabilise per-chromosome estimates

- **Recommended sliding-window sizes** for target LD thresholds ($\rho$)

- **Built-in plotting** of decay summaries, per-chromosome fits, and window recommendations

- **Two-stage LD complexity reduction** to less redundant representatives -- the same clustering feeds either a pruned marker set or, optionally, one consensus genotype per block (an eMLG)

- **Outlier regions from Stage-1 units** (`ld_outlier_test()` / `ld_outlier_perm()`) -- test Stage-1 clusters using consensus-dosage or Simes-combined marker p-values, assemble significant units into reported regions, and compare discovery counts with a study-appropriate phenotype null. A separate region-relocation test asks whether the reported intervals overlap external annotations more often than their spans predict.

- **Diagnostic plotting** comparing raw vs. consolidated clusters chromosome-by-chromosome

- **Best single-SNP proxy per block** (`eMLG_best_snp()`) -- the consensus-optimal member SNP with missing calls filled from the consensus, for signals better kept at SNP resolution

------------------------------------------------------------------------

## Installation

```
# install from GitHub
pak::pak("genomeinsights/LDscnR")

# or use devtools
devtools::install_github(repo = "genomeinsights/LDscnR")
```

------------------------------------------------------------------------

## Quick example

```
library(LDscnR)
library(data.table)

data("sim_ex")

map <- sim_ex$map
GTs <- sim_ex$GTs
colnames(GTs) <- map$marker  # the bundled matrix and map use different marker labels

# Create a GDS object from the genotype matrix
gds_path <- tempfile(fileext = ".gds")
gds <- create_gds_from_geno(GTs, map, gds_path)

# Estimate chromosome-wise LD decay
ld_decay <- compute_LD_decay(gds, keep_el = TRUE, seed = 1)

# Inspect the fitted decay parameters and window recommendations
ld_decay

# Visualise the results
plot(ld_decay, type = "summary")
plot(ld_decay, type = "recommendation", rho = 0.99)
plot(ld_decay, type = "chr", chr = "Chr2")
```

------------------------------------------------------------------------

## Conceptual overview

LD decay is modelled per chromosome as:

$$r^2(d) = b + \frac{c - b}{1 + a\,d}$$

where

- $a$ controls the rate of decay,
- $b$ is background LD (long-distance baseline),
- $c$ is short-range LD, and
- $d$ is physical distance (bp).

The workflow proceeds in the following steps:

### 1. Background LD estimation

Background LD ($b$) is estimated from inter-chromosomal SNP pairs, giving the long-distance baseline against which decay is measured.

### 2. Chromosome-wise decay fitting

Decay parameters are estimated per chromosome using sliding windows, then robustly aggregated across windows. A model relating decay rate $a$ to chromosome size is fitted so that per-chromosome estimates can be stabilised and extrapolated across heterogeneous genomic architectures.

### 3. Window-size recommendations

For each target LD threshold $\rho$, `LDscnR` derives a recommended sliding-window size (in SNP units). The helpers `d_from_rho()` and `ld_from_rho()` convert a relative threshold $\rho$ into a physical distance and an expected $r^2$, respectively.

------------------------------------------------------------------------

## Output

`compute_LD_decay()` returns an object of class `"ld_decay"` containing:

- `by_chr` — per-chromosome decay fits (and optional LD edge lists when `keep_el = TRUE`),
- `decay_sum` — chromosome-wise decay parameters and derived quantities,
- `decay_model` — model linking decay rate to chromosome size,
- `recommendation` — suggested window sizes per $\rho$ threshold,
- `params` — parameters used in the computation.

`print()` and `plot()` methods are provided for `"ld_decay"` objects.

------------------------------------------------------------------------

## From LD decay to pruned markers and eMLGs

Once you have an `ld_decay` object (built with `keep_el = TRUE`), the workflow branches into two related but distinct end uses -- pruning, and eMLG generation -- built on the **same** clustering step.

### 1. Per-marker local LD support

`compute_ld_w()` summarises, for each marker, how much LD support it has from nearby markers within the physical window implied by a relative threshold $\rho$. It accepts a vector of $\rho$ values and computes all of them in one pass per chromosome -- each chromosome's edge list is read/symmetrised once and reused, not re-read once per threshold:

```
ld_w <- compute_ld_w(ld_decay, rho = c(0.90, 0.95, 0.99), cores = 4)
map[, ld_w_095 := ld_w[, "rho_0.95"]]
```

### 2. Stage 1 -- reduce marker redundancy

`ld_complexity_reduction()` clusters markers within each chromosome (connected components, then complete-linkage refinement within each) and picks one representative marker per cluster. This single call is the shared starting point for **both** downstream uses:

```
stage1 <- ld_complexity_reduction(map = map, LD_decay = ld_decay, rho = 0.5, cores = 4)
```

- **Pruning only** (a GRM, EMMAX's `K`, BayPass's `OMEGA`): `stage1$pruned` is already a ready-to-use character vector of representative markers -- nothing further is needed.
- **Block consensus genotypes** (eMLGs, for long-range LD/Ohta's D or other analyses where a single representative SNP would throw away most of a correlated block's information): feed `stage1` into `ld_prune_and_eMLG()` (below). It's the same Stage-1 clusters either way -- the choice is only what you do with their output.

**Why complete linkage on real pairwise values?** Standard sliding-window LD pruning -- including SNPRelate's own `snpgdsLDpruning()` (the GDS backend `LDscnR` itself builds on) and PLINK's `--indep-pairwise` -- makes a single greedy, order-dependent pass along the chromosome: starting from a (by default, semi-random) position, it admits a candidate marker only if it doesn't exceed the LD threshold against any marker already retained within the window, then moves on without ever revisiting that decision as more of the chromosome is seen. Two consequences, the first checked directly on the bundled `sim_ex` data (`snpgdsLDpruning(..., ld.threshold = 0.2, method = "corr")`, three `set.seed()` values): the retained set is not deterministic -- three seeds retained 93, 96, and 95 markers respectively, with only \~39% overlap between two of those runs on *identical* input -- and the output is just a flat list of retained marker IDs, with no record of which markers were considered redundant with which survivor, so nothing downstream (an eMLG-style consensus genotype, or even just knowing a pruned marker's effective "weight") can be built from it. Neither property is a bug in `snpgdsLDpruning()` -- both are direct, expected consequences of a single greedy walk that never double-checks its own earlier decisions, which is exactly the design LDscnR's Stage 1 avoids.

Stage 1 avoids an order-dependent greedy choice: within each connected component, it refines markers by complete linkage using their pairwise r² values. This prevents a chain of pairwise connections from joining markers that are weakly correlated at opposite ends of the chain. The representative has the highest median r² to other members of its cluster. Cluster membership is retained for downstream consensus genotypes. Separate final clusters can still be correlated, however; Stage 1 reduces redundant tests without certifying their independence.

(Stage 1 alone still has one honest caveat: its edge list comes from a sliding window, so a genuinely one-block region can fragment into adjacent clusters whose representatives were simply never directly compared at all. Stage 2 below exists specifically to close that gap.)

### 3. Stage 2 (optional) -- consolidate and summarise as eMLGs

This optional full-genome use of Stage 2 supports pruning and eMLG analyses.
In the manuscript's outlier scans, Stage-1 clusters remain the test units:
Stage 2 is applied only to significant units, after testing, to assemble
reported regions.

`ld_prune_and_eMLG()` closes the sliding-window gap above, but only for the clusters that need it: those flagged by high local LD support (`ld_w_col`/`ld_w_threshold`) are re-compared directly from genotypes -- with no window restriction -- and consolidated via a distance-restricted, quality-gated dynamic cut. This produces a refined pruned marker set and an eMLG matrix from the same pass; unflagged clusters (usually the large majority) pass straight through unchanged. `distance_threshold` -- the max physical gap allowed within one mergeable, contiguous block -- defaults to a per-chromosome value derived from `rho` and `LD_decay` (`d_from_rho(a_pred, rho)`), reusing the same `rho` that `ld_w_col`'s naming already implies, rather than one fixed bp constant:

```
# Illustrative full-genome flagging threshold for this example; inspect the
# flagged fraction and sensitivity before choosing a value for new data.
result <- ld_prune_and_eMLG(
  GTs = GTs, stage1 = stage1, ld_w_col = "ld_w_095", ld_w_threshold = 0.05,
  LD_decay = ld_decay, rho = 0.95,
  score_threshold = 0.80, min_r2 = 0.2, cores = 4
)

pruned_markers <- result$pruned   # refined pruned marker set
eMLG           <- result$eMLG     # individuals x blocks consensus-genotype matrix
```

### 4. Visual diagnostic

`plot_pruning_comparison()` stacks Stage 1's raw clusters (top panel) against Stage 2's consolidated groups (bottom panel) for one chromosome, so fragmented blocks and their reunited counterparts can be compared directly. `ld_w_col`/`ld_w_threshold`/`min_n_loci_flag` default to `result$params` -- whatever `ld_prune_and_eMLG()` was actually called with above -- so the two panels can't silently end up comparing different flagging criteria; passing a value that disagrees with `result$params` warns:

```
plot_pruning_comparison(chr = "Chr3", pruned_stage1 = stage1, result = result, map = map)
```

### Speed tip: preliminary vs. final runs

`ld_prune_and_eMLG()`'s cost is dominated by an all-pairs correlation among the *flagged* clusters, which scales roughly quadratically with how many clusters get flagged (on real data: \~0.01s at 292 flagged clusters, \~31s at 15,000). `ld_w_threshold` is the lever that matters, and it should move in different directions depending on the run:

- **Preliminary/exploratory runs**: use a high `ld_w_threshold` to flag only the most obviously redundant clusters -- fast, good enough for a first look at cluster counts and eMLG behaviour.
- **Final run**: examine the fraction of clusters flagged and lower `ld_w_threshold` where broader consolidation is needed (optionally combined with `min_n_loci_flag`). The value `0.05` in this example is not a general cutoff. This optional full-genome consolidation is distinct from outlier-region assembly, which re-examines all significant units.
- **`compute_unflagged_eMLG = FALSE`** skips eMLG computation for the unflagged clusters entirely (usually the large majority) if you only need the pruned marker set, independent of `ld_w_threshold`.

### 5. Best single-SNP proxy per block (optional)

The eMLG consensus averages a block's markers -- ideal for one value per block, but it dilutes signal that is genuinely SNP-specific (differing even between strongly-linked markers). `eMLG_best_snp()` is for that case: for each block with a stored consensus it picks the member SNP most correlated with the consensus, and can return that SNP's genotype with only its *missing* calls filled from the consensus (observed calls are never overwritten, so SNP-level signal is preserved). The clustering's `representative` is chosen for cluster centrality, not to reproduce the consensus, so the two differ; the `score` (eMLG fidelity) column indicates when one SNP suffices (high score  consensus adds little) vs. when the consensus does real work (low score):

```
best <- eMLG_best_snp(result, GTs)

# per block: representative vs consensus-optimal best_marker, their |r| to the
# consensus, the fidelity score, and (with fill = TRUE) the fill counts
best$stats

# drop-in single-SNP alternative to result$eMLG: same shape, but each column is
# an observed SNP with its gaps filled from the consensus
best$geno
```

------------------------------------------------------------------------

## Outlier regions from Stage-1 units (`ld_outlier_test()`)

Everything above builds LD structure and reduces markers to Stage-1 clusters. This section
tests those clusters directly -- BH across cluster-level p-values, then assembles the
significant clusters into reported regions. It is the design behind every reported result in
the LDscnR manuscript: the simulation benchmark and both stickleback panels. See
`vignette("LDscnR_stage1_outlier_regions")` for the full, executable walkthrough this section
summarises, including the structure-aware permutation null and the separate annotation-overlap
check.

The snippets below show how the functions fit together; `K`, `y`,
`annotation`, `chrom_lengths` and the surrogate p-values must come from the
study being analysed.
For a self-contained example, run the primary vignette. The default
`size_floor = 8` is not a universal cutoff: the manuscript set empirical
floors by filtered marker density and examined their sensitivity.

It stays **engine-agnostic**: LDscnR does not fit your association model. You supply p-values
from whatever engine you like -- EMMAX, LFMM, a GLM, an $F_{ST}$ scan; the package also ships
a fast EMMAX implementation (`emmax_setup()` + `emmax_fast()`) purely for convenience, not
because anything below requires it. This does **not** include BayPass out of the box: BayPass
reports Bayes factors and XtX differentiation statistics, not p-values, so `ld_outlier_test()`
cannot take its output directly. BayPass fits into the earlier, clustering half of this
workflow instead -- Stage-1's LD-reduced marker groups (`stage1$pruned`, or an eMLG per group)
can serve as BayPass's marker set or feed its `OMEGA` covariance estimate -- not into this
p-value-testing step.

### 1. Build a per-unit variable to test

A Stage-1 cluster is not itself something an association model can take. `ld_unit_matrix()`
turns each cluster into one variable per individual:

```r
stage1 <- ld_complexity_reduction(map = map, LD_decay = ld_decay, rho = 0.5, cores = 4)

# consensus dosage -- one averaged genotype per unit, no cluster-size penalty; the arm used
# for every reported EMMAX-consensus result in the manuscript
um <- ld_unit_matrix(GTs, stage1, map, size_floor = 8L, repr = "consensus_dosage")
P     <- emmax_setup(um, K)          # K = your relationship matrix
p_obs <- emmax_fast(P, y)            # one p-value per unit
```

`repr` also accepts `"eMLG"` (`make_eMLGs()`'s own block consensus) and `"representative"`
(the cluster's most central marker). A fourth SNP-level choice, the member most correlated
with a block's consensus, is available via `ld_group_matrix()` instead -- not as a
`ld_unit_matrix()` option, because it operates on `ld_prune_and_eMLG()`'s groups, a
different (and possibly coarser) partition from this function's Stage-1 units; its output
cannot be substituted for a `ld_unit_matrix()` column when testing with `ld_outlier_test()`.
See `?ld_group_matrix`.

**Or skip `ld_unit_matrix()` entirely** and combine ordinary marker-wise p-values with Simes
instead -- the comparator arm in the manuscript, and the only option when an engine can't be
refit to a reduced matrix at all (e.g. LFMM, which estimates its latent factors from the full
marker set):

```r
Pm    <- emmax_setup(GTs, K)
p_obs <- emmax_fast(Pm, y)           # one p-value per marker
```

### 2. Test units and assemble regions

```r
test <- ld_outlier_test(
  stage1, map, p_obs,
  statistic  = "unit",                # "unit" for a matrix built above, "simes" for marker p-values
  size_floor = 8L, alpha = 0.05,
  assembly   = "stage2_discovered",   # re-examine only the significant clusters from genotypes
  GTs = GTs, LD_decay = ld_decay
)
test                                  # units tested / significant -> regions assembled
```

BH is applied to `test$units` -- one p-value per tested cluster, never per marker.
`assembly = "stage2_discovered"` then re-runs `ld_prune_and_eMLG()` over only the significant
clusters, directly from genotypes (the same Stage-2 step used for pruning/eMLGs above), so a
reported region can only span discovered signal. `assembly = "physical"`, a plain merge of
significant clusters within `gap`, is retained only as a near-free check against that
motivated rule, not as an equally preferred alternative.

### 3. Check the discovery burden and external overlap

These checks answer different questions; neither validates the other.

**Structure-aware phenotype null** -- how many discoveries the same pipeline
makes when the study's relevant population or spatial structure is preserved:

```r
# p_perm is a units-by-surrogates matrix from re-running the association
# model on study-appropriate surrogate phenotypes.
null <- ld_outlier_perm(test, stage1, map, p_perm, GTs = GTs, LD_decay = ld_decay,
                        level = "units")
null                                  # observed vs. surrogate discovery counts, one-sided p
```

Constructing `p_perm` is a study-design decision: which unit is exchangeable,
which structure is held fixed, and whether the surrogates reproduce the
confounding pattern of concern. A null that breaks that structure can give
false reassurance. The function repeats the same LDscnR testing and reporting
steps for every surrogate. The primary vignette shows a concrete example.

**Region relocation** -- whether the resulting regions overlap an external annotation (EcoPeaks, a QTL
panel, a gene list) more than a span-preserving null predicts:

```r
rot <- ld_region_rotation(test$regions, annotation, chrom_lengths, scheme = "within",
                          n_rotations = 10000)
rot                                   # observed overlaps, fold enrichment, rotation p
```

Key arguments:

| argument | what it controls |
|---|---|
| `size_floor` | minimum markers per tested Stage-1 unit; choose for the dataset and examine sensitivity |
| `statistic` | `"unit"` (pre-built matrix from `ld_unit_matrix()`) or `"simes"` (combine marker p-values per unit) |
| `assembly` | `"stage2_discovered"` (genotype-based, the inferential path) or `"physical"` (gap merge, a sanity check) |
| `alpha` | BH level for unit significance; `ld_outlier_perm()` reuses `test$params`. Region relocation is a separate annotation-overlap analysis. |
| `level` (`ld_outlier_perm()`) | count `"units"` (cheap -- skips region assembly entirely) or `"regions"` |
| `scheme` (`ld_region_rotation()`) | `"within"` (preserve each region's chromosome and randomly relocate it) or `"genome"` (allow reassignment among chromosomes) |

### Earlier method retained in the source

Earlier C-score code remains in the package for reproducibility, but it is not
the current outlier workflow and has no installed tutorial. Its former
vignettes are retained in [`vignettes-archive/`](vignettes-archive/) for
historical reference. The Stage-1 workflow above is the one used in the
current LDscnR manuscript.

------------------------------------------------------------------------

## Documentation

```
vignette("LDscnR_stage1_outlier_regions")   # PRIMARY: the current Stage-1-cluster outlier pipeline, executable end to end
vignette("LDscnR_quick_introduction")       # LD decay, ld_w, pruning, eMLGs
vignette("LDscnR_complexity_reduction")     # LD decay and complexity reduction on real stickleback data
```

`vignette("LDscnR_stage1_outlier_regions")` is the primary outlier-analysis vignette: it runs
the current pipeline (`ld_unit_matrix()` / `ld_outlier_test()` / `ld_outlier_perm()` /
`ld_region_rotation()`) end to end on the bundled `stickleback` panel,
executable as written. `ld_outlier_scan()` wraps these steps when a single
call is convenient; the component functions make each check easier to inspect.
The archived C-score vignettes are not installed or
recommended for current analyses.

------------------------------------------------------------------------

## Notes

- LD-decay estimation is performed in two steps. First on subsets of SNPs per chromosome using a large sliding window (e.g. 1000 SNPs). Based on this, a smaller sliding window is chosen such that a target fraction of the decay curve is covered, reducing the number of pairwise comparisons in the full run.

- Large data sets can generate substantial intermediate objects (LD edge lists). Edge lists can be written to files instead of held in RAM via the `el_data_folder` argument.

- Parallelization is supported via the `cores` argument (`mclapply`).

- For LD-pruning/eMLG generation (`ld_prune_and_eMLG()`), raise `ld_w_threshold` to flag fewer clusters during exploratory runs -- see the " From LD decay to pruned markers and eMLGs" section above.

- `compute_LD_decay()` takes a `seed`. It subsamples the background, thins per chromosome and samples pairs within strata, so **an unseeded refit moves every quantity derived from it** -- `ld_w`, the pruned marker set, the clustering. Set `seed` on any run whose output will be compared with another's.

- `ld_outlier_test()`/`ld_outlier_perm()` report the p-values you give them. They do not check that `p_perm` came from a valid null, because they cannot: whether a permutation scheme is admissible depends on the design, not on the numbers. A near-1.0 genomic inflation factor is a **body** statistic and is not evidence that a cluster- or region-level statistic is usable with FDR control -- FDR lives in the tail. Measure the tail.

------------------------------------------------------------------------

## Dependencies

- `data.table`

- `SNPRelate`

- `igraph` (Stage-1 complexity reduction)

- `ggplot2`, `patchwork` (`plot_pruning_comparison()`)

------------------------------------------------------------------------

## License

MIT

------------------------------------------------------------------------

## Author

Petri Kemppainen - petri\@genomeinsights.fi
