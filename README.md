# LDscnR

**Chromosome-wise LD-decay estimation, LD-based marker pruning, and eMLG generation**

`LDscnR` provides tools for estimating linkage disequilibrium (LD) decay from genotype data, and for using that decay model to reduce a marker set to LD-independent representatives. The same clustering step can be used either to produce a pruned marker set (for a kinship/relatedness matrix, EMMAX's `K`, BayPass's `OMEGA`) or, from that same clustering, one consensus genotype per LD block (an "eMLG" -- expected multi-locus genotype -- for block-level analyses such as long-range LD or Ohta's D statistics). A tutorial PDF can be found in vignettes.

------------------------------------------------------------------------

## 🚀 Key features

- 📉 **Chromosome-specific LD-decay estimation**

- 🧮 **Background LD estimation** from inter-chromosomal SNP pairs

- 📐 **Decay-rate vs. chromosome-size model** to stabilise per-chromosome estimates

- 🎯 **Recommended sliding-window sizes** for target LD thresholds ($\rho$)

- 📊 **Built-in plotting** of decay summaries, per-chromosome fits, and window recommendations

- 🧩 **Two-stage LD complexity reduction** to LD-independent representatives -- the same clustering feeds either a pruned marker set or, optionally, one consensus genotype per block (an eMLG)

- 🔎 **Diagnostic plotting** comparing raw vs. consolidated clusters chromosome-by-chromosome

------------------------------------------------------------------------

## 📦 Installation

```         
# install from GitHub
pak::pak("genomeinsights/LDscnR") 

# or (if you don't have a GitHub account)
devtools::install_github(repo = "genomeinsights/LDscnR")
```

------------------------------------------------------------------------

## 📊 Quick example

```         
library(LDscnR)
library(data.table)

data("sim_ex")

map <- sim_ex$map

# Create a GDS object from the genotype matrix
gds_path <- tempfile(fileext = ".gds")
gds <- create_gds_from_geno(sim_ex$GTs, map, gds_path)

# Estimate chromosome-wise LD decay
ld_decay <- compute_LD_decay(gds, keep_el = TRUE)

# Inspect the fitted decay parameters and window recommendations
ld_decay

# Visualise the results
plot(ld_decay, type = "summary")
plot(ld_decay, type = "recommendation", rho = 0.99)
plot(ld_decay, type = "chr", chr = "Chr2")
```

------------------------------------------------------------------------

## 🧠 Conceptual overview

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

## 📈 Output

`compute_LD_decay()` returns an object of class `"ld_decay"` containing:

- `by_chr` — per-chromosome decay fits (and optional LD edge lists when `keep_el = TRUE`),
- `decay_sum` — chromosome-wise decay parameters and derived quantities,
- `decay_model` — model linking decay rate to chromosome size,
- `recommendation` — suggested window sizes per $\rho$ threshold,
- `params` — parameters used in the computation.

`print()` and `plot()` methods are provided for `"ld_decay"` objects.

------------------------------------------------------------------------

## 🧩 From LD decay to pruned markers and eMLGs

Once you have an `ld_decay` object (built with `keep_el = TRUE`), the workflow branches into two related but distinct end uses -- pruning, and eMLG generation -- built on the **same** clustering step.

### 1. Per-marker local LD support

`compute_ld_w()` summarises, for each marker, how much LD support it has from nearby markers within the physical window implied by a relative threshold $\rho$. It accepts a vector of $\rho$ values and computes all of them in one pass per chromosome -- each chromosome's edge list is read/symmetrised once and reused, not re-read once per threshold:

```         
ld_w <- compute_ld_w(ld_decay, rho = c(0.90, 0.95, 0.99), cores = 4)
map[, ld_w_095 := ld_w[, "rho_0.95"]]
```

### 2. Stage 1 -- reduce markers to LD-independent representatives

`ld_complexity_reduction()` clusters markers within each chromosome (connected components, then complete-linkage refinement within each) and picks one representative marker per cluster. This single call is the shared starting point for **both** downstream uses:

```         
stage1 <- ld_complexity_reduction(map = map, LD_decay = ld_decay, rho = 0.5, cores = 4)
```

- **Pruning only** (a GRM, EMMAX's `K`, BayPass's `OMEGA`): `stage1$pruned` is already a ready-to-use character vector of representative markers -- nothing further is needed.
- **Block consensus genotypes** (eMLGs, for long-range LD/Ohta's D or other analyses where a single representative SNP would throw away most of a correlated block's information): feed `stage1` into `ld_prune_and_eMLG()` (below). It's the same Stage-1 clusters either way -- the choice is only what you do with their output.

**Why complete linkage on real pairwise values?** Standard sliding-window LD pruning -- including SNPRelate's own `snpgdsLDpruning()` (the GDS backend `LDscnR` itself builds on) and PLINK's `--indep-pairwise` -- makes a single greedy, order-dependent pass along the chromosome: starting from a (by default, semi-random) position, it admits a candidate marker only if it doesn't exceed the LD threshold against any marker already retained within the window, then moves on without ever revisiting that decision as more of the chromosome is seen. Two consequences, the first checked directly on the bundled `sim_ex` data (`snpgdsLDpruning(..., ld.threshold = 0.2, method = "corr")`, three `set.seed()` values): the retained set is not deterministic -- three seeds retained 93, 96, and 95 markers respectively, with only \~39% overlap between two of those runs on *identical* input -- and the output is just a flat list of retained marker IDs, with no record of which markers were considered redundant with which survivor, so nothing downstream (an eMLG-style consensus genotype, or even just knowing a pruned marker's effective "weight") can be built from it. Neither property is a bug in `snpgdsLDpruning()` -- both are direct, expected consequences of a single greedy walk that never double-checks its own earlier decisions, which is exactly the design LDscnR's Stage 1 avoids.

Stage 1 avoids this by construction: within each connected component, `hclust()` clusters on the real pairwise r² sub-matrix using complete linkage, so a marker only joins a cluster if it clears the threshold against *every* existing member, not just its nearest neighbour -- a deterministic result with no dependence on traversal order or starting position (checked directly: identical `stage1$pruned` across different seeds on the same input). That also directly rules out single linkage's "chaining" failure (A uncorrelated with C, but both pulled into one group via an intermediate B correlated with each) -- the same failure mode that undermines greedy sequential pruning, since a candidate is never re-checked once an earlier, unrelated decision has already shaped the retained set. The representative marker for each cluster is chosen the same principled way: highest median r² to the rest of its own cluster, not whichever marker the walk happened to reach first. The result is markers in *different* final clusters are, by construction, never verified as mutually redundant with each other -- about as close to "each retained marker is a genuinely independent test unit" as a threshold-based method gets. The full cluster membership is retained too (not just the representative), which is what makes the eMLG option (below) possible at all.

(Stage 1 alone still has one honest caveat: its edge list comes from a sliding window, so a genuinely one-block region can fragment into adjacent clusters whose representatives were simply never directly compared at all. Stage 2 below exists specifically to close that gap.)

### 3. Stage 2 (optional) -- consolidate and summarise as eMLGs

`ld_prune_and_eMLG()` closes the sliding-window gap above, but only for the clusters that need it: those flagged by high local LD support (`ld_w_col`/`ld_w_threshold`) are re-compared directly from genotypes -- with no window restriction -- and consolidated via a distance-restricted, quality-gated dynamic cut. This produces a refined pruned marker set and an eMLG matrix from the same pass; unflagged clusters (usually the large majority) pass straight through unchanged. `distance_threshold` -- the max physical gap allowed within one mergeable, contiguous block -- defaults to a per-chromosome value derived from `rho` and `LD_decay` (`d_from_rho(a_pred, rho)`), reusing the same `rho` that `ld_w_col`'s naming already implies, rather than one fixed bp constant:

```         
# ld_w_threshold = 0.05 here is a "final run" value (see the speed tip
# below) -- flag more generously so that as many genuine mergers as
# possible actually happen; raise it for a faster preliminary look
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

### ⚡ Speed tip: preliminary vs. final runs

`ld_prune_and_eMLG()`'s cost is dominated by an all-pairs correlation among the *flagged* clusters, which scales roughly quadratically with how many clusters get flagged (on real data: \~0.01s at 292 flagged clusters, \~31s at 15,000). `ld_w_threshold` is the lever that matters, and it should move in different directions depending on the run:

- **Preliminary/exploratory runs**: use a high `ld_w_threshold` to flag only the most obviously redundant clusters -- fast, good enough for a first look at cluster counts and eMLG behaviour.
- **Final run**: lower `ld_w_threshold` toward `~0.05` (optionally combined with `min_n_loci_flag`, to also pull in large-but-low-`ld_w` clusters) so that as many genuine mergers as possible actually happen -- slower, but this only needs to be run once.
- **`compute_unflagged_eMLG = FALSE`** skips eMLG computation for the unflagged clusters entirely (usually the large majority) if you only need the pruned marker set, independent of `ld_w_threshold`.

------------------------------------------------------------------------

## 📚 Documentation

See the vignette for a full workflow:

```         
vignette("LDscnR_quick_introduction")
```

------------------------------------------------------------------------

## ⚠️ Notes

- LD-decay estimation is performed in two steps. First on subsets of SNPs per chromosome using a large sliding window (e.g. 1000 SNPs). Based on this, a smaller sliding window is chosen such that a target fraction of the decay curve is covered, reducing the number of pairwise comparisons in the full run.

- Large data sets can generate substantial intermediate objects (LD edge lists). Edge lists can be written to files instead of held in RAM via the `el_data_folder` argument.

- Parallelization is supported via the `cores` argument (`mclapply`).

- For LD-pruning/eMLG generation (`ld_prune_and_eMLG()`), raise `ld_w_threshold` to flag fewer clusters during exploratory runs -- see the "🧩 From LD decay to pruned markers and eMLGs" section above.

------------------------------------------------------------------------

## 🔧 Dependencies

- `data.table`

- `SNPRelate`

- `igraph` (Stage-1 complexity reduction)

- `ggplot2`, `patchwork` (`plot_pruning_comparison()`)

------------------------------------------------------------------------

## 📜 License

MIT

------------------------------------------------------------------------

## 👨‍🔬 Author

Petri Kemppainen - petri\@genomeinsights.fi

------------------------------------------------------------------------

## 💡 Citation

If you use `LDscnR`, please cite: <https://www.biorxiv.org/content/10.64898/2026.01.19.700334v1>
