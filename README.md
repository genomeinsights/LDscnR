# LDscnR

**LD-aware genome scan and outlier-region detection**

`LDscnR` provides tools for analysing genome scan results while accounting for linkage disequilibrium (LD). The package combines LD decay estimation, LD-scaled statistics, and repeated outlier-region detection to identify robust signals of selection.

------------------------------------------------------------------------

## 🚀 Key features

-   📉 **LD decay estimation** (chromosome-specific)

-   ⚖️ **LD-scaled statistics** (F')

-   🧬 **Outlier-region clustering** using LD and distance thresholds

-   📊 **Consistency scores (C)** summarizing robustness

-   📈 **Manhattan-style visualization**

-   🧪 **Simulation benchmarking tools** (precision–recall, AUC)

------------------------------------------------------------------------

## 📦 Installation

```         
# install from GitHub
devtools::install_github("genomeinsights/LDscnR")
```

------------------------------------------------------------------------

## 📊 Quick example

```         
library(LDscnR)
library(data.table)

data("sim_ex")

map <- sim_ex$map

# Create GDS object
gds_path <- tempfile(fileext = ".gds")
gds <- create_gds_from_geno(sim_ex$GTs, map, gds_path)

# Estimate LD decay
ld_decay <- compute_LD_decay(gds, keep_el = TRUE)

# Generate outlier-region draws
draws <- ld_rho_draws(
  gds,
  ld_decay = ld_decay,
  F_vals = map[, .(lfmm_F, emx_F)],
  n_draws = 50
)

# Plot manhattan based on the consitency score (C) and colored by outlier region
plot_manhattan(
  map,
  gds,
  ld_decay,
  draws
)
```

------------------------------------------------------------------------

## 🧠 Conceptual overview

`LDscnR` is built around three core steps:

### 1. LD decay estimation

LD decay is estimated per chromosome to define a biologically meaningful scale for clustering SNPs into regions. Decay rate is regressed against chromosome size and predicted rates are used to reduce chromosome level variation caused by large haplotype blocks etc (e.g. inversions).

### 2. Repeated outlier-region detection

Based on F-outlier statistics e.g. from genotype-environment association analyses, scaled by local LD ($ld_w$ the median LD between a focal SNPs and other snips with window size $w$), outlier regions are detected repeatedly across a range of:

-   LD window sizes

-   distance thresholds for outlier clustering

-   significance thresholds

This avoids reliance on a single arbitrary cutoff.

### 3. Consistency scoring

Each SNP receives a **consistency score (C)** reflecting how often it appears in detected regions across parameter space.

------------------------------------------------------------------------

## 📈 Output interpretation

-   **High C-score SNPs** → robust signals across parameter space

-   **Outlier regions (ORs)** → clusters of SNPs supported by LD

-   **Joint mode** → combines evidence across methods

-   **Per-method mode** → method-specific detection

------------------------------------------------------------------------

## 🧪 Simulation support

The package includes tools to evaluate performance in simulated data:

-   Precision–recall summaries (`get_PR`)

-   AUC of cumulative PR curves (`get_AUC_OR`)

-   True/false positive OR classification

------------------------------------------------------------------------

## 📚 Documentation

See the vignette for a full workflow:

```         
vignette("LDscnR_quick_introduction")
```

------------------------------------------------------------------------

## ⚠️ Notes

-   LD-decay estimation is performed in two steps. First based on subsets of SNPs from chromosomes, but large sliding windows (1000 bps). Based on this, a new (smaller) sliding window is determined such that 99% of the decay curve is covered to reduce the number of pairwise comparisons for subsequent downstream analyses.

-   However, large data sets can still generate substantial intermediate objects (LD edge lists, OR draws)

-   For heavy workflows, consider saving intermediate results using `saveRDS()`

-   Parallelization is supported via the `cores` argument (`mclapply`)

------------------------------------------------------------------------

## 🔧 Dependencies

-   `data.table`

-   `ggplot2`

-   `patchwork`

-   `SNPRelate`

-   `future.apply`

------------------------------------------------------------------------

## 📜 License

MIT

------------------------------------------------------------------------

## 👨‍🔬 Author

Petri Kemppainen - petri\@genomeinsights.fi

------------------------------------------------------------------------

## 💡 Citation

If you use `LDscnR`, please cite: <https://www.biorxiv.org/content/10.64898/2026.01.19.700334v1>
