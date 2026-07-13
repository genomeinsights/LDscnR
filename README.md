# LDscnR

**Chromosome-wise LD-decay estimation**

`LDscnR` provides tools for estimating linkage disequilibrium (LD) decay from genotype data. It fits chromosome-specific decay curves, relates decay rate to chromosome size, and derives recommended LD window sizes for user-defined relative LD thresholds ($\rho$). A tutorial PDF can be found in vignettes.

------------------------------------------------------------------------

## 🚀 Key features

- 📉 **Chromosome-specific LD-decay estimation**

- 🧮 **Background LD estimation** from inter-chromosomal SNP pairs

- 📐 **Decay-rate vs. chromosome-size model** to stabilise per-chromosome estimates

- 🎯 **Recommended sliding-window sizes** for target LD thresholds ($\rho$)

- 📊 **Built-in plotting** of decay summaries, per-chromosome fits, and window recommendations

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

------------------------------------------------------------------------

## 🔧 Dependencies

- `data.table`

- `SNPRelate`

------------------------------------------------------------------------

## 📜 License

MIT

------------------------------------------------------------------------

## 👨‍🔬 Author

Petri Kemppainen - petri\@genomeinsights.fi

------------------------------------------------------------------------

## 💡 Citation

If you use `LDscnR`, please cite: <https://www.biorxiv.org/content/10.64898/2026.01.19.700334v1>
