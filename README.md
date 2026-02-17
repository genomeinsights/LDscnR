---
editor_options: 
  markdown: 
    wrap: 72
---

# Statistical Assumptions and Input Requirements

## Overview

LDscnR performs LD-aware scaling and quantile correction of genome scan
statistics. The F-distribution provides a convenient and flexible
reference family for one-degree-of-freedom tests commonly used in genome
scans. The method is, however, agnostic to the specific test statistic
used, and only assumes:

1.  **Monotonicity.** Larger values of the statistic $S_i$ correspond to
    stronger evidence against the null hypothesis.

2.  **Known null distribution.** Under the null hypothesis, $S_i$ can be
    mapped to a theoretical distribution with known quantiles.

3.  **Independence of null distribution from LD correction.** The
    theoretical null distribution is defined independently of the
    LD-weighting procedure.

Thus, the LD-scaling algorithm does not depend on the interpretation of
$S_i$ itself, but rather on its null quantile structure.

------------------------------------------------------------------------

## Transformation to the F-scale

LDscnR expects statistics expressed as:

$$
  F_i \sim F(df_1, df_2)
$$

under the null hypothesis.

Many commonly used genome scan statistics satisfy this condition either
directly or via transformation:

| Original statistic   | Transformation to F-scale             |
|----------------------|---------------------------------------|
| p-value              | $F = F^{-1}(1 - p)$                   |
| $\chi^2_1$ statistic | $\chi^2_1 = F(1, \infty)$             |
| $t$-statistic        | $t^2 \sim F(1, df)$                   |
| $z$-score            | $z^2 \sim \chi^2_1 \sim F(1, \infty)$ |
| Likelihood ratio     | asymptotically $\chi^2 \rightarrow F$ |

Thus, for most genome scan statistics, the transformation route is:

$$
  \text{statistic} \rightarrow p\text{-value} \rightarrow F\text{-statistic}.
$$

Once transformed, the statistic can be passed to `ld_scan()` or
`LDscn-pipeline()`.

## LD-decay model with short-range plateau

In many empirical datasets, linkage disequilibrium (LD) remains close to
its maximum value over short genomic distances before beginning to
decline. To accommodate this behavior, LDscnR extends the classical
hyperbolic LD-decay model by introducing a short-range plateau parameter
$d_0$:

$$
r^2(d) = b + \frac{c - b}{1 + a \max(d - d_0, 0)}.
$$

Here,

-   $b$ denotes the background LD level,
-   $c$ is the short-range LD plateau (bounded by 1),
-   $a$ controls the rate of decay,
-   $d_0$ represents the distance up to which LD remains approximately
    constant.

The term $\max(d - d_0, 0)$ ensures that LD remains constant at level
$c$ for distances $d \le d_0$, and then decays smoothly according to
rate parameter $a$ beyond this point. This formulation avoids
discontinuities that could hinder numerical optimization, while allowing
greater flexibility in modeling short-range LD patterns.

When SNP density is high, $d_0$ captures the observed short-distance LD
plateau explicitly and stabilizes estimation of $c$. In sparse marker
datasets, $d_0$ typically converges to zero, reducing the model to the
standard hyperbolic decay formulation.

## Edge-effect correction for local LD weights (`ld_w`)

Local LD weights (`ld_w`) are defined as the median $r^2$ between a
focal SNP and all neighboring SNPs within a distance threshold $d_\rho$.
However, SNPs located near chromosome boundaries or near large physical
gaps naturally have fewer neighbors on one side. Without correction,
this asymmetry leads to systematically downward-biased LD weights at
genomic edges.

To remove this bias, LDscnR applies a symmetric edge correction.

For each focal SNP, we:

1.  Compute distances to neighboring SNPs on both the left and right.
2.  Determine the maximum observed distance on each side: $$
    d_L^{\max}, \quad d_R^{\max}.
    $$
3.  Define the effective symmetric radius: $$
    S = \min(d_L^{\max}, d_R^{\max}).
    $$

This ensures that only the symmetric region around the SNP contributes
fully to the LD estimate.

If one side extends further than the other (the “long side”), pairs
beyond the symmetric radius $S$ on the long side are duplicated. This
duplication balances the contribution of left and right neighborhoods
when computing the median LD:

$$
ld_w = \operatorname{median}(r^2).
$$

This approach:

-   Uses *all available pairwise information* within the chosen window.
-   Maintains symmetry in the effective window.
-   Avoids artificial down-weighting of SNPs near chromosome ends.
-   Does not introduce smoothing or additional model assumptions.

In dense genomic datasets, the correction has minimal impact in interior
regions but substantially improves stability near chromosome boundaries.

## LDscnR Pipeline Overview

`ld_pipeline()` provides a streamlined workflow for LD-aware genome
scans. Starting from a genotype matrix and SNP annotation table, it
performs:

1.  **LD structure computation**\
    Efficient calculation of pairwise LD within chromosomes.

2.  **LD-decay estimation**\
    Fits a flexible LD-decay model:

    $$
    r^2(d) = b + \frac{c - b}{1 + a \max(d - d_0, 0)},
    $$

    where:

    -   $b$ = background LD\
    -   $c$ = short-range plateau\
    -   $a$ = decay rate\
    -   $d_0$ = optional plateau distance

3.  **LD-scaled scan (F′ statistic)**\
    Adjusts test statistics using local LD weights to correct for
    LD-induced inflation while preserving selection-consistent signals.

4.  **Outlier region (OR) detection**\
    Identifies clusters of significant SNPs based on LD connectivity and
    distance thresholds expressed relative to the empirical decay curve.

5.  **Consistency scoring (C-score)**\
    Aggregates evidence across multiple LD-window parameterisations,
    yielding a stability-based measure of signal robustness.

------------------------------------------------------------------------

### Minimal example

``` r
LDscn <- LDscn_pipeline(
  geno   = GTs,
  map    = map,
  F_cols = c("emx_F", "lfmm_F"),
  n_inds = nrow(GTs)
)

plot(LDscn$scan, method = "lfmm_F")
```

------------------------------------------------------------------------

### Why use the pipeline?

-   Fully LD-aware genome scan workflow
-   Relative thresholding across heterogeneous LD landscapes
-   Plateau-aware LD-decay modelling
-   Robust signal ranking via consistency scoring
-   Modular — components can be reused independently

------------------------------------------------------------------------

### Memory note

The LD structure object can be large. If needed, save and remove it
after computation:

``` r
saveRDS(pipe$ld_struct, "ld_structure.rds")
pipe$ld_struct <- NULL
gc()
```

For large datasets, computing the LD structure once and reusing it
across analyses is recommended.
