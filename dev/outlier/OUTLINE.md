# Outline: LD-aware outlier-analysis vignette (draft)

**Title:** LD-aware outlier analysis on stickleback data
**Data:** bundled `stickleback` (Chr1 inversion, Chr4 *Eda*, Chr6 neutral).
Chr1 + Chr4 are established marine–freshwater loci = **positive controls**;
Chr6 = **negative control**. So the whole vignette doubles as a control test.
**Builds on:** the complexity-reduction vignette (recomputes the same pipeline).
**No new deps:** Zissou1 colours hardcoded (no `wesanderson` needed).

## 0. Overview + setup
- One paragraph: from a pruned set + eMLG (complexity-reduction vignette) to
  LD-aware association. State the positive/negative controls up front.
- Setup chunk: pipeline (ld_decay → ld_w → stage1 → result), ecotype phenotype,
  **LOCO** GRMs from the pruned representatives.

## 1. The cluster-consensus scan (Li et al. 2018) — and where it falls short
- Method: test each LD cluster's **eMLG consensus** (the PC1 analog) with
  `emmax()` (LOCO), FDR over clusters.
- **Fig D** (LOCO Manhattan): consensus-significant clusters coloured (Chr1
  inversion); *Eda* appears as black triangles = single-SNP-significant but
  **not** cluster-flagged.
- Result: finds the inversion, **misses Chr4 *Eda*** (consensus q just above FDR).
- Beat: averaging into a consensus + few-test FDR trades away localized signals.
  Motivates keeping SNP resolution but prioritizing by local LD.

## 2. Local-LD prioritization: `ld_w`
- Premise: genuine signal sits in locally-coherent LD; isolated significant SNPs
  with low `ld_w` are likely false positives.
- **Fig A / A2** (single-SNP Manhattan, colour+alpha = `ld_w`, Zissou): the peaks
  (inversion, *Eda*) are the high-`ld_w` structure; Chr6 fades to nothing.
- **Bimodality** (table): single-SNP-only-significant hits split into
  in-`ld_w`-cluster (median `ld_w` 0.61, real) vs isolated (0.01, noise).
  The quantitative anchor.

## 3. A data-driven `ld_w` threshold (independent filtering)
- `ld_w` is a filter statistic ~independent of association p under the null
  (LOCO). Choose the cutoff maximizing FDR discoveries (Bourgon et al. 2010).
- **Fig E** (rejection-maximization curve): single peak → **q\* = 0.87**.
- Prioritized scan (q ≥ q\*): 1445 discoveries vs 56 all-SNP, and it **catches
  *Eda*** that the consensus-cluster scan missed.
- Caveats (state plainly): independence assumption; rejection count is
  inversion-dominated (a threshold-*selection* device, not an effect size).

## 4. Validating the premise: F vs `ld_w` quantile
- **Fig B** (F vs genome-wide `ld_w` quantile per chr; rolling median; threshold
  line): Chr1/Chr4 (positive controls) rise past q\*; Chr6 (negative control) is
  flat with essentially nothing beyond the threshold.
- Beat: signal concentrates above the data-driven threshold on signal
  chromosomes, and q\* lands at the F-rise elbow — data and biology agree.
  Honest note: it's a threshold/peak effect, so `ld_w` is a **prioritization**,
  not proof.

## 5. Region-centric resolution: overlapping LD components
- `ld_w`-filtered clusters on the 500 kb Stage-2 clustering: the inversion
  resolves into **two physically-interspersed LD components** (the two
  arrangements) — LD topology, not physical distance. A window method would lump
  them.
- Beat: LDscnR resolves overlapping components a physical-window scan cannot.

## 6. Summary
- Prune → LOCO GRM → single-SNP scan prioritized by a data-driven `ld_w`
  threshold: catches localized (*Eda*) and block (inversion) signals the
  cluster-consensus (Li et al.) approach misses; `ld_w` screens isolated false
  positives; region-centric clustering resolves overlapping components.

## Figures (all from dev/outlier/analysis.R)
D, A, A2, E, B, + the 500 kb `ld_w`-cluster figure (§5). Stats: strategy table,
bimodality table, per-chr.

## Open questions for review
1. Keep §5 (overlapping components) here, or leave it in the complexity-reduction
   vignette's scope? (It's arguably a clustering point, but it lands hardest in
   the outlier context.)
2. A/A2 both, or just A (Zissou) — A2 (alpha-only) is the cleaner "what remains".
3. Vignette recomputes the pipeline (~40s build). Acceptable, or precompute?
