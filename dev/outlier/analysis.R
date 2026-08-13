## dev/outlier/analysis.R
## Loads dev/outlier/cache.rds (built by build_cache.R) and produces the
## LD-aware outlier-analysis figures/stats WITHOUT recomputing the pipeline.
## Iterate here; when settled we design the outlier vignette from this.
## Run from the LDscnR repo root:  Rscript dev/outlier/analysis.R

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
o <- readRDS("dev/outlier/cache.rds")
map <- copy(o$map); g <- copy(o$g)
rg  <- o$result_groups; rg500 <- o$result_500_groups
dir.create("dev/outlier/fig", showWarnings = FALSE, recursive = TRUE)
chrs <- unique(map$Chr)
ZISSOU <- c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")   # wesanderson Zissou1
sv <- function(name, p, w = 200, h = 80) ggsave(file.path("dev/outlier/fig", name), p, width = w, height = h, units = "mm", dpi = 130)

## ---- FDR + categories -----------------------------------------------------
map[, q_all := p.adjust(p_loco, "fdr")]
map[, top_ldw := ld_w_095 >= quantile(ld_w_095, 0.90, na.rm = TRUE)]
map[top_ldw == TRUE, q_top := p.adjust(p_loco[top_ldw], "fdr")]
g[, q_clust := p.adjust(p_clust_loco, "fdr")]
map[, neglogq := -log10(q_all)]

cat("=== strategy comparison (LOCO, q<0.05) ===\n")
print(data.frame(strategy = c("all SNPs","top 10% ld_w","consensus clusters"),
                 n_tests = c(nrow(map), sum(map$top_ldw), nrow(g)),
                 n_sig   = c(sum(map$q_all<0.05,na.rm=T), sum(map$q_top<0.05,na.rm=T), sum(g$q_clust<0.05))))

## =========================================================================
## (A) NEW: Manhattan coloured by ld_w (Zissou)
## =========================================================================
pA <- ggplot(map, aes(Pos/1e6, neglogq, colour = ld_w_095)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey60") +
  geom_point(size = 0.8) +
  scale_colour_gradientn(colours = ZISSOU, name = expression(ld[w])) +
  facet_wrap(~ Chr, scales = "free_x", nrow = 1) +
  labs(x = "position (Mbp)", y = expression(-log[10](q)~", LOCO single-SNP"),
       title = "Single-SNP scan coloured by ld_w (Zissou)") +
  theme_bw(base_size = 9)
sv("A_manhattan_ldw_zissou.png", pA, h = 80)

## =========================================================================
## (B) NEW: F vs ld_w quantile, per chromosome (visual) + correlations
## =========================================================================
map[, ldw_q := ecdf(ld_w_095)(ld_w_095), by = Chr]
map[, is_rep := marker %in% o$pruned]
pB <- ggplot(map, aes(ldw_q, F_loco)) +
  geom_point(size = 0.4, alpha = 0.3, colour = "grey40") +
  geom_smooth(method = "loess", se = TRUE, colour = "#F21A00", linewidth = 0.6) +
  facet_wrap(~ Chr, nrow = 1) +
  labs(x = "ld_w quantile (within chromosome)", y = "EMMAX F (LOCO single-SNP)",
       title = "F vs ld_w quantile per chromosome") +
  theme_bw(base_size = 9)
sv("B_F_vs_ldwquantile.png", pB, h = 75)

cat("\n=== per-chr Spearman(F, ld_w): all-SNP | representatives ===\n")
for (ch in chrs) { d <- map[Chr==ch]; dr <- d[is_rep==TRUE]
  cat(sprintf("%-5s all=%+.3f (n=%d) | rep=%+.3f p=%.2g (n=%d)\n", ch,
      cor(d$F_loco, d$ld_w_095, method="spearman"), nrow(d),
      suppressWarnings(cor.test(dr$F_loco, dr$ld_w_095, method="spearman"))$estimate,
      suppressWarnings(cor.test(dr$F_loco, dr$ld_w_095, method="spearman"))$p.value, nrow(dr))) }

## =========================================================================
## (C) consensus F vs cluster size, per chromosome
## =========================================================================
pC <- ggplot(g, aes(n_loci, Fc_loco)) +
  geom_point(size = 0.8, alpha = 0.5, colour = "#3B9AB2") +
  geom_smooth(method = "loess", se = TRUE, colour = "#F21A00", linewidth = 0.6) +
  scale_x_log10() + facet_wrap(~ Chr, nrow = 1) +
  labs(x = "cluster size (n markers, log)", y = "consensus EMMAX F (LOCO)",
       title = "Consensus F vs cluster size per chromosome") +
  theme_bw(base_size = 9)
sv("C_consensusF_vs_size.png", pC, h = 75)

cat("\n=== per-chr Spearman(consensus F, cluster size) ===\n")
for (ch in chrs) { d <- g[Chr==ch]; ct <- suppressWarnings(cor.test(d$Fc_loco, d$n_loci, method="spearman"))
  cat(sprintf("%-5s rho=%+.3f p=%.2g (n=%d)\n", ch, ct$estimate, ct$p.value, nrow(d))) }

## =========================================================================
## (D) reference: LOCO Manhattan with significant consensus clusters + black triangles
## =========================================================================
pal <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF","#1B9E77")
sig <- rg[group_id %in% g[q_clust<0.05, group_id]]
if (nrow(sig)) {
  sig[, pos0 := vapply(members, function(mk) min(map$Pos[match(mk, map$marker)]), numeric(1))]
  setorder(sig, Chr, pos0); sig[, col := pal[(seq_len(.N)-1)%%length(pal)+1]]
  s2c <- rbindlist(lapply(seq_len(nrow(sig)), function(i) data.table(marker=sig$members[[i]], col=sig$col[i])))
  map[, sig_col := s2c$col[match(marker, s2c$marker)]]
} else map[, sig_col := NA_character_]
map[, cat := fifelse(!is.na(sig_col), "cluster", fifelse(q_all<0.05, "snp_only", "ns"))]
pD <- ggplot() +
  geom_hline(yintercept=-log10(0.05), linetype=2, colour="grey60") +
  geom_point(data=map[cat=="ns"], aes(Pos/1e6, neglogq), colour="grey78", size=0.5, alpha=0.5) +
  geom_point(data=map[cat=="snp_only"], aes(Pos/1e6, neglogq), colour="black", shape=17, size=1.4) +
  geom_point(data=map[cat=="cluster"], aes(Pos/1e6, neglogq, colour=sig_col), size=1.6) +
  scale_colour_identity() + facet_wrap(~Chr, scales="free_x", nrow=1) +
  labs(x="position (Mbp)", y=expression(-log[10](q)), title="LOCO single-SNP; sig consensus clusters coloured, black = SNP-only") +
  theme_bw(base_size=9)
sv("D_loco_manhattan.png", pD, h = 80)

cat("\nWrote figures to dev/outlier/fig/: A (ldw Zissou), B (F~ld_w quantile), C (consensus F~size), D (LOCO manhattan)\n")
