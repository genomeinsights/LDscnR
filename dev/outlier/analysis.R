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
g[, q_clust := p.adjust(p_clust_loco, "fdr")]
map[, neglogq := -log10(q_all)]
map[, ldw_rank := ecdf(ld_w_095)(ld_w_095)]   # genome-wide ld_w rank (0-1)

## data-driven ld_w threshold: rejection-maximizing (independent filtering; Bourgon 2010)
thr_qs  <- seq(0, 0.98, by = 0.01)
thr_rej <- vapply(thr_qs, function(q) {
  keep <- map$ld_w_095 >= quantile(map$ld_w_095, q, na.rm = TRUE)
  sum(p.adjust(map$p_loco[keep], "fdr") < 0.05, na.rm = TRUE)
}, integer(1))
best_thr <- thr_qs[which.max(thr_rej)]
map[, top_ldw := ld_w_095 >= quantile(ld_w_095, best_thr, na.rm = TRUE)]
map[top_ldw == TRUE, q_top := p.adjust(p_loco[top_ldw], "fdr")]

cat(sprintf("=== data-driven ld_w threshold: q* = %.2f (rejections = %d) ===\n", best_thr, max(thr_rej)))
cat("=== strategy comparison (LOCO, q<0.05) ===\n")
print(data.frame(strategy = c("all SNPs", sprintf("ld_w q>=%.2f", best_thr), "consensus clusters"),
                 n_tests = c(nrow(map), sum(map$top_ldw), nrow(g)),
                 n_sig   = c(sum(map$q_all<0.05,na.rm=T), sum(map$q_top<0.05,na.rm=T), sum(g$q_clust<0.05))))

## (E) threshold-selection curve
pE <- ggplot(data.table(q = thr_qs, rej = thr_rej), aes(q, rej)) +
  geom_line(colour = "grey45") + geom_point(size = 0.7) +
  geom_vline(xintercept = best_thr, linetype = 2, colour = "#F21A00") +
  annotate("text", x = best_thr, y = max(thr_rej), hjust = -0.1,
           label = sprintf("q* = %.2f", best_thr), colour = "#F21A00", size = 3) +
  labs(x = "ld_w quantile threshold (genome-wide)", y = "FDR-significant SNPs",
       title = "Data-driven ld_w threshold: rejection maximization (independent filtering)") +
  theme_bw(base_size = 9) + theme(strip.background = element_blank())
sv("E_threshold_selection.png", pE, w = 130, h = 80)

## =========================================================================
## (A) Manhattan: colour AND alpha = ld_w (Zissou) so high-ld_w peaks pop
## =========================================================================
pA <- ggplot(map, aes(Pos/1e6, neglogq, colour = ld_w_095, alpha = ld_w_095)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey60") +
  geom_point(size = 0.9) +
  scale_colour_gradientn(colours = ZISSOU, name = expression(ld[w])) +
  scale_alpha(range = c(0, 1), guide = "none") +
  facet_wrap(~ Chr, scales = "free_x", nrow = 1) +
  labs(x = "position (Mbp)", y = expression(-log[10](q)~", LOCO single-SNP"),
       title = "Single-SNP scan: colour + alpha = ld_w (Zissou)") +
  theme_bw(base_size = 9) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
sv("A_manhattan_ldw_zissou.png", pA, h = 80)

## (A2) alpha-only on RAW ld_w (single colour) -- what peaks remain when low-ld_w fades
pA2 <- ggplot(map, aes(Pos/1e6, neglogq, alpha = ld_w_095)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey60") +
  geom_point(size = 0.9, colour = "black") +
  scale_alpha(range = c(0, 1), name = expression(ld[w])) +
  facet_wrap(~ Chr, scales = "free_x", nrow = 1) +
  labs(x = "position (Mbp)", y = expression(-log[10](q)~", LOCO single-SNP"),
       title = "Single-SNP scan: alpha = ld_w (raw) only") +
  theme_bw(base_size = 9) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
sv("A2_manhattan_alpha_only.png", pA2, h = 80)

## =========================================================================
## (B) F vs GENOME-WIDE ld_w quantile; rolling median + data-driven threshold line
## =========================================================================
setorder(map, Chr, ldw_rank)
map[, F_roll := frollapply(F_loco, 201, median, align = "center"), by = Chr]
pB <- ggplot(map, aes(ldw_rank, F_loco)) +
  geom_point(size = 0.4, alpha = 0.25, colour = "grey55") +
  geom_line(aes(y = F_roll), colour = "#F21A00", linewidth = 0.7, na.rm = TRUE) +
  geom_vline(xintercept = best_thr, linetype = 2, colour = "#3B9AB2") +
  facet_wrap(~ Chr, nrow = 1) +
  labs(x = "ld_w quantile (genome-wide)", y = "EMMAX F (LOCO single-SNP)",
       title = sprintf("F vs genome-wide ld_w quantile; rolling median; data-driven q* = %.2f", best_thr)) +
  theme_bw(base_size = 9) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
sv("B_F_vs_ldwquantile.png", pB, h = 75)

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
  theme_bw(base_size = 9) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
sv("D_loco_manhattan.png", pD, h = 80)

cat("\nWrote: A, A2, B (F~ld_w + data-driven thr), D (LOCO manhattan), E (threshold selection)\n")
