## dev/outlier/cluster_on_candidates.R
## Prototype of the INVERTED (candidate-first) LD-aware outlier pipeline:
##   1. RMSC on (ld_w, single-SNP F/p) -> data-driven ld_w threshold q*
##   2. candidates = SNPs with ld_w >= q* quantile   (only these can be outliers)
##   3. cluster ONLY the candidates by LD            (cheap: small set)
##   4. single-SNP FDR within candidates -> "outliers" (significant SNPs)
##   5. a cluster is FLAGGED if it contains >=1 significant SNP; report the
##      whole cluster (significant SNPs + SNPs in LD with them) as the unit.
## In sim data each cluster is the test unit, scored TP/FP by whether ANY member
## is causal. Here (stickleback) we have no causal labels -- structure only.
## Run from repo root:  Rscript dev/outlier/cluster_on_candidates.R

suppressMessages({ library(data.table); library(ggplot2) })
o   <- readRDS("dev/outlier/cache.rds"); map <- copy(o$map)
data("stickleback", package = "LDscnR"); GTs <- stickleback$GTs
stopifnot(all(map$marker %in% colnames(GTs)))
ZISSOU <- c("#3B9AB2","#78B7C5","#EBCC2A","#E1AF00","#F21A00")

## ---- 1. RMSC -> q* --------------------------------------------------------
thr_qs <- seq(0, 0.98, by = 0.01)
rej <- vapply(thr_qs, function(q) {
  k <- map$ld_w_095 >= quantile(map$ld_w_095, q, na.rm = TRUE)
  sum(p.adjust(map$p_loco[k], "fdr") < 0.05, na.rm = TRUE)
}, integer(1))
best_thr <- thr_qs[which.max(rej)]

## ---- 2. candidates + 4. outliers (FDR within candidates) ------------------
map[, cand := ld_w_095 >= quantile(ld_w_095, best_thr, na.rm = TRUE)]
map[, q_top := NA_real_]; map[cand == TRUE, q_top := p.adjust(p_loco, "fdr")]
map[, sig := !is.na(q_top) & q_top < 0.05]
cat(sprintf("q* = %.2f | candidates = %d | outlier SNPs (q<0.05) = %d\n",
            best_thr, sum(map$cand), sum(map$sig)))

## ---- 3. cluster ONLY the candidates, per chromosome -----------------------
## single-linkage connected components at r^2 >= r2_thresh: a SNP joins a cluster
## if it is in LD (r^2 >= thresh) with ANY member -> "SNPs in LD with each other".
cluster_candidates <- function(r2_thresh) {
  cl_col <- rep(NA_character_, nrow(map))
  for (ch in unique(map$Chr)) {
    idx <- which(map$cand & map$Chr == ch)
    if (!length(idx)) next
    mk <- map$marker[idx]
    if (length(mk) == 1) { cl_col[idx] <- paste0(ch, "_c1"); next }
    R  <- cor(GTs[, mk])^2
    cl <- cutree(hclust(as.dist(1 - R), method = "single"), h = 1 - r2_thresh)
    cl_col[idx] <- paste0(ch, "_c", cl)
  }
  cl_col
}

## ---- sensitivity across the LD-linkage threshold --------------------------
sens <- rbindlist(lapply(c(0.2, 0.5, 0.8), function(rt) {
  map[, clust := cluster_candidates(rt)]
  cl <- map[cand == TRUE, .(n = .N, n_sig = sum(sig),
                            lo = min(Pos)/1e6, hi = max(Pos)/1e6), by = .(Chr, clust)]
  flagged <- cl[n_sig > 0]
  data.table(r2_thresh = rt,
             n_clusters      = nrow(cl),
             n_flagged       = nrow(flagged),
             snps_in_flagged = flagged[, sum(n)],
             recruited       = flagged[, sum(n)] - sum(map$sig),  # non-sig members pulled in
             max_flagged_sz  = flagged[, max(n)])
}))
cat("\n=== LD-linkage sensitivity (cluster = test unit) ===\n"); print(sens)

## ---- fix r2_thresh = 0.5 for the reported result --------------------------
RT <- 0.5
map[, clust := cluster_candidates(RT)]
cl  <- map[cand == TRUE, .(n = .N, n_sig = sum(sig),
                           lo = min(Pos)/1e6, hi = max(Pos)/1e6), by = .(Chr, clust)]
flagged <- cl[n_sig > 0][order(Chr, lo)]
cat(sprintf("\n=== flagged clusters (r2>=%.1f): %d clusters, %d SNPs (%d significant + %d recruited) ===\n",
            RT, nrow(flagged), flagged[, sum(n)], sum(map$sig), flagged[, sum(n)] - sum(map$sig)))
print(flagged)

## ---- Manhattan: colour SNPs by FLAGGED cluster (sig + LD-linked) -----------
pal <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF","#1B9E77",
         "#66C2A5","#FC8D62","#8DA0CB","#E78AC3")
flagged[, col := pal[(seq_len(.N) - 1) %% length(pal) + 1]]
map[, flag_col := flagged$col[match(clust, flagged$clust)]]
map[, in_flagged := !is.na(flag_col)]

p <- ggplot() +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey60") +
  geom_point(data = map[in_flagged == FALSE], aes(Pos/1e6, -log10(p_loco)),
             colour = "grey82", size = 0.5, alpha = 0.5) +
  geom_point(data = map[in_flagged == TRUE], aes(Pos/1e6, -log10(p_loco), colour = flag_col),
             size = 1.5) +
  scale_colour_identity() +
  facet_wrap(~ Chr, scales = "free_x", nrow = 1) +
  labs(x = "position (Mbp)", y = expression(-log[10](p)~", LOCO single-SNP"),
       title = sprintf("Candidate-first clustering: flagged clusters (sig SNP + LD-linked), r2>=%.1f", RT)) +
  theme_bw(base_size = 9) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
ggsave("dev/outlier/fig/H_candidate_clusters.png", p, width = 200, height = 80, units = "mm", dpi = 130)
cat("\nwrote dev/outlier/fig/H_candidate_clusters.png\n")
