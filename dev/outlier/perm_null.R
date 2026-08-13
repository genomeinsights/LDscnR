## dev/outlier/perm_null.R
## Permutation null for the cluster-level filter of the candidate-first pipeline.
## Clusters are genotype-only (fixed). We permute the phenotype B times, re-run
## the LOCO single-SNP scan on the candidates, and count "significant" SNPs
## (p <= observed BH cutoff) per FIXED cluster. This gives a size-aware null for
## n_sig per cluster (LD-clumping captured, unlike an analytic binomial), so each
## observed cluster gets an empirical p-value.
## Run from repo root:  Rscript dev/outlier/perm_null.R

suppressMessages({ library(LDscnR); library(data.table); library(ggplot2); library(SNPRelate) })
o   <- readRDS("dev/outlier/cache.rds"); map <- copy(o$map)
data("stickleback", package = "LDscnR"); GTs <- stickleback$GTs; pheno <- stickleback$pheno
stopifnot(all(map$marker %in% colnames(GTs)))
eco  <- as.integer(pheno$ecotype == "Marine")
chrs <- unique(map$Chr)
B    <- 1000
set.seed(1)

## ---- RMSC -> q*, candidates, observed FDR significance ---------------------
thr_qs <- seq(0, 0.98, by = 0.01)
rej <- vapply(thr_qs, function(q) { k <- map$ld_w_095 >= quantile(map$ld_w_095, q, na.rm = TRUE)
  sum(p.adjust(map$p_loco[k], "fdr") < 0.05, na.rm = TRUE) }, integer(1))
best_thr <- thr_qs[which.max(rej)]
map[, cand := ld_w_095 >= quantile(ld_w_095, best_thr, na.rm = TRUE)]
map[, q_top := NA_real_]; map[cand == TRUE, q_top := p.adjust(p_loco, "fdr")]
pcut <- map[cand == TRUE & q_top < 0.05, max(p_loco)]   # fixed per-SNP significance cutoff
cat(sprintf("q* = %.2f | candidates = %d | BH p-cutoff = %.4g\n", best_thr, sum(map$cand), pcut))

## ---- fixed clusters: single-linkage r^2 >= 0.5 on candidates ---------------
RT <- 0.5
map[, clust := NA_character_]
for (ch in chrs) {
  idx <- which(map$cand & map$Chr == ch)
  if (!length(idx)) next
  mk <- map$marker[idx]
  if (length(mk) == 1) { map$clust[idx] <- paste0(ch, "_c1"); next }
  R  <- cor(GTs[, mk])^2
  map$clust[idx] <- paste0(ch, "_c", cutree(hclust(as.dist(1 - R), method = "single"), h = 1 - RT))
}
cand_dt <- map[cand == TRUE, .(marker, Chr, Pos, clust, p_loco)]
cl <- cand_dt[, .(n = .N, lo = min(Pos)/1e6, hi = max(Pos)/1e6,
                  n_sig_obs = sum(p_loco <= pcut)), by = clust]

## ---- LOCO GRMs (rebuilt from the cached pruned set) ------------------------
gp <- tempfile(fileext = ".gds"); gds <- create_gds_from_geno(GTs, map, gp)
pruned_chr <- map$Chr[match(o$pruned, map$marker)]
Kloco <- setNames(lapply(chrs, function(ch)
  snpgdsGRM(gds, snp.id = o$pruned[pruned_chr != ch], method = "GCTA",
            autosome.only = FALSE, verbose = FALSE)$grm), chrs)
snpgdsClose(gds); unlink(gp)

## ---- permutation loop ------------------------------------------------------
cand_mk_by_chr <- setNames(lapply(chrs, function(ch) cand_dt[Chr == ch, marker]), chrs)
null_mat <- matrix(0L, nrow = B, ncol = nrow(cl), dimnames = list(NULL, cl$clust))
for (b in seq_len(B)) {
  yb <- sample(eco)
  pv <- setNames(rep(NA_real_, nrow(cand_dt)), cand_dt$marker)
  for (ch in chrs) {
    mk <- cand_mk_by_chr[[ch]]; if (!length(mk)) next
    pv[mk] <- emmax(yb, GTs[, mk, drop = FALSE], Kloco[[ch]])$pval
  }
  ns <- cand_dt[, .(k = sum(pv[.I] <= pcut)), by = clust]   # n_sig per cluster this perm
  null_mat[b, ns$clust] <- ns$k
}

## ---- empirical p-values + size-aware null band -----------------------------
cl[, `:=`(null_mean = colMeans(null_mat)[clust],
          null_q95  = apply(null_mat, 2, quantile, 0.95)[clust])]
cl[, emp_p := vapply(seq_len(.N), function(i)
      (1 + sum(null_mat[, clust[i]] >= n_sig_obs[i])) / (B + 1), numeric(1))]
cl[, emp_q := p.adjust(emp_p, "fdr")]
setorder(cl, emp_p, -n_sig_obs)
cat(sprintf("\n=== cluster-level permutation null (B=%d), r2>=%.1f ===\n", B, RT))
print(cl[, .(clust, n, n_sig_obs, null_mean = round(null_mean, 2),
             null_q95, emp_p = round(emp_p, 4), emp_q = round(emp_q, 4),
             region = sprintf("%.2f-%.2f", lo, hi))])

## ---- figure: observed n_sig vs cluster size, with null 95th-pct band -------
cl[, signif := emp_p < 0.05]
pf <- cl[n_sig_obs > 0]   # the flagged clusters
pf[, lab := sprintf("%s (%.1f Mb)", sub("_c.*", "", clust), (lo + hi)/2)]
p <- ggplot(pf, aes(n, n_sig_obs)) +
  geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey70") +
  geom_linerange(aes(ymin = n_sig_obs, ymax = null_q95), colour = "grey80") +
  geom_point(aes(y = null_q95), colour = "grey45", shape = 4, size = 2.4) +
  geom_point(aes(colour = signif), size = 3) +
  ggrepel::geom_text_repel(aes(label = lab), size = 2.6, seed = 1) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c(`TRUE` = "#F21A00", `FALSE` = "grey55"),
                      name = "cluster sig\n(emp p<0.05)") +
  labs(x = "cluster size (candidates, log)", y = "observed n_sig (log)",
       title = sprintf("Cluster-level permutation null: obs vs null 95%%-pct (x), B=%d", B),
       subtitle = "grey x = null 95th-pctile for that cluster; dotted = n_sig = size") +
  theme_bw(base_size = 9) + theme(strip.background = element_blank())
ggsave("dev/outlier/fig/I_perm_null.png", p, width = 155, height = 100, units = "mm", dpi = 130)
cat("\nwrote dev/outlier/fig/I_perm_null.png\n")
