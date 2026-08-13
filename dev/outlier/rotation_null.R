## dev/outlier/rotation_null.R
## Method-agnostic cluster-level null from the BELOW-q* SNPs as genomic
## background, compared to the EMMAX permutation null. The association method is
## run ONCE (cached p_loco); no refit -> works with ANY per-SNP statistic.
##
## Schemes compared per flagged cluster:
##   perm  : EMMAX phenotype permutation (ref, from perm_null.R, B=1000)
##   rotB  : rotate significance WITHIN candidates, per chromosome  (fails: no
##           null background -- candidates are pre-enriched for signal)
##   bg    : draw null n_sig from the BELOW-q* background field, size-matched,
##           rotated to preserve background autocorrelation  (user's idea)
## Run from repo root:  Rscript dev/outlier/rotation_null.R

suppressMessages({ library(data.table); library(ggplot2) })
o   <- readRDS("dev/outlier/cache.rds"); map <- copy(o$map)
data("stickleback", package = "LDscnR"); GTs <- stickleback$GTs
stopifnot(all(map$marker %in% colnames(GTs)))
B <- 5000; set.seed(1)

## permutation emp_p from perm_null.R (B=1000), for reference
perm_emp <- c(Chr4_c1 = 0.0010, Chr1_c3 = 0.0080, Chr4_c3 = 0.0120,
              Chr6_c1 = 0.1219, Chr4_c4 = 0.2478)
rotate <- function(v, k) { n <- length(v); v[((seq_len(n) - 1 + k) %% n) + 1] }

## ---- candidates, clusters (RT=0.5), observed significance ------------------
thr_qs <- seq(0, 0.98, by = 0.01)
rej <- vapply(thr_qs, function(q) { k <- map$ld_w_095 >= quantile(map$ld_w_095, q, na.rm = TRUE)
  sum(p.adjust(map$p_loco[k], "fdr") < 0.05, na.rm = TRUE) }, integer(1))
best_thr <- thr_qs[which.max(rej)]
map[, cand := ld_w_095 >= quantile(ld_w_095, best_thr, na.rm = TRUE)]
map[, q_top := NA_real_]; map[cand == TRUE, q_top := p.adjust(p_loco, "fdr")]
pcut <- map[cand == TRUE & q_top < 0.05, max(p_loco)]
RT <- 0.5; map[, clust := NA_character_]
for (ch in unique(map$Chr)) {
  idx <- which(map$cand & map$Chr == ch); if (!length(idx)) next
  mk <- map$marker[idx]
  if (length(mk) == 1) { map$clust[idx] <- paste0(ch, "_c1"); next }
  R <- cor(GTs[, mk])^2
  map$clust[idx] <- paste0(ch, "_c", cutree(hclust(as.dist(1 - R), method = "single"), h = 1 - RT))
}
map[, sig := as.integer(p_loco <= pcut)]
obs <- map[cand == TRUE, .(n = .N, n_sig_obs = sum(sig)), by = clust]

## background significance rate (below-q* SNPs)
bg_rate <- map[cand == FALSE, mean(sig)]
cat(sprintf("q*=%.2f | pcut=%.4g | candidates=%d | below-q* background rate=%.4f\n",
            best_thr, pcut, sum(map$cand), bg_rate))

## ---- rotB: rotate significance within candidates, per chromosome -----------
co <- map[cand == TRUE][order(Chr, Pos)]
nullB <- setNames(vector("list", nrow(obs)), obs$clust)
idx_by_chr <- split(seq_len(nrow(co)), co$Chr)
for (b in seq_len(B)) {
  sr <- co$sig
  for (ix in idx_by_chr) if (length(ix) > 1) sr[ix] <- rotate(co$sig[ix], sample.int(length(ix) - 1, 1))
  t <- tapply(sr, co$clust, sum)
  for (cl in names(t)) nullB[[cl]][b] <- t[[cl]]
}

## ---- bg: null n_sig drawn from BELOW-q* background, size-matched -----------
nullBG <- setNames(vector("list", nrow(obs)), obs$clust)
for (ch in unique(map$Chr)) {
  sbg <- map[Chr == ch & cand == FALSE][order(Pos), sig]; nb <- length(sbg)
  for (cl in obs[startsWith(clust, paste0(ch, "_")), clust]) {
    ni <- obs[clust == cl, n]
    nullBG[[cl]] <- vapply(seq_len(B), function(b) {
      st <- sample.int(nb, 1); sum(sbg[((st - 1 + seq_len(ni) - 1) %% nb) + 1])
    }, numeric(1))
  }
}

## ---- empirical p-values ----------------------------------------------------
emp <- function(nulls, cl) (1 + sum(nulls[[cl]] >= obs[clust == cl, n_sig_obs])) / (B + 1)
res <- obs[n_sig_obs > 0][order(-n_sig_obs)]
res[, `:=`(perm_p = perm_emp[clust],
           rotB_p = vapply(clust, function(c) emp(nullB,  c), numeric(1)),
           bg_p   = vapply(clust, function(c) emp(nullBG, c), numeric(1)),
           bg_q95 = vapply(clust, function(c) quantile(nullBG[[c]], 0.95), numeric(1)))]
cat("\n=== permutation vs background (below-q*) vs candidate-rotation ===\n")
print(res[, .(clust, n, n_sig_obs, bg_q95,
              perm_p = round(perm_p, 4), bg_p = round(bg_p, 4), rotB_p = round(rotB_p, 4))])

## ---- figure ---------------------------------------------------------------
m <- melt(res[, .(clust, perm_p, bg_p, rotB_p)], id.vars = "clust",
          variable.name = "method", value.name = "emp_p")
m[, method := factor(method, c("perm_p","bg_p","rotB_p"),
                     c("permutation", "below-q* background", "candidate rotation"))]
m[, reg := sub("_c.*", "", clust)]
p <- ggplot(m, aes(method, pmax(emp_p, 1/(B+1)), group = clust, colour = reg)) +
  geom_hline(yintercept = 0.05, linetype = 2, colour = "grey60") +
  geom_line(alpha = 0.5) + geom_point(size = 2.6) +
  ggrepel::geom_text_repel(data = m[method == "below-q* background"],
                           aes(label = clust), size = 2.4, seed = 1) +
  scale_y_log10() +
  labs(x = NULL, y = "empirical p (log)", colour = "region",
       title = sprintf("Cluster-level null: permutation vs below-q* background (B=%d)", B),
       subtitle = "dashed = 0.05; lines join the same cluster across methods") +
  theme_bw(base_size = 9) + theme(strip.background = element_blank())
ggsave("dev/outlier/fig/J_perm_vs_background.png", p, width = 165, height = 100, units = "mm", dpi = 130)
cat("\nwrote dev/outlier/fig/J_perm_vs_background.png\n")
