## dev/outlier/build_cache.R
## Build the expensive pipeline + EMMAX objects ONCE and cache them, so
## dev/outlier/analysis.R can iterate on figures/stats without recomputing.
## Run from the LDscnR repo root:  Rscript dev/outlier/build_cache.R
## (cache.rds and fig/ are gitignored; dev/ is .Rbuildignore'd.)

suppressMessages({
  devtools::load_all(".")
  library(data.table); library(SNPRelate)
})

data(stickleback); GTs <- stickleback$GTs; map <- copy(stickleback$map); pheno <- stickleback$pheno
gp <- tempfile(fileext = ".gds"); gds <- create_gds_from_geno(GTs, map, gp)
set.seed(1)

## ---- pipeline -------------------------------------------------------------
ld <- compute_LD_decay(gds, keep_el = TRUE, slide = 200, ld_method = "corr",
                       n_win_decay = 50, cores = 4)
lw <- compute_ld_w(ld, rho = 0.95, cores = 4)
ids <- unlist(lapply(ld$by_chr, function(o) o$snp_ids), use.names = FALSE)
map[, ld_w_095 := lw[match(marker, ids)]]
s1 <- ld_complexity_reduction(map = map, LD_decay = ld, rho = 0.5, cores = 4)

COMM <- list(GTs = GTs, stage1 = s1, ld_w_col = "ld_w_095", ld_w_threshold = 0.05,
             LD_decay = ld, rho = 0.95, score_threshold = 0.80, min_r2 = 0.2,
             min_n_loci_flag = 5, min_n_loci_eMLG = 5, compute_unflagged_eMLG = TRUE, cores = 4)
result     <- do.call(ld_prune_and_eMLG, COMM)
result_500 <- do.call(ld_prune_and_eMLG, c(COMM, list(distance_threshold = 5e5)))
best <- eMLG_best_snp(result, GTs, round_fill = FALSE)

eco  <- as.integer(pheno$ecotype == "Marine")
chrs <- unique(map$Chr)
pruned_chr <- map$Chr[match(result$pruned, map$marker)]

g <- as.data.table(result$groups)[has_eMLG == TRUE]
g[, best_marker := best$stats$best_marker[match(group_id, best$stats$group_id)]]

## ---- fixed-GRM representation comparison (rep / consensus / max) ----------
K <- snpgdsGRM(gds, snp.id = result$pruned, method = "GCTA", autosome.only = FALSE, verbose = FALSE)$grm
g[, `:=`(Fr_fix = emmax(eco, GTs[, representative], K)$F,
         Fc_fix = emmax(eco, result$eMLG[, group_id], K)$F,
         Fm_fix = emmax(eco, GTs[, best_marker],     K)$F)]

## ---- LOCO GRMs + scans ----------------------------------------------------
Kloco <- setNames(lapply(chrs, function(ch)
  snpgdsGRM(gds, snp.id = result$pruned[pruned_chr != ch], method = "GCTA",
            autosome.only = FALSE, verbose = FALSE)$grm), chrs)

# single-SNP LOCO (F + p per marker)
p_snp <- setNames(rep(NA_real_, nrow(map)), map$marker); F_snp <- p_snp
for (ch in chrs) { mk <- map[Chr == ch, marker]
  e <- emmax(eco, GTs[, mk], Kloco[[ch]]); p_snp[mk] <- e$pval; F_snp[mk] <- e$F }
map[, `:=`(p_loco = p_snp[marker], F_loco = F_snp[marker])]

# consensus LOCO (F + p per eMLG cluster)
g[, `:=`(p_clust_loco = NA_real_, Fc_loco = NA_real_)]
for (ch in chrs) { gi <- which(g$Chr == ch)
  if (length(gi)) { e <- emmax(eco, result$eMLG[, g$group_id[gi]], Kloco[[ch]])
    g$p_clust_loco[gi] <- e$pval; g$Fc_loco[gi] <- e$F } }

snpgdsClose(gds); unlink(gp)

saveRDS(list(
  map               = map,                              # ld_w_095, p_loco, F_loco
  g                 = g,                                # cluster F/p (fixed + LOCO), n_loci, ...
  result_groups     = as.data.table(result$groups),    # members (default distance)
  result_500_groups = as.data.table(result_500$groups),# members (500 kb)
  pruned            = result$pruned,
  eco               = eco
), "dev/outlier/cache.rds")

cat("cache written: dev/outlier/cache.rds\n")
cat(sprintf("  SNPs %d | eMLG clusters %d | pruned %d\n", nrow(map), nrow(g), length(result$pruned)))
