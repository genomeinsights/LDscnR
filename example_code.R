library(LDscnR)
library(data.table)

# ------------------------------------------------------------
# 0. (One-off) how the bundled example data was created
# ------------------------------------------------------------
# tmp <- readRDS("../LD-scaling-genome-scans/results_sim/c1_V0.5_rep4.rds")
# map <- tmp$SNP_res[,.(V,c,rep,Chr,Pos,marker,type, Chr_type,max_LD_with_QTN, bp_to_focal_QTN, focal_QTN, p_Va,emx_F,lfmm_F,emx_q,lfmm_q )]
# GTs <- tmp$GTs
#
# sim_ex <- list(map = map, GTs = GTs)
#
# # Save in package format
# usethis::use_data(sim_ex, overwrite = TRUE)
# usethis::use_r("data_sim_ex")

# ------------------------------------------------------------
# 1. Create the GDS object used throughout
# ------------------------------------------------------------

data("sim_ex")

map <- sim_ex$map
GTs <- sim_ex$GTs
colnames(GTs) <- map$marker

gds_path <- tempfile(fileext = ".gds")
gds <- create_gds_from_geno(geno = sim_ex$GTs, map, gds_path)

# ------------------------------------------------------------
# 2. Estimate chromosome-wise LD decay
# ------------------------------------------------------------

# number of cores used throughout
cores <- 8

## Preliminary LD-decay on a subset of SNPs, to choose a good sliding window
## for the full run. n_win_decay is 5 (not the default 20) because the
## simulated chromosomes are only ~10 cM; slide is large because we analyse a
## subset of SNPs only.
ld_decay_pre <- compute_LD_decay(
  gds,
  n_win_decay    = 5,
  max_SNPs_decay = 5000,
  slide          = 1000,
  cores          = cores
)
ld_decay_pre

## Recommended new slide window for the full run (mean across chromosomes at rho = 0.99)
new_slide_window <- mean(ld_decay_pre$decay_sum$slide_snp_rho_0.99)

## Full LD-decay using all SNPs, keeping the edge lists
ld_decay <- compute_LD_decay(
  gds,
  n_win_decay    = 5,
  max_SNPs_decay = Inf,   # all SNPs used now
  keep_el        = TRUE,  # keep edge lists
  slide          = new_slide_window,
  cores          = cores,
  ld_method      = "r"
)

## Alternative: write edge lists to disk instead of RAM (large data sets)
# ld_decay <- compute_LD_decay(
#   gds, el_data_folder = "./el/",
#   n_win_decay = 5, max_SNPs_decay = Inf, keep_el = TRUE,
#   slide = new_slide_window, cores = 1
# )

# ------------------------------------------------------------
# 3. Inspect and visualise the decay fit
# ------------------------------------------------------------

## Chromosome-wise decay parameters (a = decay rate, b = background LD, c = short-range LD)
ld_decay$decay_sum

## Recommended sliding-window sizes (in SNP units) per rho threshold
ld_decay$recommendation

## Diagnostic plots
plot(ld_decay, type = "summary")
plot(ld_decay, type = "recommendation")
plot(ld_decay, type = "recommendation", rho = 0.99)
plot(ld_decay, type = "chr", chr = "Chr2")

# ------------------------------------------------------------
# 4. Convert between rho, physical distance, and expected r^2
# ------------------------------------------------------------

chr1 <- ld_decay$decay_sum[1]

## physical distance (bp) at which rho = 0.99 is reached
d_from_rho(a = chr1$a, rho = 0.99)

## expected r^2 at rho = 0.99
ld_from_rho(b = chr1$b, c = chr1$c, rho = 0.99)

# ------------------------------------------------------------
# 5. Clean up
# ------------------------------------------------------------

SNPRelate::snpgdsClose(gds)
unlink(gds_path)
