## data-raw/stickleback.R
## Builds LDscnR::stickleback -- a real-data LD example from a three-spined
## stickleback (Gasterosteus aculeatus) marine-freshwater ecotype data set,
## chosen to contrast three LD regimes across three chromosomes:
##   * Chr4 ~12.8 Mb : the Eda marine-freshwater adaptation locus (high-LD outlier)
##   * Chr1 ~21   Mb : the chromosome-1 inversion (structural, high long-range LD)
##   * Chr6 ~8    Mb : a background region with no marine-freshwater outliers (low LD)
## (>= 2 chromosomes are required for compute_LD_decay()'s inter-chromosomal
## background-LD estimate; the three regimes also make the per-chromosome decay
## differences visible.)
##
## SOURCE: Fang B, Kemppainen P, Momigliano P, Merila J (2021). "Population
##   structure limits parallel evolution in sticklebacks." Molecular Biology and
##   Evolution 38(10):msab144. doi:10.1093/molbev/msab144.
##   Three-spined stickleback (Gasterosteus aculeatus) marine-freshwater ecotype.
##   ANGSD genotype dosages in [0, 2] (no missing data), 117 individuals
##   (27 Marine, 90 Freshwater). Inversion/background positions located from a
##   local-LD scan of this data set (Chr1 20-22 Mb mean r2 ~0.079 vs ~0.017
##   background; Chr6 flat at background).
##
## SUBSET: each region is all MAF > 0.1 markers within +/- 1.25 Mb of the centre
##   (~15,900 markers total). Re-running requires the source .RData, which is NOT
##   distributed with the package (data-raw/ is .Rbuildignore'd).

library(data.table)

SRC <- path.expand("~/gitlab/LD-scaling-genome-scans/empirical_data/3sp/3sp_data.RData")
e <- new.env(); load(SRC, envir = e)          # GTs_3sp, map_3sp, pheno_3sp
GTs <- e$GTs_3sp
map <- as.data.table(e$map_3sp)
ph  <- as.data.table(e$pheno_3sp)

regions <- list(
  list(chr = "Chr1", center = 21.0e6),   # inversion
  list(chr = "Chr4", center = 12.8e6),   # Eda
  list(chr = "Chr6", center =  8.0e6)    # background
)
halfwidth <- 1.25e6

## marker indices per region, each ordered by position; regions concatenated
sel_idx <- unlist(lapply(regions, function(r) {
  ii <- which(map$Chr == r$chr & map$maf > 0.1 &
              abs(map$Pos - r$center) <= halfwidth)
  ii[order(map$Pos[ii])]
}))

GTs_sub <- GTs[, sel_idx]
map_sub <- map[sel_idx, .(Chr, Pos, marker, maf)]
## store Chr as character (a factor Chr is coerced to its integer level code
## when a list is indexed as by_chr[[Chr]], which breaks ld_complexity_reduction)
map_sub[, Chr := as.character(Chr)]
dimnames(GTs_sub) <- list(ph$ID, map_sub$marker)

pheno_sub <- ph[, .(ID, ecotype = as.character(ecotype), pop_ID, lineage, pop_locality)]

stickleback <- list(GTs = GTs_sub, map = map_sub, pheno = pheno_sub)

usethis::use_data(stickleback, overwrite = TRUE, compress = "xz")
