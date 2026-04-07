library(LDscnR)
devtools::load_all()
devtools::document()
library(ggplot2)
library(data.table)
library(SNPRelate)
#devtools::load_all()

# tmp <- readRDS("../LD-scaling-genome-scans/results_sim/c1_V0.5_rep4.rds")
# map <- tmp$SNP_res[,.(V,c,rep,Chr,Pos,marker,type, Chr_type,max_LD_with_QTN, bp_to_focal_QTN, focal_QTN, p_Va,emx_F,lfmm_F,emx_q,lfmm_q )]
# GTs <- tmp$GTs
#
#
# sim_ex <- list(
#   map = map,
#   GTs = GTs
# )
#
# # Save in package format
# usethis::use_data(sim_ex, overwrite = TRUE)
# usethis::use_r("data_sim_ex")

message("Creating gds file")

## the map we will need later
library(LDscnR)
library(ggplot2)
library(data.table)
library(patchwork)

data("sim_ex")

map <- sim_ex$map

## create gds object used throghout
gds_path <- tempfile(fileext = ".gds")
gds <- create_gds_from_geno(geno=sim_ex$GTs,map, gds_path)

# ------------------------------------------------------------
# 2. LD decay + structure
# ------------------------------------------------------------

# specify number of cores used throughout
cores = 8

## estimate preliminary LD-decay to be able to run full analyses on smaller sliding window
ld_decay_pre <- compute_LD_decay(
  gds,
  n_win_decay = 5, ## Number of windows for estimating LD-decay; default is 20 but since the simulated chromosomes are only ~10cm, only 5 is used here
  max_SNPs_decay = 5000, ## maximum number of SNPs per chromosome for initial decay analyses. For this simulated example, average number of SNPs is 1500 so all are used, but makes a difference for empirical data!
  slide=1000, ## Slide window large because we are analysing a susbset of SNPs only
  cores = cores
)
#ld_decay_pre

## specify new slide window. Here it will rougly same as above, since chromosomes are only ~10cM (ld decays slower)
new_slide_window <- mean(ld_decay_pre$decay_sum$slide_snp_rho_0.99)
ld_decay <- compute_LD_decay(
  gds,
  ## for LD-decay and bg
  n_win_decay = 5,
  max_SNPs_decay = Inf, ## all SNPs are now used...
  keep_el = TRUE, ## now we keep the edge lists for downstream analysse
  slide=2000, ## .. but a smaller window is used
  cores = 8,
  ld_method = "r"
)


## with saving EL
ld_decay <- compute_LD_decay(
  gds,el_data_folder = "./el/",
  ## for LD-decay and bg
  n_win_decay = 5,
  max_SNPs_decay = Inf, ## all SNPs are now used...
  keep_el = TRUE, ## now we keep the edge lists for downstream analysse
  slide=1000, ## .. but a smaller window is used
  cores = 1
)

# ------------------------------------------------------------
# 3. Generate draws
# ------------------------------------------------------------

#pre-calulate ld_w's
ld_ws <- precalculate_ld_w(c(seq(0.75,0.95,by=0.05),0.99),ld_decay)

#plot(ld_ws[,1])

draws <- ld_rho_draws(gds,
                      ld_decay  = ld_decay,
                      F_vals     = map[,.(lfmm_F,emx_F)], ## the F values are pre-calculated and must be provided, these will be used for getting F_prime.
                      q_vals     = NULL, ## q-values from original analyses can also be used but not necessary
                      n_draws    = 100,
                      stat_type  = "q",
                      rho        = NULL, ## grid of rho values used for ld_w
                      ld_ws      = ld_ws,
                      rho_d_lim  = list(min=0.9,max=0.99),
                      rho_ld_lim = list(min=0.9,max=0.99),
                      alpha_lim  = list(min=0.3,max=2), ## lowest alpha is 1/10^0.6=0.25; this range is defined in the -log10 scale!
                      lmin_lim   = list(min=1,max=10),
                      cores      = cores,
                      mode       = "joint"
)

# ------------------------------------------------------------
# 4. Plot Manhattan
# ------------------------------------------------------------

## add log value for LFMM
map[,log_lfmm_q := -log10(lfmm_q)]

plot_manhattan(map,
               gds,
               ld_decay,
               draws,
               ## for defining ORs
               sign_th=0.05,
               mode="joint",
               sign_if ="greater",
               rho_d = 0.9,
               rho_ld = 0.9,
               ## outlier regiosn are defined based on the first in y_vars
               y_vars = c("Joint_C","log_lfmm_q"),
               y_labels = c("C (Joint analyses)","-log10(q) | LFMM"),
               titles = c("a) Joint analyses - C","b) Latent factor mixed model (LFMM) analyses"),
               thresholds = c(0.05,-log10(0.05)),
               col_var = "OR_id")


# ------------------------------------------------------------
# 5. Explore true and false positives in simulated data (using same classification as in the original paper)
# ------------------------------------------------------------

# maximum LD between a SNP and its QTN_focal relative to LD-decay
ld_test_th = 0.75
# maximum physical distance between a SNP and its QTN_focal relative to LD-decay
d_test_th  = 0.95
# Mimumum proportion of variance explaiend by QTN to be considered a true positive
p_Va_th    = 0.05

# Map decay parameters to map
as <- setNames(ld_decay$decay_sum$a, ld_decay$decay_sum$Chr)
bs <- setNames(ld_decay$decay_sum$b, ld_decay$decay_sum$Chr)
cs <- setNames(ld_decay$decay_sum$c, ld_decay$decay_sum$Chr)
C <- cs[map$Chr]
a <- as[map$Chr]
b <- bs[map$Chr]

# add distance and LD with QTN_focal relative to LD-decay
map[,rho_ld := 1 - (max_LD_with_QTN - b) / (C - b)]
map[,rho_d := a*bp_to_focal_QTN/(a*bp_to_focal_QTN+1)]
map[type=="QTN", rho_ld:=1]

# Keep only potential true positive SNPs (those that are correlated and close enough to at least one QTN)
map_filt <- map[rho_ld > ld_test_th & rho_d < d_test_th]

# true positive QTN focal (those that we have a chance to find based on p_Va)
true_QTN_focal <- unique(map_filt[p_Va > p_Va_th, focal_QTN])

# ---- add outlier regions to map data
map_manh <- add_consistency_to_map(map, consistency_obj = consistency_score(draws$draws))
map_manh <- add_ORs(gds, ld_decay, map_manh, stat="Joint_C", sign_th=0.05,sign_if="greater",mode="joint",rho_d=0.99, rho_ld=0.99,l_min=1)

# Unique outlier regions
unique_ORs <- unique(na.omit(map_manh$OR_id))

# Split outlier regions into clusters
ORs <- split(map_manh[OR_id!="ns"]$marker,map_manh[OR_id!="ns"]$OR_id)

# Get QTN_focal for each OR (only the one that has the highest LD with any SNP in that cluster)
OR_focals <- vapply(ORs, function(or) {

  sub <- map_filt[marker %in% or]

  if (nrow(sub) == 0L) return(NA_character_)

  sub[which.max(rho_ld), focal_QTN]

}, character(1))

OR_focals <- unique(na.omit(OR_focals))

# Not all OR focals are true
true_pos_focals <- OR_focals[OR_focals %in% true_QTN_focal]

# These are the true outlier regions - each true positive focal can be counted once (the OR with highest LD with QTN_focal wins)
true_ORs <- sapply(true_pos_focals,function(t_focal){
  map_manh[focal_QTN==t_focal][order(-max_LD_with_QTN)][1,OR_id]
})

# exclude non-significant SNPs
true_ORs <- true_ORs[true_ORs != "ns"]

# Now map true and false positive outlier regiosn to map
map_manh[,OR_test := ifelse(OR_id %in% true_ORs,"True pos","False pos") ]


# plot manhattan
layout <- prep_manhattan(
  map_manh[, .(
    Pos,
    Chr,
    marker,
    OR_test,
    log_lfmm_q,
    Joint_C,
    # type is "QTN" or "ntrl"
    type
  )]
)

plot_manhattan_gg(
  layout,
  y_vars = c("Joint_C","log_lfmm_q"),
  y_labels = c("C (Joint analyses)","-log10(q) | LFMM"),
  titles = c("a) Joint analyses - C","b) Latent factor mixed model (LFMM) analyses"),
  thresholds = c(0.05,-log10(0.05)),
  col_var = "OR_test",
  shape_var = "type",
  ncol=1,
  col_vector=c("firebrick4","steelblue")
)

# ------------------------------------------------------------
# 6. Area under the PR curve for simulated data
# ------------------------------------------------------------

PR_dat <- cbind(method="joint",get_PR(draws$draws$OR, map = map,map_filt, p_Va_th=0.05,cores = cores))

get_AUC_OR(PR_dat)$AUC$AUC_norm

C_scores <- add_consistency_to_map(map,consistency_score(draws$draws))[, .(marker = marker, Joint_C = Joint_C)]


draws_C <- ld_rho_draws(
  gds,
  ld_decay  = ld_decay,
  F_vals    = NULL,
  C_scores  = C_scores[, .(Joint_C)],
  q_vals    = NULL,
  n_draws   = 1000,
  stat_type = "C",
  rho_d_lim = list(min = 0.9, max = 0.99),
  rho_ld_lim = list(min = 0.9, max = 0.99),
  alpha_lim = NULL,
  lmin_lim  = NULL,
  C_lim     = list(min = 0, max = 0.5),
  cores     = cores,
  mode      = "joint"
)

PR_C <- cbind(method="joint",get_PR(draws_C$draws$OR, map = map,map_filt,p_Va_th = 0.05, cores = cores))

get_AUC_OR(PR_C)$AUC$AUC_norm
