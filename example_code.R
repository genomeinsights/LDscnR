#library(LDscnR)
library(ggplot2)
library(data.table)
library(SNPRelate)

devtools::load_all()


tmp <- readRDS("../LD-scaling-genome-scans/results_sim/c1_V0.5_rep4.rds")
map_org <- tmp$SNP_res[,.(V,c,rep,Chr,Pos,marker,type, Chr_type,max_LD_with_QTN, bp_to_focal_QTN, focal_QTN, p_Va,emx_F,lfmm_F,emx_q,lfmm_q,a,b,c )]
map <- tmp$SNP_res[,.(V,c,rep,Chr,Pos,marker,type, Chr_type,max_LD_with_QTN, bp_to_focal_QTN, focal_QTN, p_Va,emx_F,lfmm_F,emx_q,lfmm_q )]
geno <- tmp$GTs




#data(sim_ex)

#geno = sim_ex$GTs
#map = sim_ex$map

## generate ld_w for draws

#
# rho_w     = 0.9
# full = TRUE
# n_rep = 10
# q = 0.95
# n_sub = 5000
# window_size = 1e6
# step_size = 5e5
# dist_unit = 5000
# K_target = 1000
# n_win = 1000
# prob_robust = 0.5
# use="robust"
# n_rho = 40
# n_or_draws = 25
# rho_w_lim = list(min=0.5,max=0.95)
# rho_d_lim=list(min=0.9,max=0.999)
# rho_ld_lim=list(min=0.5,max=0.999)
# alpha_lim=list(min=1.31,max=4)
# lmin_lim=list(min=1,max=10)
# cores=8
# n_inds = 125
# ## output
# verbose=TRUE
#
# max_rho = 0.9
# t1 <- Sys.time()
#
# SNP_ids <- map$marker
# n_inds <- nrow(geno)

message("Creating gds file")

gds_path <- tempfile(fileext = ".gds")
gds <- create_gds_from_geno(geno, map, gds_path)
#on.exit({ SNPRelate::snpgdsClose(gds); unlink(gds_path) }, add = TRUE)

# ------------------------------------------------------------
# 2. LD decay + structure
# ------------------------------------------------------------
t1 <- Sys.time()
ld_struct <- compute_ld_structure(
  gds,
  ## for LD-decay and bg
  q = 0.95,
  ## for bg
  n_sub_bg = 5000,
  ## for decay
  n_win_decay = 20,
  max_pairs = 5000,
  n_strata = 20,
  overlap = 0.5,
  prob_robust = 0.95,
  target_dist_bins_for_decay = 60,
  keep_el = TRUE,
  n_bins_ld_int = 20,
  rho = 0.995,
  k_max=1000,
  cores = 8
)
t2 <- Sys.time()
difftime(t2,t1)
#plot(-log(ld_struct$decay_sum$a),map[,max(Pos),by=Chr][,V1])
#rm(ld_struct)



map[,ld_int:=unlist(lapply(ld_struct$by_chr,function(x)x$ld_int))]

map[,cor(ld_int*lfmm_F,max_LD_with_QTN)^2]

lfmm_F_prime_full <- .compute_Fprime(map$lfmm_F,map$ld_int,n_rep=10,n_inds = nrow(geno),full=TRUE)
map[,lfmm_q_prime := -log10(lfmm_F_prime_full[[1]]$q_prime)]

p1 <- ggplot(map, aes(ld_int,lfmm_F,col=max_LD_with_QTN,shape=type,size=type)) +
  geom_point() +
  scale_color_viridis_c(option="turbo")+
  scale_shape_manual(values=c(20,3)) +
  scale_size_manual(values=c(1,3))

p2 <- ggplot(map, aes(ld_int,lfmm_F,col=(ld_int*lfmm_F)^0.5,shape=type,size=type)) +
  geom_point() +
  scale_color_viridis_c(option="turbo") +
  scale_shape_manual(values=c(20,3)) +
  scale_size_manual(values=c(1,3))

p3 <- ggplot(map, aes(1:nrow(map),lfmm_q_prime,col=max_LD_with_QTN,shape=type,size=type)) +
  geom_point() +
  scale_color_viridis_c(option="turbo") +
  scale_shape_manual(values=c(20,3)) +
  scale_size_manual(values=c(1,3)) +
  geom_hline(yintercept = 1.31)

p4 <- ggplot(map, aes(1:nrow(map),lfmm_F,col=max_LD_with_QTN,shape=type,size=type)) +
  geom_point() +
  scale_color_viridis_c(option="turbo") +
  scale_shape_manual(values=c(20,3)) +
  scale_size_manual(values=c(1,3))

p5 <- ggplot(map, aes(1:nrow(map),ld_int,col=max_LD_with_QTN,shape=type,size=type)) +
  geom_point() +
  scale_color_viridis_c(option="turbo") +
  scale_shape_manual(values=c(20,3)) +
  scale_size_manual(values=c(1,3))

(p1 | p2) / p3 / p4 / p5


as <- setNames(ld_struct$decay_sum$a,ld_struct$decay_sum$Chr)
bs <- setNames(ld_struct$decay_sum$b,ld_struct$decay_sum$Chr)
cs <- setNames(ld_struct$decay_sum$c, ld_struct$decay_sum$Chr)
C <- cs[map$Chr]
a <- as[map$Chr]
b <- bs[map$Chr]

map[,rho_ld := 1 - (max_LD_with_QTN - b) / (C - b)]
map[,rho_d := a*bp_to_focal_QTN/(a*bp_to_focal_QTN+1)]
map[type=="QTN", rho_ld:=1]


F_cols=c("emx_F","lfmm_F")
q_cols=c("emx_q","lfmm_q")

draws_ld_w <- ld_rho_draws(gds,
                      ld_struct  = ld_struct,
                      F_vals     = map[,..F_cols],
                      q_vals     = NULL,
                      n_rho_w    = 25,
                      n_draws    = 25,
                      ld_w_int   = NULL,
                      stat_type  = "q",
                      rho_w_lim  = list(min=0.75,max=0.9999),
                      rho_d_lim  = list(min=0.5,max=0.99),
                      rho_ld_lim = list(min=0.5,max=0.99),
                      alpha_lim  = list(min=1.31,max=4),
                      lmin_lim   = list(min=1,max=10),
                      cores      = 8,
                      mode       = c("joint","per_method")
)

draws_ld_int <- ld_rho_draws(gds,
                           ld_struct  = ld_struct,
                           F_vals     = map[,..F_cols],
                           q_vals     = NULL,
                           n_rho_w    = 1,
                           n_draws    = 25*25,
                           ld_w_int   = map$ld_int,
                           stat_type  = "q",
                           rho_w_lim  = NULL,
                           rho_d_lim  = list(min=0.5,max=0.99),
                           rho_ld_lim = list(min=0.5,max=0.99),
                           alpha_lim  = list(min=1.31,max=4),
                           lmin_lim   = list(min=1,max=10),
                           cores      = 8,
                           mode       = c("joint","per_method")
)


ld_test_th = 0.75
d_test_th  = 0.95
p_Va_th    = 0.05
map_filt <- map[rho_ld > ld_test_th & rho_d < d_test_th]

PR_ld_w <- get_PR(draws_ld_w$draws$OR,map_filt)
PR_ld_w <- cbind(map[1,.(c,V,rep)],draws_ld_w$draws[,.(method, rho_w, rho_d, rho_ld, alpha, l_min,n_ORs=OR_size)],PR_ld_w)


PR_int <- get_PR(draws_ld_int$draws$OR,map_filt)
PR_int <- cbind(map[1,.(c,V,rep)],draws_ld_int$draws[,.(method, rho_w, rho_d, rho_ld, alpha, l_min,n_ORs=OR_size)],PR_int)
#PR_int[method=="Joint",method:="Joint_int"]

PR_int[,ld:="int"]
PR_int[,method := gsub("_int","",method)]
PR_ld_w[,ld:="rho_w"]
PR_data <- rbind(PR_ld_w,PR_int)

PR_data[,median(PR)/mad(PR, constant = 1),by=.(method,ld)]

ggplot(PR_data,
       aes(x = recall,
           y = precision,
           colour = method)) +
  facet_grid(ld~.)+
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE) +
  theme_bw() +
  labs(title = "Precision–Recall tradeoff across parameter draws")

ggplot(PR_data,
       aes(x = method,
           y = PR,
           fill = ld)) +
  #facet_grid(ld~.)+
  #geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distribution of PR score across parameter draws")


#map2[,plot(Joint_int_C)]

map2 <- add_consistency_to_map(map, consistency_obj = consistency_score(draws_ld_int$draws))


p1 <- ggplot(map2, aes(1:nrow(map),Joint_C,col=max_LD_with_QTN,shape=type,size=type)) +
  geom_point() +
  scale_color_viridis_c(option="turbo") +
  scale_shape_manual(values=c(20,3)) +
  scale_size_manual(values=c(1,3)) +
  geom_hline(yintercept = 0.03)

map2 <- add_consistency_to_map(map, consistency_obj = consistency_score(draws_ld_w$draws))


p2 <- ggplot(map2, aes(1:nrow(map),Joint_C,col=max_LD_with_QTN,shape=type,size=type)) +
  geom_point() +
  scale_color_viridis_c(option="turbo") +
  scale_shape_manual(values=c(20,3)) +
  scale_size_manual(values=c(1,3)) +
  geom_hline(yintercept = 0.03)

p1 / p2

library(dplyr)
library(tidyr)
#PR_data[,length(table(unique(paste0(rho_w,rho_d,rho_ld,alpha,l_min))))/length(unique(method))]

best <- PR_data %>%
  group_by(rho_d, rho_ld, alpha, l_min) %>%
  slice_max(PR, n = 1, with_ties = FALSE) %>%
  ungroup()
#table(best$method)

prop_df <- best %>%
  count(method,ld) %>%
  mutate(prop = n / sum(n))
# prop_df <- prop_df %>%
#   mutate(prop = wins / total_draws)

ggplot(prop_df, aes(method, prop, fill = ld)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    y = "Proportion of parameter draws won",
    x = "Method",
    title = "Dominance across OR parameter space"
  )


better_df <- PR_data %>%
  tidyr::pivot_wider(
    id_cols = c(rho_w, rho_d, rho_ld, alpha, l_min),
    names_from = method,
    values_from = PR
  ) %>%
  mutate(better = lfmm_F_prime > lfmm_q)

mean(better_df$better)


ggplot(prop_df, aes(x = method, y = prop, fill = method)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw() +
  labs(
    y = "Proportion of draws won",
    x = "Method",
    title = "Dominance frequency across parameter draws"
  ) +
  theme(legend.position = "none")


ggplot(best, aes(method)) +
  geom_bar() +
  theme_bw() +
  labs(title = "Which method wins per draw?")

#nrow(draws_joint$draws[[1]])
PR <- rbindlist(lapply(draws_joint$draws,function(draws){
  rbindlist(parallel::mclapply(seq_len(nrow(draws)),function(draw){
    get_PR(draws[draw,],map = map)
  },mc.cores=cores))
}))
t2 <- Sys.time()
difftime(t2,t1)

tmp <- do.call(rbind,strsplit(PR$rho_w,"_"))
PR[,stat:=tmp[,1]]
PR[,rho_w:=as.numeric(tmp[,3])]

saveRDS("./PR_out/")PR

#PR[,plot(-log10(alpha),PR)]
#PR[,plot(l_min,PR)]

ggplot(PR,aes(factor(rho_w),PR,fill=stat)) +
  geom_violin() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))

PR[which.max(PR)]

draws_joint$draws$joint_ldw <- rbindlist(draws_joint$draws[grep("ldw",names(draws_joint$draws))])
draws_joint$draws$joint_ldexc <- rbindlist(draws_joint$draws[grep("ldexc",names(draws_joint$draws))])


draws <- draws_joint$draws[[1]]
x <- 1
rho_name <- names(draws_joint$draws)[1]

draws_joint$draws$joint <- rbindlist(draws_joint$draws)

plots <- lapply(names(draws_joint$draws)[],function(rho_name){
  map2 <- add_consistency_to_map(map, consistency_obj = consistency_score(draws_joint$draws[[rho_name]]))
  #map2 <- add_consistency_to_map(map2, consistency_obj = consistency_score(draws_per_method))
  #ldw = unlist(ld_struct$ld_w_dt[,..rho_name])
  #map2[,ld_w:=ldw]
  #map2[,plot(LD_AUC)]
  #plot(unlist(lapply(ld_struct$by_chr, function(x) x$LD_AUC$LD_AUC)))
  #map2[,plot(ld_w,max_LD_with_QTN )]
  #map2[,plot(Joint_C)]

  #map2[!is.na(Joint_C)]

  #map2[!is.na(Joint_C),plot(Joint_C)]
  sign_th = 0.01
  idx <- which(map2[,Joint_C>sign_th])
  el  <- get_el(gds, idx, slide_win_ld = -1,by_chr = TRUE)


  ORs_tbl <- detect_or(el,
                       q_vals=map2[,.(Joint_C)],
                       SNP_ids = map$marker,
                       SNP_chr = map$Chr,
                       ld_struct=ld_struct,
                       sign_th=sign_th,
                       sign_if = "greater",
                       rho_d=0.99,
                       rho_ld=0.99,
                       l_min = 1,
                       ret_table = TRUE)


  #ORs_tbl[,table(do.call(rbind, strsplit(SNP,":"))[,1],OR_id)]
  map2 <- merge(map2,
                ORs_tbl,
                by.x = "marker",
                by.y = "SNP",
                all.x = TRUE)
  map2 <-  map2[order(as.numeric(gsub("Chr","",map2$Chr), map2$Pos))]
  map2[is.na(OR_id),OR_id:="ns"]
  map2[,table(OR_id,Chr)]

  unique_ORs <- unique(na.omit(map2$OR_id))
  #unique_ORs <- unique_ORs[unique_ORs != "Chr12_OR10"]

  map2[, OR_factor := factor(OR_id, levels = sample(unique_ORs))]

  map2[!is.na(OR_factor),.(max(max_LD_with_QTN),.N,focal_QTN[which.max(max_LD_with_QTN)]),by=OR_factor]

  layout <- prep_manhattan(
    map2[, .(
      bp = Pos,
      Chr,
      marker,
      Joint_C,
      lfmm_q = -log10(lfmm_q),
      ld_w,
      OR_factor,
      type
    )],spacer = 1e6
  )


  plot_manhattan_gg(
    layout,
    y_vars = c("Joint_C"),
    y_labels = c("Joint_C"),
    thresholds = c(0.05),
    or_var = "OR_factor",
    type_var = "type",
    ncol=1
  )+ggplot2::ggtitle(rho_name)
})

wrap_plots(plots, ncol = 3)

(plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]]) /
  (plots[[5]] | plots[[6]] | plots[[7]])

(plots[[8]] | plots[[9]] | plots[[10]] | plots[[11]]) /
  (plots[[12]] | plots[[13]] | plots[[14]])

map2[OR_factor=="Chr13_OR10"]


#### compression #####
build_2d_ld_hist <- function(el, focal,
                             dist_unit = 1e6,
                             r2_unit   = 0.001) {

  library(data.table)

  # 1️⃣ Symmetric neighbor extraction
  dt <- rbind(
    el[SNP1 == focal,.(pos = pos1, pos_other = pos2, r2)],
    el[SNP2 == focal,.(pos = pos2, pos_other = pos1, r2)]
  )

  if (nrow(dt) == 0)
    return(NULL)

  # 2️⃣ Distance + side
  dt[, d := abs(pos_other - pos)]
  dt[, side := ifelse(pos_other > pos, "R", "L")]

  # 3️⃣ Compute maxL / maxR
  maxL <- if (any(dt$side == "L")) max(dt$d[dt$side == "L"]) else 0
  maxR <- if (any(dt$side == "R")) max(dt$d[dt$side == "R"]) else 0

  S <- min(maxL, maxR)

  # 4️⃣ Identify long side
  if (maxL > maxR) {
    long_side <- "L"
  } else if (maxR > maxL) {
    long_side <- "R"
  } else {
    long_side <- NA_character_
  }

  # 5️⃣ Reweight: duplicate tail beyond S on long side
  if (!is.na(long_side) && S > 0) {

    tail_dt <- dt[side == long_side & d > S]

    if (nrow(tail_dt) > 0)
      dt <- rbind(dt, tail_dt)
  }

  # 6️⃣ Bin distance + r2
  dt[, dist_bin := floor(d / dist_unit)]
  dt[, r2_bin   := round(r2 / r2_unit) * r2_unit]

  # 7️⃣ Build 2D histogram
  hist2d <- dt[, .N, by = .(dist_bin, r2_bin)]

  hist_mat <- dcast(hist2d,
                    dist_bin ~ r2_bin,
                    value.var = "N",
                    fill = 0)

  mat <- as.matrix(hist_mat[, -1])
  rownames(mat) <- hist_mat$dist_bin

  # Ensure sorted rows
  ord <- order(as.numeric(rownames(mat)))
  mat <- mat[ord, , drop = FALSE]

  # 8️⃣ Convert to cumulative along distance
  cum_mat <- apply(mat, 2, cumsum)

  rownames(cum_mat) <- rownames(mat)
  colnames(cum_mat) <- colnames(mat)

  return(cum_mat)

}


cum_histograms2d <- lapply(
  ld_struct_w1000$by_chr$Chr1$snp_ids,
  function(focal) {
    build_cum_ld_hist2d_reweighted(
      el = el,
      focal = focal,
      dist_unit = 1e6,
      r2_unit = 0.001
    )
  }
)

cum_hist2d <- cum_histograms2d[[1]]
ld_median_from_hist <- function(cum_hist2d, W_bin) {

  dist_bins <- as.numeric(rownames(cum_hist2d))
  idx <- max(which(dist_bins <= W_bin))
  if (length(idx) == 0 || is.infinite(idx))
    return(NA_real_)

  r2_counts <- cum_hist2d[idx, ]
  total_n <- sum(r2_counts)
  if (total_n == 0)
    return(NA_real_)

  cum_r2 <- cumsum(r2_counts)
  median_pos <- ceiling(total_n / 2)
  j <- which(cum_r2 >= median_pos)[1]

  as.numeric(colnames(cum_hist2d))[j]
}

ld_w_from_hist <- sapply(cum_histograms2d,ld_median_from_cum,22)

ld_w_full <- .unbiased_ld_w(ld_struct_w1000$by_chr$Chr1$edges)

res <- data.table::data.table(
  SNP = ld_struct_w1000$by_chr$Chr1$snp_ids
)
ld_w_full <- ld_w_full[res, on="SNP"]
cor(ld_w_from_hist,ld_w_full$ld_w)^2



build_cum_ld_hist2d <- function(dt,
                                dist_unit = 1e6,
                                r2_unit   = 0.05) {

  library(data.table)
  dt <- copy(dt)

  # Bin distance and r2
  dt[, dist_bin := floor(d / dist_unit)]
  dt[, r2_bin   := round(r2 / r2_unit) * r2_unit]

  # Count
  hist2d <- dt[, .N, by = .(dist_bin, r2_bin)]

  # Cast to wide
  hist_mat <- dcast(hist2d,
                    dist_bin ~ r2_bin,
                    value.var = "N",
                    fill = 0)

  # Convert to matrix
  mat <- as.matrix(hist_mat[, -1])
  rownames(mat) <- hist_mat$dist_bin

  # Ensure rows sorted by distance
  ord <- order(as.numeric(rownames(mat)))
  mat <- mat[ord, , drop = FALSE]

  # Convert to cumulative along distance
  cum_mat <- apply(mat, 2, cumsum)

  # Preserve rownames and colnames
  rownames(cum_mat) <- rownames(mat)
  colnames(cum_mat) <- colnames(mat)

  return(cum_mat)
}

W_bin <- 4
cum_hist2d <- cum_histograms2d[[1]]
ld_median_from_cum <- function(cum_hist2d, W_bin) {

  dist_bins <- as.numeric(rownames(cum_hist2d))

  idx <- max(which(dist_bins <= W_bin))
  if (length(idx) == 0 || is.infinite(idx))
    return(NA_real_)

  r2_counts <- cum_hist2d[idx, ]

  total_n <- sum(r2_counts)
  if (total_n == 0)
    return(NA_real_)

  cum_r2 <- cumsum(r2_counts)
  median_pos <- ceiling(total_n / 2)

  j <- which(cum_r2 >= median_pos)[1]

  r2_vals <- as.numeric(colnames(cum_hist2d))

  r2_vals[j]
}

cum_histograms2d <- lapply(ld_struct_w1000$by_chr$Chr1$snp_ids,function(focal){
  dt <- rbind(
    el[SNP1 == focal, .(d, r2)],
    el[SNP2 == focal, .(d, r2)]
  )

  dt[, d := abs(pos_other - pos)]
  dt[, side := ifelse(pos_other > pos, "R", "L")]

  dt[, `:=`(
    maxL = if (any(side == "L")) max(d[side == "L"]) else 0,
    maxR = if (any(side == "R")) max(d[side == "R"]) else 0
  ), by = SNP]

  dt[, S := pmin(maxL, maxR)]

  dt[, long_side := data.table::fifelse(
    maxL > maxR, "L",
    data.table::fifelse(maxR > maxL, "R", NA_character_)
  )]

  dt[, tail := (side == long_side & d > S)]
  dt[is.na(tail), tail := FALSE]

  dt <- data.table::rbindlist(list(
    dt[, .(SNP, r2)],
    dt[tail == TRUE, .(SNP, r2)]
  ))

  build_cum_ld_hist2d(dt, dist_unit = 1e6,r2_unit = 0.001)
})
ld_w_from_hist <- sapply(cum_histograms2d,ld_median_from_cum,W_bin)


