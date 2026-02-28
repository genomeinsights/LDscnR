#library(LDscnR)
library(ggplot2)
library(data.table)
library(SNPRelate)

tmp <- readRDS("../LD-scaling-genome-scans/results_sim/c1_V0.5_rep4.rds")
map_org <- tmp$SNP_res[,.(V,c,rep,Chr,Pos,marker,type, Chr_type,max_LD_with_QTN, bp_to_focal_QTN, focal_QTN, p_Va,emx_F,lfmm_F,emx_q,lfmm_q,a,b,c )]
map <- tmp$SNP_res[,.(V,c,rep,Chr,Pos,marker,type, Chr_type,max_LD_with_QTN, bp_to_focal_QTN, focal_QTN, p_Va,emx_F,lfmm_F,emx_q,lfmm_q )]
geno <- tmp$GTs
data(sim_ex)

#geno = sim_ex$GTs
#map = sim_ex$map

## generate ld_w for draws


F_cols=c("emx_F","lfmm_F")
q_cols=c("emx_q","lfmm_q")
F_vals     = map[,..F_cols]
q_vals     = map[,..q_cols]
rho_w     = 0.9
full = TRUE
n_rep = 10
q = 0.95
n_sub = 5000
window_size = 1e6
step_size = 5e5
dist_unit = 5000
K_target = 1000
n_win = 1000
prob_robust = 0.5
use="robust"
n_rho = 40
n_or_draws = 25
rho_w_lim = list(min=0.5,max=0.95)
rho_d_lim=list(min=0.9,max=0.999)
rho_ld_lim=list(min=0.5,max=0.999)
alpha_lim=list(min=1.31,max=4)
lmin_lim=list(min=1,max=10)
cores=8
n_inds = 125
## output
verbose=TRUE

max_rho = 0.9
t1 <- Sys.time()

SNP_ids <- map$marker
n_inds <- nrow(geno)

message("Creating gds file")

gds_path <- tempfile(fileext = ".gds")
gds <- create_gds_from_geno(geno, map, gds_path)
#on.exit({ SNPRelate::snpgdsClose(gds); unlink(gds_path) }, add = TRUE)

# ------------------------------------------------------------
# 2. LD decay + structure
# ------------------------------------------------------------
t1 <- Sys.time()
ld_struct <-  compute_ld_structure(
  gds,
  ## for LD-decay and bg
  q = 0.95,
  ## for bg
  n_sub_bg = 5000,
  ## for decay
  n_win_decay = 20,
  overlap = 0.5,
  prob_robust = 0.95,
  target_dist_bins_for_decay = 40,
  n_snps_for_decay = 500,
  ## for histograms
  n_dist_target_for_hist = 100,
  eps = 0.005,
  r2_unit = 0.001,
  ## cores
  cores = 1
)

#ld_struct$histograms$Chr1
#plot(ld_struct)
t2 <- Sys.time()
difftime(t2,t1)
ld_struct

cors <- rbindlist(lapply(seq(0.5,1,by=0.05)-0.000001,function(rho_w){
  d_star <- derive_ld_radius(ld_struct$by_chr$Chr1$decay_sum$a, 1-rho_w)
  ld_w <- compute_ld_summary(ld_structure=ld_struct,
                             method = c("rho_w"),
                             eps = 0.005,
                             d_window = d_star,
                             shell_type = "median")

  data.table(rho_w=rho_w,
             cor_ld=cor(ld_w,map$max_LD_with_QTN,use="pair")^2,
             cor_d=cor(ld_w,map$bp_to_focal_QTN,use="pair")^2 ,
             mean_QTN= mean((ld_w * map$lfmm_F)[map$type=="QTN"],na.rm=TRUE)/mean((ld_w * map$lfmm_F)[map$type!="QTN"],na.rm=TRUE))
}))

# LD_int <- compute_ld_summary(ld_structure=ld_struct,
#                              method = c("ld_int"),
#                              eps = 0.005,
#                              d_window = NULL,
#                              shell_type = "excess")

par(mfcol=c(1,3))
plot(seq(0.5,1,by=0.05),cors$cor_d,type="l")
abline(h=cor(LD_int_median,map$bp_to_focal_QTN)^2)

plot(seq(0.5,1,by=0.05),cors$cor_ld,type="l")
abline(h=cor(LD_int_median,map$max_LD_with_QTN)^2)

plot(seq(0.5,1,by=0.05),cors$mean_QTN,type="l")
abline(h=mean((LD_int_median * map$lfmm_F)[map$type=="QTN"])/mean((LD_int_median * map$lfmm_F)[map$type!="QTN"]))

plot(seq(0.5,1,by=0.05),cors$cor_ld,type="l")
abline(h=mean((LD_int_median * map$lfmm_F)[map$type=="QTN"])/mean((LD_int * map$lfmm_F)[map$type!="QTN"]))


plot(cors$cor_ld,cors$cor_d)

cor(LD_int_median,map$lfmm_F)^2
cor(LD_int,map$lfmm_F)^2
LD_int_median <- compute_ld_summary(ld_structure=ld_struct,
                                     method = c("ld_int"),
                                     eps = 0.005,
                                     d_window = NULL,
                                     shell_type = "median",
                                    cores = 8)
cor(LD_int_mean,map$max_LD_with_QTN)^2
cor(LD_int_median,map$max_LD_with_QTN)^2
cor(LD_int_median_int,map$max_LD_with_QTN)^2
cor(LD_int,map$max_LD_with_QTN)^2

cor(LD_int_median, LD_int)^2
sd(LD_int_median)
sd(LD_int)

ld_w <- compute_ld_summary(ld_structure=ld_struct,
                           method = c("rho_w"),
                           eps = 0.005,
                           d_window = derive_ld_radius(ld_struct$by_chr$Chr1$decay_sum$a, 1-0.01),
                           shell_type = "excess")

par(mfcol=c(2,1))
plot(LD_int*map$lfmm_F,col=ifelse(map$type=="QTN","red","black"),pch=ifelse(map$type=="QTN",3,16),cex=ifelse(map$type=="QTN",3,1))
plot(ld_w*map$lfmm_F,col=ifelse(map$type=="QTN","red","black"),pch=ifelse(map$type=="QTN",3,16),cex=ifelse(map$type=="QTN",3,1))

#plot(LD_int*map$lfmm_F,ld_w*map$lfmm_F,pch=ifelse(map$type=="QTN",3,16),cex=ifelse(map$type=="QTN",3,1))

mean((LD_int * map$lfmm_F)[map$type=="QTN"])/mean((LD_int * map$lfmm_F)[map$type!="QTN"])
mean((ld_w * map$lfmm_F)[map$type=="QTN"])/mean((ld_w * map$lfmm_F)[map$type!="QTN"])

tapply(LD_int * map$lfmm_F, map$type, mean)
tapply(ld_w  * map$lfmm_F, map$type, mean)

ld_w <- compute_ld_summary(ld_structure=ld_struct,
                           method = c("rho_w"),
                           eps = 0.005,
                           d_window = derive_ld_radius(ld_struct$by_chr$Chr1$decay_sum$a, 1-0.9),
                           shell_type = "median")

dt_plot <- rbindlist(list(
  data.table(
    index = seq_along(LD_int),
    value = LD_int_median,
    type  = map$type,
    method = "LD_int"
  ),
  data.table(
    index = seq_along(ld_w),
    value = ld_w,
    type  = map$type,
    method = "rho_w"
  )
))

ggplot(dt_plot, aes(x = index, y = value)) +
  geom_point(aes(color = type,
                 shape = type,
                 size = type),
             alpha = 0.8) +
  scale_color_manual(values = c("neutral" = "black",
                                "QTN" = "red")) +
  scale_shape_manual(values = c(16,3)) +
  scale_size_manual(values = c(1,3)) +
  facet_wrap(~method, ncol = 1, scales = "free_y") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "SNP index",
    y = "LD-scaled LFMM statistic",
    color = "SNP type",
    shape = "SNP type",
    size  = "SNP type"
  )



ld_Ws <- do.call(cbind,lapply(seq(0.1,1,by=0.1)-0.000001,function(rho_w){
  d_star <- derive_ld_radius(ld_struct$by_chr$Chr1$decay_sum$a, 1-rho_w)
  ld_w <- compute_ld_summary(ld_structure=ld_struct,
                             method = c("rho_w"),
                             eps = 0.005,
                             d_window = d_star,
                             shell_type = "median")


}))

mean(cor(ld_Ws,use="pair")^2)
mean(cor(ld_Ws*map$lfmm_F,use="pair")^2)

pairs((ld_Ws)[sample(1:nrow(ld_Ws),10000),])

pairs(ld_Ws*map$lfmm_F[1:10000,])

plot(ld_w*map$lfmm_F)

dt <- data.table(LD_int,ld_w,max_LD_with_QTN=map$max_LD_with_QTN,bp_to_focal_QTN=map$bp_to_focal_QTN)

ggplot(dt,aes(LD_int,ld_w,col=max_LD_with_QTN)) +
  geom_point() +
  scale_color_viridis_c(option="turbo")

ggplot(dt,aes(LD_int,ld_w,col=bp_to_focal_QTN)) +
  geom_point() +
  scale_color_viridis_c(option="turbo")


cor(LD_int,map$bp_to_focal_QTN)^2
cor(ld_w,map$bp_to_focal_QTN)^2

abline(h=cor(LD_int,map$max_LD_with_QTN,use="pair")^2,lty=2,col="red")ld_w
plot(LD_int)

a <- ld_struct$decay_summary[Chr=="Chr1",a]

rho_max <- 0.5

d_min <- d_from_rho(a, rho_min)
d_max <- d_from_rho(a, rho_max)

par(mfcol=c(4,2))

LD_star <- do.call(cbind,lapply(seq(0.25,0.95,by=0.1),function(rho_max){
  LD_star <- sapply(ld_struct$histograms$Chr1,function(hist_mat){
    dist_bin <- as.numeric(rownames(hist_mat))
    keep_rows <- dist_bin <= d_from_rho(a, rho_max)
    shell_mat <- hist_mat <- hist_mat[keep_rows, , drop = FALSE]
    shell_mat[-1, ] <- hist_mat[-1, ] - hist_mat[-nrow(hist_mat), ]
    r2_vals <- as.numeric(colnames(shell_mat))
    excess <- pmax(r2_vals - b, 0)
    sum(shell_mat * matrix(excess, nrow=nrow(shell_mat),ncol=ncol(shell_mat),byrow=TRUE)) / sum(shell_mat)
  })
  plot(LD_star,main=paste("rho_max =",rho_max),pch=20)
  return(LD_star)
}))

colnames(LD_rho) <- seq(0.25,0.95,by=0.1)
pairs(LD_star)

LD_rho <- do.call(cbind,
                  lapply(seq(0.25, 0.95, by = 0.001), function(rho_w){

                    LD_val <- sapply(ld_struct$histograms$Chr1, function(hist_mat){

                      if (is.null(hist_mat)) return(NA_real_)

                      dist_bin <- as.numeric(rownames(hist_mat))

                      # distance threshold for this rho
                      d_th <- d_from_rho(a, rho_w)

                      keep_rows <- dist_bin <= d_th
                      if (!any(keep_rows)) return(NA_real_)

                      sub_mat <- hist_mat[keep_rows, , drop = FALSE]

                      # convert cumulative to shell counts
                      shell_mat <- sub_mat
                      if (nrow(shell_mat) > 1) {
                        shell_mat[-1, ] <- sub_mat[-1, ] - sub_mat[-nrow(sub_mat), ]
                      }

                      r2_vals <- as.numeric(colnames(shell_mat))
                      excess  <- pmax(r2_vals - b, 0)

                      weighted_sum <- sum(
                        shell_mat *
                          matrix(excess,
                                 nrow = nrow(shell_mat),
                                 ncol = ncol(shell_mat),
                                 byrow = TRUE)
                      )

                      weighted_sum / sum(shell_mat)
                    })

                    # plot(LD_val,
                    #      main = paste("rho_w =", rho_w),
                    #      pch = 20)

                    LD_val
                  })
)
LD_rho
colnames(LD_rho) <- seq(0.25, 0.95, by = 0.01)
pairs(LD_rho)


cor(LD_rho[,1],LD_star[,1],use = "pair")^2

plot(LD_star_075_075)
par(mfcol=c(4,1))
plot(LD_star_0_1,pch=20)
plot(LD_star_05_1,pch=20)
plot(LD_star_09_1,pch=20)
plot(LD_star_05_075,pch=20)

plot(LD_star_0_1,LD_star_0_075)

abline(0,1)

shell_mat

par(mfcol=c(4,1))
plot(ld_struct$ld_w_dt$ldw_rho_0.7)
plot(ld_struct$ld_w_dt$ldw_rho_0.95)
plot(ld_struct$ld_w_dt$ldexc_rho_0.7)
plot(ld_struct$ld_w_dt$ldexc_rho_0.95)

plot(ld_struct$ld_w_dt$ldw_rho_0.7,ld_struct$ld_w_dt$ldw_rho_0.95)
plot(ld_struct$ld_w_dt$ldexc_rho_0.7,ld_struct$ld_w_dt$ldexc_rho_0.95)


draws_joint <- ld_rho_draws(gds,
                            ld_struct  = ld_struct,
                            F_vals     = map[,..F_cols],
                            q_vals     = map[,..q_cols],
                            n_inds     = nrow(geno),
                            n_draws    = 1000,
                            rho_d_lim  = list(min=0.9,max=0.99),
                            rho_ld_lim = list(min=0.5,max=0.99),
                            alpha_lim  = list(min=1.31,max=4),
                            lmin_lim   = list(min=1,max=10),
                            cores      = 8,
                            mode       = "joint"

)
#draws_joint$draws
#or_dt <- draws_joint$
#draws_joint$draws[2,OR]
#draws_joint$draws
## add decay info to map

as <- setNames(ld_struct$decay_summary$a,ld_struct$decay_summary$Chr)
bs <- setNames(ld_struct$decay_summary$b,ld_struct$decay_summary$Chr)
cs <- setNames(ld_struct$decay_summary$c, ld_struct$decay_summary$Chr)
C <- cs[map$Chr]
a <- as[map$Chr]
b <- bs[map$Chr]

map[,rho_ld := 1 - (max_LD_with_QTN - b) / (C - b)]
map[,rho_d := a*bp_to_focal_QTN/(a*bp_to_focal_QTN+1)]
abline(0,1)

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


