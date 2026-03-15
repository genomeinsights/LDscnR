library(ggplot2)
library(data.table)
library(SNPRelate)


devtools::load_all()

all_files <- list.files("../LD-scaling-genome-scans/results_sim/",full.names = TRUE)
all_files <- all_files[!grepl("PR",all_files)]
files <- c(all_files[grep("c2_V0.5",all_files)],all_files[grep("c1.5_V1",all_files)],all_files[grep("c1_V2",all_files)])
file <- files[1]

ld_test_th = 0.75
d_test_th  = 0.95
p_Va_th    = 0.05
cores = 8

done <- list.files("./PR_out/")
todo <- do.call(rbind,strsplit(files,"/"))
todo <- todo[,ncol(todo)]

for(file in files[which(!todo %in% done)][-1]){
  tmp <- readRDS(file)
  map <- tmp$SNP_res[,.(V,c,rep,Chr,Pos,marker,type, Chr_type,max_LD_with_QTN, bp_to_focal_QTN, focal_QTN, p_Va,emx_F,lfmm_F,emx_q,lfmm_q)]
  geno <- tmp$GTs
  tmp <- strsplit(file,"/")[[1]]
  out_file <- tmp[length(tmp)]

  gds_path <- tempfile(fileext = ".gds")
  gds <- create_gds_from_geno(geno, map, gds_path)

  on.exit({ SNPRelate::snpgdsClose(gds); unlink(gds_path) }, add = TRUE)


  ld_decay <- compute_LD_decay(
    gds,
    ## for LD-decay and bg
    q = 0.95,
    ## for bg
    n_sub_bg = 5000,
    ## for decay
    n_win_decay = 5,
    max_pairs = 5000,
    max_SNPs_decay = Inf,
    n_strata = 10,
    overlap = 0.5,
    prob_robust = 0.95,
    keep_el = TRUE,
    slide=2000,
    cores = 10
  )


  as <- setNames(ld_decay$decay_sum$a_pred,ld_decay$decay_sum$Chr)
  bs <- setNames(ld_decay$decay_sum$b,ld_decay$decay_sum$Chr)
  cs <- setNames(ld_decay$decay_sum$c_pred, ld_decay$decay_sum$Chr)
  C <- cs[map$Chr]
  a <- as[map$Chr]
  b <- bs[map$Chr]

  map[,rho_ld := 1 - (max_LD_with_QTN - b) / (C - b)]
  map[,rho_d := a*bp_to_focal_QTN/(a*bp_to_focal_QTN+1)]
  map[type=="QTN", rho_ld:=1]

  map_filt <- map[rho_ld > ld_test_th & rho_d < d_test_th]

  #min_alpha==0.75
  rho_bin_int_data <- rbindlist(lapply(-log10(c(0.05,0.25,0.5,0.75)),function(min_alpha){

    message("working on min_alpha=", 1/10^min_alpha,"\n")
    message("Computing data for rho_w\n")
    draws_ld_w <- ld_rho_draws(gds,
                               ld_decay  = ld_decay,
                               F_vals     = map[,.(emx_F,lfmm_F)],
                               q_vals     = map[,.(emx_q,lfmm_q)],
                               n_rho_w    = 40,
                               n_draws    = 25,
                               ld_w_int   = NULL,
                               stat_type  = "q",
                               rho_w_lim  = list(min=0.75,max=1-1e-10),
                               rho_d_lim  = list(min=0.9,max=0.99),
                               rho_ld_lim = list(min=0.5,max=0.99),
                               alpha_lim  = list(min=min_alpha,max=4),
                               lmin_lim   = list(min=1,max=10),
                               cores      = 10,
                               mode       = c("joint")
    )

    PR_w <- get_PR(draws_ld_w$draws$OR,map_filt,cores=10)


    C_scores <- data.table(add_consistency_to_map(map, consistency_obj = consistency_score(draws_ld_w$draws[method=="Joint"]))[,.(marker=marker,Joint_C=Joint_C)])

    PR_w <- cbind(draws_ld_w$draws[,.(method, rho_w, rho_d, rho_ld, alpha, l_min,n_ORs=OR_size)],PR_w)

    #PR_w$PR_dat
    PR_w <- PR_w[,.(AUC_PR="PR*",
                    C_scores=list(C_scores[!is.na(Joint_C)]),
                    PR_dat=list(PR_w),
                    AUC=get_AUC_OR(PR_w)$AUC$AUC_norm)]
    #PR_w
    draws_ld_w_C <- ld_rho_draws(gds,
                                 ld_decay  = ld_decay,
                                 F_vals     = NULL,
                                 C_scores   = C_scores[,.(Joint_C)],
                                 q_vals     = NULL,
                                 n_rho_w    = NULL,
                                 n_draws    = 1000,
                                 ld_w_int   = NULL,
                                 stat_type  = "C",
                                 rho_w_lim  = NULL,
                                 rho_d_lim  = list(min=0.9,max=0.99),
                                 rho_ld_lim = list(min=0.5,max=0.99),
                                 alpha_lim  = NULL,
                                 lmin_lim   = NULL,
                                 C_lim      = list(min=0,max=0.5),
                                 cores      = 10,
                                 mode       = c("joint")
    )

    PR_C <- get_PR(draws_ld_w_C$draws$OR,map_filt,cores=10)


    PR_C <- cbind(draws_ld_w_C$draws[,.(method, rho_w, rho_d, rho_ld, alpha, l_min,n_ORs=OR_size)],PR_C)

    PR_C <- PR_C[,.(AUC_PR="PR-C",
                    C_scores=list(NULL),
                    PR_dat=list(PR_C),
                    AUC=get_AUC_OR(PR_C)$AUC$AUC_norm)]



    PR_ld_w <- rbind(PR_w,PR_C)

    PR_ld_w[,ld:="rho_w"]

    bins_grid = c(5,10,20,40,80,160,320,640,1280)
    rhos = c(0.5,0.6,0.7,0.8,0.90,0.99)
    # rho = 0.9
    # bins = 20
    message("Computing data for DPI \n")
    PR_ld_int <- rbindlist(lapply(rhos,function(rho){
      rbindlist(lapply(bins_grid,function(bins){

        message("Working on rho=", rho, " and n_bins=",bins)

        ldint <- unlist(parallel_apply(ld_decay$by_chr, function(chr_obj) {

          a <- chr_obj$decay_sum$a

          compute_ld_int_fill(chr_obj$el, snp_ids=chr_obj$snp_ids, a,n_bins = bins,rho_max = rho,edge_mode = "fill")

        }, cores = cores))


        map[,ld_int := ldint]

        draws_ld_int <- ld_rho_draws(gds,
                                     ld_decay  = ld_decay,
                                     F_vals     = map[,.(emx_F,lfmm_F)],
                                     q_vals     = map[,.(emx_q,lfmm_q)],
                                     C_scores   = NULL,
                                     n_rho_w    = 1,
                                     n_draws    = 1000,
                                     ld_w_int   = map$ld_int,
                                     stat_type  = "q",
                                     rho_w_lim  = NULL,
                                     rho_d_lim  = list(min=0.9,max=0.99),
                                     rho_ld_lim = list(min=0.5,max=0.99),
                                     alpha_lim  = list(min=min_alpha,max=4),
                                     lmin_lim   = list(min=1,max=10),
                                     C_lim      = NULL,
                                     cores      = 10,
                                     mode       = c("joint")
        )

        PR_int <- get_PR(draws_ld_int$draws$OR,map_filt,cores=cores)

        #map$Joint_C <- NULL
        C_scores <- data.table(add_consistency_to_map(map, consistency_obj = consistency_score(draws_ld_int$draws[method=="Joint"]))[,.(marker=marker,Joint_C=Joint_C)])

        PR_int <- cbind(draws_ld_int$draws[,.(method, rho_w, rho_d, rho_ld, alpha, l_min,n_ORs=OR_size)],PR_int)

        PR_int <- PR_int[,.(AUC_PR="PR*",
                            rho = rho,
                            bins=bins,
                            C_scores=list(C_scores[!is.na(Joint_C)]),
                            ld_int=list(map$ld_int),
                            PR_dat=list(PR_int),
                            AUC=get_AUC_OR(PR_int)$AUC$AUC_norm)]

        #PR_int$PR_dat

        draws_ld_int_C <- ld_rho_draws(gds,
                                       ld_decay  = ld_decay,
                                       F_vals     = NULL,
                                       C_scores   = C_scores[,.(Joint_C)],
                                       q_vals     = NULL,
                                       n_rho_w    = NULL,
                                       n_draws    = 1000,
                                       ld_w_int   = NULL,
                                       stat_type  = "C",
                                       rho_w_lim  = NULL,
                                       rho_d_lim  = list(min=0.9,max=0.99),
                                       rho_ld_lim = list(min=0.5,max=0.99),
                                       alpha_lim  = NULL,
                                       lmin_lim   = NULL,
                                       C_lim      = list(min=0,max=0.5),
                                       cores      = 10,
                                       mode       = c("joint")
        )

        PR_C <- get_PR(draws_ld_int_C$draws$OR,map_filt,cores=cores)


        PR_C <- cbind(draws_ld_int_C$draws[,.(method, rho_w, rho_d, rho_ld, alpha, l_min,n_ORs=OR_size)], PR_C)

        PR_C <- PR_C[,.(AUC_PR="PR-C",
                        rho = rho,
                        bins=bins,
                        C_scores=list(NULL),
                        ld_int=list(NULL),
                        PR_dat=list(PR_C),
                        AUC=get_AUC_OR(PR_C)$AUC$AUC_norm)]


        #PR_C
        return(rbind(PR_int,PR_C))

      }))
    }),fill=TRUE,use.names = TRUE)
    PR_ld_int[,ld:="DPI"]


    cbind(rbind(PR_ld_w,PR_ld_int,fill=TRUE,use.names=TRUE),alpha_min=min_alpha)

  }))

  saveRDS(cbind(map[1,.(c,V,rep)],rho_bin_int_data),paste0("./PR_out/",out_file))

}


q("no")
PR_ld_int_fill[,edge:="fill"]
PR_ld_int[,edge:="none"]
dt <- rbind(PR_ld_int_fill,PR_ld_int,fill=TRUE,use.names=TRUE)


ggplot(dt[,], aes(bins,AUC,col=factor(rho))) +
  geom_line() +
  facet_grid(AUC_PR~edge) +
  geom_hline(yintercept = PR_ld_w$AUC)

ggplot(dt, aes(edge,AUC)) +
  geom_boxplot() +
  facet_grid(AUC_PR~.) +
  geom_hline(yintercept = PR_ld_w$AUC)

PR_ld_int[AUC_PR=="PR-C",plot(rho,AUC,type="l")]

t2 <- Sys.time()
difftime(t2,t1)
file <- files[8]

done <- list.files("./PR_out_robust_fill/")
rho_bin_int_data <- rbindlist(lapply(files[which(todo %in% done)],function(file){
  tmp <- strsplit(file,"/")[[1]]
  out_file <- tmp[length(tmp)]
  rho_bin_int_data <- readRDS(paste0("./PR_out_robust_fill/",out_file))
}))


dt <- rho_bin_int_data[,
  .(
    AUC_mean = mean(AUC, na.rm = TRUE),
    AUC_sd   = sd(AUC, na.rm = TRUE),
    n        = .N,
    AUC_se   = sd(AUC, na.rm = TRUE) / sqrt(.N),
    AUC_lwr  = mean(AUC, na.rm = TRUE) - 1.96 * sd(AUC, na.rm = TRUE) / sqrt(.N),
    AUC_upr  = mean(AUC, na.rm = TRUE) + 1.96 * sd(AUC, na.rm = TRUE) / sqrt(.N)
  ),
  by = .(AUC_PR, alpha_min, c, V, ld, rho, bins)
]

dt[, alpha := paste0("'alpha=", 1 / 10^(alpha_min),"'")]
dt[, rho := paste0("'rho=", rho,"'")]

rho_levels <- unique(dt[AUC_PR == "PR-C" & ld == "DPI", rho])

baseline_rhow <- dt[
  AUC_PR == "PR-C" & ld == "rho_w",
  .(AUC_mean = unique(AUC_mean)),
  by = .(alpha,c)
]

baseline_rhow <- baseline_rhow[
  , .(rho = rho_levels), by = .(alpha, c,AUC_mean)
]

p1 <- ggplot(
  dt[AUC_PR == "PR-C" & ld == "DPI" & c=="c2"],
  aes(log(bins), AUC_mean)
) +
  geom_line() +
  geom_hline(
   yintercept = baseline_rhow[1,AUC_mean], col="salmon"
  )+
  geom_point(size = 2, shape = 1) +
  #geom_errorbar(aes(ymin = AUC_lwr, ymax = AUC_upr), width = 0.02) +
  facet_grid(rho ~ alpha,labeller = label_parsed) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab("log(bins)") +
  ylab("AUC") +
  geom_ribbon(aes(ymin = AUC_lwr, ymax = AUC_upr, group = 1), alpha = 0.15, fill = "grey20") +
  geom_line() +
  #geom_point(size = 1.5, shape = 1)+
  geom_hline(
    data = baseline_rhow[c=="c2"],
    aes(yintercept = AUC_mean),
    linetype = 2
  ) +
  theme(
    aspect.ratio = 1,
    #strip.text = element_text(margin = margin(1, 1, 1, 1)),
    #plot.margin = margin(1, 1, 1, 1),
    #panel.spacing = unit(0.05, "lines"),
    #axis.text.x = element_text(angle = 90, hjust = 1),
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    #legend.position = "right"
  )

p1

p2 <- ggplot(
  dt[AUC_PR == "PR-C" & ld == "DPI" & c=="c1.5"],
  aes(log(bins), AUC_mean)
) +
  geom_line() +
  geom_hline(
    yintercept = baseline_rhow[1,AUC_mean], col="salmon"
  )+
  geom_point(size = 2, shape = 1) +
  #geom_errorbar(aes(ymin = AUC_lwr, ymax = AUC_upr), width = 0.02) +
  facet_grid(rho ~ alpha,labeller = label_parsed) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab("log(bins)") +
  ylab("AUC") +
  geom_ribbon(aes(ymin = AUC_lwr, ymax = AUC_upr, group = 1), alpha = 0.15, fill = "grey20") +
  geom_line() +
  #geom_point(size = 1.5, shape = 1)+
  geom_hline(
    data = baseline_rhow[c=="c1.5"],
    aes(yintercept = AUC_mean),
    linetype = 2
  ) +
  theme(
    aspect.ratio = 1,
    #strip.text = element_text(margin = margin(1, 1, 1, 1)),
    #plot.margin = margin(1, 1, 1, 1),
    #panel.spacing = unit(0.05, "lines"),
    #axis.text.x = element_text(angle = 90, hjust = 1),
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    #legend.position = "right"
  )

p2


dt <- rho_bin_int_data[bins>500  | is.na(bins),.(AUC=median(AUC)),by=.(AUC_PR,alpha_min,c,V,ld,rho)]

p1 <- ggplot(dt[AUC_PR=="PR-C"& ld=="DPI"],aes(rho,AUC,col=factor(1/10^(alpha_min)))) +
  geom_line() +
  geom_point(size=2,shape=1) +
  facet_grid(c~.,scales="free_y") +
  scale_color_viridis_d(option="turbo",name="alpha")+
  theme_bw() +
  theme(aspect.ratio = 1)+
  xlab("Max rho")  +
  geom_hline(data =dt[AUC_PR=="PR-C" & ld=="rho_w" ],linetype=2,aes(yintercept = AUC,col=factor(1/10^(alpha_min))))



dt <- rho_bin_int_data[rho>0.9 | is.na(rho),.(AUC=median(AUC)),by=.(AUC_PR,alpha_min,c,V,ld,bins)]

p2 <- ggplot(dt[AUC_PR=="PR-C" & ld=="DPI"],aes(log(bins),AUC,col=factor(1/10^(alpha_min)))) +
  geom_line() +
  geom_point(size=2,shape=1) +
  facet_grid(c~.,scales="free_y") +
  theme_bw() +
  scale_color_viridis_d(option="turbo",name="alpha")+
  theme(aspect.ratio = 1) +
  xlab("No. bins")+
  geom_hline(data =dt[AUC_PR=="PR-C" & ld=="rho_w" ],linetype=2,aes(yintercept = AUC,col=factor(1/10^(alpha_min))))


p1 | p2
#+
  #facet_grid(alpha_min~.)+
  geom_hline(data=dt[is.na(n_bins)& AUC_PR=="PR-C" ],aes(yintercept = AUC, col=factor(alpha_min),linetype=AUC_PR))


ggplot(rho_bin_int_data[],aes(factor(signif(1/10^(alpha_min),1)),AUC,fill=factor(ld))) +
  geom_boxplot() +
  facet_grid(AUC_PR~c+V) +
  theme_bw() +
  xlab("Minimum alpha")

ggplot(rho_bin_int_data[ld=="int"],aes(factor(n_bins),AUC,fill=factor(alpha_min))) +
  geom_boxplot() +
  facet_wrap(c~V+alpha_min)

ggplot(rho_bin_int_data[AUC_PR=="PR-C"],aes(n_bins,AUC,col=factor(-log10(alpha_min)))) +
  geom_line() +
  #facet_grid(alpha_min~.)+
  geom_hline(data=rho_bin_int_data[ld=="rho_w" & AUC_PR=="PR-C"],aes(yintercept = AUC, col=factor(alpha_min)))


rho_bin_int_data[ld=="rho_w" & AUC_PR=="PR-C",plot(alpha_min,AUC,type="l")]
rho_bin_int_data[is.na(n_bins),]


rho_bin_int_data[,cor.test(R,AUC)]

rho_bin_int_data[is.na(n_bins)][which.max(R),PR_dat[[1]]][,boxplot(PR)]

map2 <- add_consistency_to_map(map, consistency_obj = rho_bin_int_data[!is.na(n_bins)][which.max(R),C_scores][[1]])

p1 <- ggplot(map2, aes(1:nrow(map),Joint_C,col=max_LD_with_QTN,shape=type,size=type)) +
  geom_point() +
  scale_color_viridis_c(option="turbo") +
  scale_shape_manual(values=c(20,3)) +
  scale_size_manual(values=c(1,3)) +
  geom_hline(yintercept = 0.03)

rho_bin_int_data[is.na(n_bins)][which.max(R)]

map2 <- add_consistency_to_map(map, consistency_obj = rho_bin_int_data[is.na(n_bins)][which.max(R),C_scores][[1]])

p2 <- ggplot(map2, aes(1:nrow(map),Joint_C,col=max_LD_with_QTN,shape=type,size=type)) +
  geom_point() +
  scale_color_viridis_c(option="turbo") +
  scale_shape_manual(values=c(20,3)) +
  scale_size_manual(values=c(1,3)) +
  geom_hline(yintercept = 0.03)

p1 / p2



plot(sapply(rho_bin_int_data[!is.na(n_bins),ld_int],function(ldint) cor(map$max_LD_with_QTN,ldint)^2),rho_bin_int_data[!is.na(n_bins),AUC],type="l")


plot(rho_bin_int_data[!is.na(n_bins),n_bins],rho_bin_int_data[!is.na(n_bins),AUC]/1.5,type="l")
lines(rho_bin_int_data[!is.na(n_bins),n_bins],sapply(rho_bin_int_data[!is.na(n_bins),ld_int],function(ldint) cor(map$max_LD_with_QTN,ldint)^2),type = "l")

#1sapply(rho_bin_int_data[is.na(n_bins),ld_int],function(ldint) cor(map$max_LD_with_QTN,ldint)^2)


rho_bin_int_data[!is.na(n_bins),]

plot(sapply(rho_bin_int_data[!is.na(n_bins),ld_int],function(ldint) cor(map$max_LD_with_QTN,ldint)^2),rho_bin_int_data[!is.na(n_bins),R])
#saveRDS(rho_bin_int_data,"rho_bin_int_data.rds")


#save.image()
# PR_int[,ld:="int"]
# PR_int[,method := gsub("_int","",method)]
# PR_ld_w[,ld:="rho_w"]
# PR_data <- rbind(PR_ld_w,PR_int)



ggplot(PR_data,
       aes(x = recall,
           y = precision,
           colour = ld)) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE) +
  theme_bw() +
  labs(title = "Precision–Recall tradeoff across parameter draws")

ggplot(PR_data,
       aes(x = ld,
           y = PR,
           fill = ld)) +
  geom_boxplot(width = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distribution of PR score across parameter draws")

PR_data[method=="Joint" & ld=="int"][which.max(PR)]
PR_data[method=="Joint" & ld=="rho_w"][which.max(PR)]

# PR_data[method=="Joint" & ld=="int"][which.max(PR)]
# PR_data[method=="Joint" & ld=="rho_w"][which.max(PR)]

map2 <- add_consistency_to_map(map, consistency_obj = consistency_score(draws_ld_int$draws[method=="Joint"]))


