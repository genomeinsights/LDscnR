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

for(file in files[which(!todo %in% done)]){
  tmp <- readRDS(file)
  map <- tmp$SNP_res[,.(V,c,rep,Chr,Pos,marker,type, Chr_type,max_LD_with_QTN, bp_to_focal_QTN, focal_QTN, p_Va,emx_F,lfmm_F,emx_q,lfmm_q)]
  geno <- tmp$GTs
  tmp <- strsplit(file,"/")[[1]]
  out_file <- tmp[length(tmp)]
  message("Creating gds file")

  gds_path <- tempfile(fileext = ".gds")
  gds <- create_gds_from_geno(geno, map, gds_path)

  on.exit({ SNPRelate::snpgdsClose(gds); unlink(gds_path) }, add = TRUE)

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
    n_bins_ld_int = 20, ## does not matter, ld_int is recalculated below for different bin sizes
    rho = 1-1e-10,
    k_max=1000,
    cores = 8
  )


  ld_struct$decay_sum$chr_size <- map[,max(Pos),by=Chr][,V1]

  d_mod <- rlm(log(a) ~ log(chr_size), data = ld_struct$decay_sum)

  ld_struct$decay_sum[,a_pred := exp(predict(d_mod,ld_struct$decay_sum))]
  ld_struct$decay_sum[,plot(a_pred,a)]
  abline(0,1)

  as <- setNames(ld_struct$decay_sum$a_pred,ld_struct$decay_sum$Chr)
  bs <- setNames(ld_struct$decay_sum$b,ld_struct$decay_sum$Chr)
  cs <- setNames(ld_struct$decay_sum$c, ld_struct$decay_sum$Chr)
  C <- cs[map$Chr]
  a <- as[map$Chr]
  b <- bs[map$Chr]

  map[,rho_ld := 1 - (max_LD_with_QTN - b) / (C - b)]
  map[,rho_d := a*bp_to_focal_QTN/(a*bp_to_focal_QTN+1)]
  map[type=="QTN", rho_ld:=1]

  map_filt <- map[rho_ld > ld_test_th & rho_d < d_test_th]

  rho_bin_int_data <- rbindlist(lapply(c(1.31,1,0.5,0.1),function(min_alpha){

    draws_ld_w <- ld_rho_draws(gds,
                               ld_struct  = ld_struct,
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
                               cores      = 8,
                               mode       = c("joint")
    )

    PR_w <- get_PR(draws_ld_w$draws$OR,map_filt,cores=8)
    C_scores <- data.table(add_consistency_to_map(map, consistency_obj = consistency_score(draws_ld_w$draws[method=="Joint"]))[,.(marker=marker,Joint_C=Joint_C)])

    PR_w <- cbind(draws_ld_w$draws[,.(method, rho_w, rho_d, rho_ld, alpha, l_min,n_ORs=OR_size)],PR_w)


    PR_w <- PR_w[,.(AUC_PR="PR*",
                    C_scores=list(C_scores[!is.na(Joint_C)]),
                    PR_dat=list(PR_w),
                    AUC=get_AUC_OR(PR_w)$AUC$AUC_norm)]

    draws_ld_w_C <- ld_rho_draws(gds,
                                 ld_struct  = ld_struct,
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
                                 cores      = 8,
                                 mode       = c("joint")
    )

    PR_C <- get_PR(draws_ld_w_C$draws$OR,map_filt,cores=8)


    PR_C <- cbind(draws_ld_w_C$draws[,.(method, rho_w, rho_d, rho_ld, alpha, l_min,n_ORs=OR_size)],PR_C)

    PR_C <- PR_C[,.(AUC_PR="PR-C",
                    C_scores=list(NULL),
                    PR_dat=list(PR_C),
                    AUC=get_AUC_OR(PR_C)$AUC$AUC_norm)]



    PR_ld_w <- rbind(PR_w,PR_C)

    PR_ld_w[,ld:="rho_w"]


    bin_grid = seq(20,50,by=2)
    #bins =30
    PR_ld_int <- rbindlist(lapply(bin_grid,function(bins){

      cat(bins," -- ")
      ldint <- unlist(parallel_apply(ld_struct$by_chr, function(chr_obj) {

        a <- chr_obj$decay_sum$a
        b <- chr_obj$decay_sum$b

        #ld_int_vect <- compute_ld_int(chr_obj$el, snp_ids=chr_obj$snp_ids, a, n_bins=bins)

        # bin width based on LD decay scale
        dist_unit <- delta / a

        el <- chr_obj$el[d<max(chr_obj$el$d)/100]
        #assign linear bins
        el[, dist_idx := floor(d / dist_unit)]
        el[,length(unique(dist_idx))]

        #max_d <- max(chr_obj$el$d)/10
        #n_bins <- ceiling(chr_obj$decay_sum$a * max_d / 2)
        #

        #el <- chr_obj$el[d<max_d/10]
        ld_int_vect_adaptive <- compute_ld_int_adaptive(chr_obj$el[d<max(chr_obj$el$d)/100], snp_ids=chr_obj$snp_ids, a, delta = 5,min_pairs = 0)

      }, cores = cores))

      map[,ld_int := ldint]

      draws_ld_int <- ld_rho_draws(gds,
                                   ld_struct  = ld_struct,
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
                                   cores      = 8,
                                   mode       = c("joint")
      )

      PR_int <- get_PR(draws_ld_int$draws$OR,map_filt,cores=cores)
      C_scores <- data.table(add_consistency_to_map(map, consistency_obj = consistency_score(draws_ld_int$draws[method=="Joint"]))[,.(marker=marker,Joint_C=Joint_C)])

      PR_int <- cbind(draws_ld_int$draws[,.(method, rho_w, rho_d, rho_ld, alpha, l_min,n_ORs=OR_size)],PR_int)

      PR_int <- PR_int[,.(AUC_PR="PR*",
                          n_bins=bins,
                          C_scores=list(C_scores[!is.na(Joint_C)]),
                          ld_int=list(map$ld_int),
                          PR_dat=list(PR_int),
                          AUC=get_AUC_OR(PR_int)$AUC$AUC_norm)]

      PR_int
      draws_ld_int_C <- ld_rho_draws(gds,
                                     ld_struct  = ld_struct,
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
                                     cores      = 8,
                                     mode       = c("joint")
      )

      PR_C <- get_PR(draws_ld_int_C$draws$OR,map_filt,cores=cores)


      PR_C <- cbind(draws_ld_int$draws[,.(method, rho_w, rho_d, rho_ld, alpha, l_min,n_ORs=OR_size)], PR_C)

      PR_C <- PR_C[,.(AUC_PR="PR-C",
                      n_bins= bins,
                      C_scores=list(NULL),
                      ld_int=list(NULL),
                      PR_dat=list(PR_C),
                      AUC=get_AUC_OR(PR_C)$AUC$AUC_norm)]



      return(rbind(PR_int,PR_C))

    }))

    PR_ld_int[,ld:="int"]

    cbind(rbind(PR_ld_w,PR_ld_int,fill=TRUE,use.names=TRUE),alpha_min=min_alpha)

  }))

  saveRDS(cbind(map[1,.(c,V,rep)],rho_bin_int_data),paste0("./PR_out/",out_file))

}
#rm(bins)
q("no")
t2 <- Sys.time()
difftime(t2,t1)
file <- files[8]
rm(rbindlist)
rho_bin_int_data <- rbindlist(lapply(files[which(todo %in% done)],function(file){
  tmp <- strsplit(file,"/")[[1]]
  out_file <- tmp[length(tmp)]
  rho_bin_int_data <- readRDS(paste0("./PR_out/",out_file))
}))
rho_bin_int_data

dt <- rho_bin_int_data[,.(AUC=mean(AUC)),by=.(AUC_PR,alpha_min,c,V)]

ggplot(dt[!is.na(n_bins) & AUC_PR=="PR-C"],aes(n_bins,AUC,col=factor(alpha_min),linetype=AUC_PR)) +
  geom_line() +
  facet_wrap(c~V)+
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


