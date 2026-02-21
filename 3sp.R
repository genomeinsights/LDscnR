
load("../LD-scaling-genome-scans/empirical_data/3sp/3sp_data.RData") ## contains SNP_res_3sp, GTs_3sp and pheno_3sp

##filter by maf
keep <- map_3sp$maf>0.1
GTs_3sp <- GTs_3sp[,map_3sp$maf>0.1]
map_3sp <- map_3sp[maf>0.1]
#if(any(ls() %in% "gds_3sp")) snpgdsClose(gds_3sp); unlink(gds_3sp)

gds_3sp <- create_gds_from_geno(geno = GTs_3sp, map=map_3sp,"gds_3sp.gds")

#------------------------------------------------------------
t1 <- Sys.time()
max_rho = 0.999
ld_struct_3sp <-  compute_ld_structure(gds_3sp,
                                       use         = "robust",
                                       max_rho     = max_rho,
                                       cores       = 8
)
t2 <- Sys.time()
difftime(t2,t1)
plot(ld_struct_3sp)


map_3sp[,emx_F:=readRDS("../LD-scaling-genome-scans/empirical_data/3sp/emx_3sp.rds")$F] ## add to map
emx_gif = map_3sp[,median(emx_F)/qf(0.5,1,115,lower.tail = FALSE)] ## inflation factor
map_3sp[,emx_F_GC:=emx_F/emx_gif]  ## genomic control
map_3sp[,emx_p_GC:=pf(emx_F_GC,1,115,lower.tail = FALSE)] ## p-value
map_3sp[,emx_q:=p.adjust(emx_p_GC,"fdr")] ## fdr correction

map_3sp[,lfmm_F:=readRDS("../LD-scaling-genome-scans/empirical_data/3sp/lfmm_F.rds")[keep]]
map_3sp[,lfmm_P:=pf(lfmm_F,1,115,lower.tail = FALSE)]
map_3sp[,lfmm_q:=p.adjust(lfmm_P,"fdr")]

SNP_ids <- map_3sp$marker
n_inds <- nrow(GTs_3sp)
#SNP_ids == map_3sp$marker
table(ld_struct_3sp$by_chr$Chr1$edges$SNP1 %in% map_3sp$marker)
#saveRDS(draws_joint,"draws_joint.rds")

t1 <- Sys.time()
max_rho = 0.999
ld_struct_3sp <-  compute_ld_structure(gds_3sp,
                                       use         = "robust",
                                       max_rho     = max_rho,
                                       cores       = 8
)


draws_3sp <- ld_rho_draws(gds = gds_3sp,
                            ld_struct  = ld_struct_3sp,
                            F_vals     = map_3sp[,.(emx_F_GC,lfmm_F)],
                            q_vals     = map_3sp[,.(emx_q,lfmm_q)],
                            n_inds     = nrow(GTs_3sp),
                            n_rho      = 40,
                            n_or_draws = 25,
                            rho_w_lim  = list(min=0.75,max=max_rho),
                            rho_d_lim  = list(min=0.9,max=max_rho),
                            rho_ld_lim = list(min=0.5,max=max_rho),
                            alpha_lim  = list(min=1.31,max=4),
                            lmin_lim   = list(min=1,max=30),
                            cores      = 8,
                            use        = "robust",
                            mode       = "joint"

)
t2 <- Sys.time()
difftime(t2,t1)
# draws_3sp_per_meth <- ld_rho_draws(
#   gds        = gds_3sp,
#   ld_struct  = ld_struct_3sp,
#   F_vals     = map_3sp[,.(emx_F_GC,lfmm_F)],
#   q_vals     = map_3sp[,.(emx_q,lfmm_q)],
#   n_inds     = n_inds,
#   SNP_ids    = SNP_ids,
#   n_rho      = 100,
#   n_or_draws = 25,
#   rho_w_lim  = list(min=0.75,max=max_rho),
#   rho_d_lim  = list(min=0.9,max=max_rho),
#   rho_ld_lim = list(min=0.5,max=max_rho),
#   alpha_lim  = list(min=1.31,max=4),
#   lmin_lim   = list(min=1,max=30),
#   cores      = 8,
#   mode       = "per_method"
# )

#draws_3sp
saveRDS(draws_3sp,"draws_3sp.rds")
#draws_wide$draws[OR_size>0,hist(l)]
#tmp <- consistency_score(draws_3s)
table(unique(unlist(draws_3sp$draws$OR)) %in% map_3sp$marker)

# table(tmp$consistency[,SNP] %in% map_3sp$marker)
# tmp$consistency$C[match(map_3sp$marker,tmp$consistency[,SNP])])


map2 <- add_consistency_to_map(map_3sp, consistency_obj = consistency_score(draws_3sp))
map2[,ld_w:=compute_ld_w(ld_struct_3sp,0.95)]

map2[,plot(lfmm_F,ld_w)]
#map2[,plot(ld_w)]

#map2[!is.na(Joint_C),plot(Joint_C)]
sign_th = 0.2
map2[,]
idx <- which(map2[,Joint_C>sign_th])
el  <- get_el(gds_3sp, idx, slide_win_ld = -1,by_chr = TRUE)


ORs_tbl <- detect_or(el,
                     q_vals=map2[,.(Joint_C)],
                     SNP_ids = map_3sp$marker,
                     SNP_chr = map_3sp$Chr,
                     ld_struct=ld_struct_3sp,
                     sign_th=sign_th,
                     sign_if = "greater",
                     rho_d=0.95,
                     rho_ld=0.95,
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

map2[, OR_factor := factor(OR_id, levels = sample(unique_ORs))]

#map2[!is.na(OR_factor),.(max(max_LD_with_QTN),.N,focal_QTN[which.max(max_LD_with_QTN)]),by=OR_factor]

tmp <- map2[, .(
  bp = Pos,
  Chr,
  marker,
  lfmm_q=-log10(lfmm_q),
  Joint_C,
  ld_w,
  OR_factor
)]
#tmp[Joint_C>0.05, Joint_C:=0.05]
layout <- prep_manhattan(
  tmp
)


plot_manhattan_gg(
  layout,
  y_vars = c("lfmm_q","Joint_C","ld_w"),
  y_labels = c("lfmm_q","Joint_C","ld_w"),
  thresholds = c(1.31,sign_th, NA),
  or_var = "OR_factor",
  ncol=1
)
#
# ## Analysis pipeline for three spined sticklebacks.
# if(!file.exists("./empirical_data/3sp/SNP_res_3sp.rds")){
#   #### ==== Estimate LD-decay for chromosomes ==== ####
#   if(!file.exists("./empirical_data/3sp/LD_decay_3sp.rds")){
#
#     b <- get_bg_ld(gds_3sp, n_sub = 5000, q = 0.95)
#
#     LD_decay_3sp <- ld_decay(gds=gds_3sp,
#                              q = 0.95,
#                              b = b,
#                              n_sub = 5000,
#                              slide_win_ld = 1000,
#                              window_size = 1e6,
#                              step_size = 5e5,
#                              dist_unit = 5000,
#                              n_cores_ld =  8)
#
#     LD_decay_3sp$summary[,plot(Chr_size,1/a)]
#     saveRDS(LD_decay_3sp,"./empirical_data/3sp/LD_decay_3sp.rds")
#
#   }else{
#     LD_decay_3sp <- readRDS("./empirical_data/3sp/LD_decay_3sp.rds")
#   }
#
#
#   #### ==== Get LD-pruned data using single linkage clustering with threshold based on 75% quantile of LD-decay rate ==== ####
#   if(!file.exists("./empirical_data/3sp/pruning.rds")){
#
#
#     cls_05_3sp <- SLC(gds_3sp,decay_tbl = LD_decay_3sp,idx = NULL,slide_win_ld = 1000,q = 0.5)$CLS
#
#     length(cls_05_3sp)/length(unlist(cls_05_3sp)) ## data compression
#
#     pruned_SNPs <- sapply(cls_05_3sp,sample,1) ## from clusters, pick one randomly
#
#     saveRDS(list(cls=cls_05_3sp,pruned_SNPs=pruned_SNPs),"./empirical_data/3sp/pruning.rds")
#
#   }else{
#     pruned_SNPs <- readRDS("./empirical_data/3sp/pruning.rds")$pruned_SNPs
#   }
#
#   #### ==== EMMAX analyses based on LD-pruned data ==== ####
#
#
#   if(!file.exists("./empirical_data/3sp/emx_3sp.rds")){
#     ## GRM based on pruned SNPs
#     GRM <- snpgdsGRM(gds_3sp,method = "GCTA",snp.id = pruned_SNPs,verbose = FALSE,autosome.only = FALSE)$grm
#
#     ## the binary phenotype
#     eco_bin  <- as.numeric(as.factor(pheno_3sp$ecotype))
#
#     ## EMMAX does not expect a file in 012 format so the maximum likelihood genotypes can be used (without rounding)
#     emx <- emmax(eco_bin,GTs_3sp,K = GRM) ## function lives in ./R/emmax.R
#
#     map_3sp[,emx_F:=emx$F] ## add to map
#     emx_gif = map_3sp[,median(emx_F)/qf(0.5,1,115,lower.tail = FALSE)] ## inflation factor
#     map_3sp[,emx_F_GC:=emx_F/1]  ## genomic control
#     map_3sp[,emx_p_GC:=pf(emx_F_GC,1,115,lower.tail = FALSE)] ## p-value
#     map_3sp[,emx_q:=p.adjust(emx_p_GC,"fdr")] ## fdr correction
#
#     saveRDS(emx,"./empirical_data/3sp/emx_3sp.rds")
#   }else{
#     map_3sp[,emx_F:=readRDS("./empirical_data/3sp/emx_3sp.rds")$F] ## add to map
#     emx_gif = map_3sp[,median(emx_F)/qf(0.5,1,115,lower.tail = FALSE)] ## inflation factor
#     map_3sp[,emx_F_GC:=emx_F/emx_gif]  ## genomic control
#     map_3sp[,emx_p_GC:=pf(emx_F_GC,1,115,lower.tail = FALSE)] ## p-value
#     map_3sp[,emx_q:=p.adjust(emx_p_GC,"fdr")] ## fdr correction
#   }
#
#   #### ==== Latent factor mixed model analyses (LFMM) ==== ####
#   #lfmm takes very long time
#   if(!file.exists("./empirical_data/3sp/lfmm_F.rds")){
#     library(LEA)
#     phe <- as.numeric(as.factor(pheno_3sp$ecotype))
#
#     write.lfmm(GTs_3sp, "./tmp/genotypes.lfmm")
#     write.env(phe, "./tmp/gradients.env")
#     project = NULL
#
#     project = lfmm2("./tmp/genotypes.lfmm", "./tmp/gradients.env", K=7) # K is user defined
#     pv = lfmm2.test(project, "./tmp/genotypes.lfmm", "./tmp/gradients.env",genomic.control = TRUE,full = TRUE)
#
#     saveRDS(pv$f,"./empirical_data/3sp/lfmm_F.rds") # only F-value is needed
#
#   }else{
#     map_3sp[,lfmm_F:=readRDS("./empirical_data/3sp/lfmm_F.rds")[keep]]
#     map_3sp[,lfmm_P:=pf(lfmm_F,1,115,lower.tail = FALSE)]
#     map_3sp[,lfmm_q:=p.adjust(lfmm_P,"fdr")]
#   }
#
#
#   #### ==== Draw 100 rho_w values and get ld_w vector for each of them ==== ####
#
#   if(!file.exists("./empirical_data/3sp/ld_w_draws_3sp_200.rds")){
#     ld_w_draws_3sp <- get_ld_w_draws(gds_3sp,decay_tbl = LD_decay_3sp,slide_win_ld = 1000,n_draws = 200,rho_min =  0.75,rho_max =  1)
#     saveRDS(ld_w_draws_3sp,"./empirical_data/3sp/ld_w_draws_3sp_200.rds")
#   }else{
#     ld_w_draws_3sp <- readRDS("./empirical_data/3sp/ld_w_draws_3sp_200.rds")
#   }
#
#   #### ==== For each of the 100 draws, draw 25 values for alpha, rho_OR and l_min and get C-scores ==== ####
#   if(!file.exists("./empirical_data/3sp/OR_draws_3sp_lmin1_30.rds")){
#
#     df2 = nrow(GTs_3sp)-2
#     F_vals <- map_3sp[,.(emx=emx_F_GC,lfmm=lfmm_F)]
#
#     q_orgs <- apply(F_vals,2,function(Fval){
#       p <- pf(Fval, df1=1, df2, lower.tail = FALSE)
#       q <- p.adjust(p,"fdr")
#     })
#
#
#     OR_draws_3sp <- get_ORs_for_draw(gds=gds_3sp,
#                                      F_vals = F_vals,
#                                      q_orgs = q_orgs,
#                                      decay_tbl = LD_decay_3sp,
#                                      ld_w_draws = ld_w_draws_3sp,
#                                      n_inds=nrow(GTs_3sp),
#                                      n_draws = 25,
#                                      rho_OR_lim = NULL,
#                                      rho_ld_lim = list(min=0.5, max=1.0),
#                                      rho_d_lim = list(min=0.90, max=1.0),
#                                      alpha_lim = list(min=-log10(0.05), max=4),
#                                      lmin_lim = list(min=1, max=30),
#                                      n_cores = 1)
#
#
#     saveRDS(OR_draws_3sp,"./empirical_data/3sp/OR_draws_3sp_lmin1_30.rds")
#
#   }
#
#   saveRDS(map_3sp,"./empirical_data/3sp/SNP_res_3sp.rds")
#
# }
#
# ## --------------------------
# ## Manhattan for 3sp stickleback data
# ## --------------------------
#
# ## data created above and/or available from Zenondo
# OR_draws_3sp <- readRDS("./empirical_data/3sp/OR_draws_3sp_lmin1_30.rds")
# SNP_res_3sp <- readRDS("./empirical_data/3sp/SNP_res_3sp.rds")
# SNP_res_3sp <- cbind(SNP_res_3sp,get_C(OR_draws_3sp,markers=SNP_res_3sp$marker))
# LD_decay_3sp <- readRDS("./empirical_data/3sp/LD_decay_3sp.rds")
#
# ## specify tau_C
# tau_C <- 0.2
#
# ## get ORs with rho_ld=0.999 and rho_d=0.999
# ORs <- find_ORs(gds_3sp,LD_decay_3sp,outliers = SNP_res_3sp[C_joint>tau_C,marker],rho_ld=0.999,rho_d=0.999)
#
# ## add ORs to data
# SNP_res_3sp[,OR:=ORs$OR[match(marker, ORs$marker)]]
# SNP_res_3sp[,n_loci:=ORs$n_loci[match(marker, ORs$marker)]]
#
#
# plot_data_manh <- prep_manhattan(SNP_res_3sp[,.(bp=Pos,Chr,marker,C_joint, OR = as.character(OR),     # <- Joint_C ORs
#                                                 LFMM_log = -log10(lfmm_q),
#                                                 lfmm_q)],chr_cols = c("white","grey80"),spacer =0)
#
#
#
# ## manhattan plot for lfmm
# p1 <- ggplot(plot_data_manh$data, aes(BPcum, LFMM_log, col = OR)) +
#   geom_rect(
#     data = plot_data_manh$rect,
#     aes(xmin = x1, xmax = x2, ymin = y1, ymax = Inf),
#     fill = plot_data_manh$rect$col,
#     alpha = 0.5,
#     linewidth = 0.25,
#     inherit.aes = FALSE
#   ) +
#   geom_point(size = 1) +
#   geom_point(data = plot_data_manh$data[!is.na(OR)],size = 1) +
#   geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.5) +
#   scale_color_manual(values = rep(col_vector, 4), guide = "none") +
#   scale_x_continuous(
#     label = plot_data_manh$axis$Chr,
#     breaks = plot_data_manh$axis$center
#   ) +
#   theme_bw() +
#   theme(
#     aspect.ratio = 0.25,
#     axis.text.x = element_text(angle = 90),
#     panel.grid = element_blank()
#   ) +
#   xlab(NULL) +
#   ylab(expression(-log[10](q)))
#
# ## manhattan plot for Joint-C
# p2 <- ggplot(plot_data_manh$data, aes(BPcum, C_joint, col = OR)) +
#   geom_rect(
#     data = plot_data_manh$rect,
#     aes(xmin = x1, xmax = x2, ymin = y1, ymax = Inf),
#     fill = plot_data_manh$rect$col,
#     alpha = 0.5,
#     linewidth = 0.25,
#     inherit.aes = FALSE
#   ) +
#   geom_point(size = 1) +
#   geom_point(data = plot_data_manh$data[!is.na(OR)],size = 1) +
#   geom_hline(yintercept = tau_C, linetype = 2, linewidth = 0.5) +
#   scale_color_manual(values = rep(col_vector, 4), guide = "none") +
#   scale_x_continuous(
#     label = plot_data_manh$axis$Chr,
#     breaks = plot_data_manh$axis$center
#   ) +
#   theme_bw() +
#   theme(
#     aspect.ratio = 0.25,
#     axis.text.x = element_text(angle = 90),
#     panel.grid = element_blank()
#   ) +
#   xlab(NULL) +
#   ylab(expression(Joint[C]))
#
#
#
#
# p_3sp <- grid.arrange(p1+ggtitle(expression("e) Three-spined sticklebacks | LFMM" )),
#                       p2+ggtitle(expression("f) Three-spined sticklebacks | Joint" ))
#                       ,ncol=1)
# saveRDS(p_3sp,"./figures/p_3sp.rds")
#
# nrow(SNP_res_3sp) # markers with MAF>0.1
# SNP_res_3sp[!is.na(OR),length(unique(OR))] # number of ORs
# SNP_res_3sp[!is.na(OR),length(unique(OR)),by=Chr] # number of ORs by Chr
#
# all_chr <- SNP_res_3sp[, .(Chr = unique(Chr))]
#
# # count ORs per chromosome
# cnt <- SNP_res_3sp[!is.na(OR),
#                    .(n_OR = uniqueN(OR)),
#                    by = Chr]
#
# # merge and fill zeros
# res <- merge(all_chr, cnt, by = "Chr", all.x = TRUE)
# res[is.na(n_OR), n_OR := 0]
# res[,mean(n_OR)]
