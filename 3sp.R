
load("../LD-scaling-genome-scans/empirical_data/3sp/3sp_data.RData") ## contains SNP_res_3sp, GTs_3sp and pheno_3sp

##filter by maf
keep <- map_3sp$maf>0.1
GTs_3sp <- GTs_3sp[,map_3sp$maf>0.1]
map_3sp <- map_3sp[maf>0.1]
#if(any(ls() %in% "gds_3sp")) snpgdsClose(gds_3sp); unlink(gds_3sp)

gds_3sp <- create_gds_from_geno(geno = GTs_3sp, map=map_3sp,"gds_3sp.gds")

#------------------------------------------------------------
t1 <- Sys.time()
max_rho = 0.9
ld_struct_3sp_w1000 <-  compute_ld_structure(gds_3sp,
                                       use         = "robust",
                                       max_rho     = max_rho,
                                       adapt_thin = FALSE,
                                       cores       = 8
)
t2 <- Sys.time()
difftime(t2,t1)
plot(ld_struct_3sp_w1000)
#rm(ld_struct_3sp_w1000)
gc()

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
table(ld_struct_3sp_w1000$by_chr$Chr1$edges$SNP1 %in% map_3sp$marker)
#saveRDS(draws_joint,"draws_joint.rds")

# max_rho = 0.99
# t1 <- Sys.time()
# ld_struct_3sp_w1000 <-  compute_ld_structure(gds_3sp,
#                                          use         = "robust",
#                                          max_rho     = max_rho,
#                                          adapt_thin  = FALSE
# )
# t2 <- Sys.time()
# difftime(t2,t1)
# gc()
# t1 <- Sys.time()
# max_rho = 0.999
# ld_struct_3sp <-  compute_ld_structure(gds_3sp,
#                                        use         = "robust",
#                                        max_rho     = max_rho,
#                                        cores       = 8
# )

t1 <- Sys.time()
draws_3sp_000 <- ld_rho_draws(gds = gds_3sp,
                            ld_struct  = ld_struct_3sp_w1000,
                            F_vals     = map_3sp[,.(emx_F_GC,lfmm_F)],
                            q_vals     = map_3sp[,.(emx_q,lfmm_q)],
                            n_inds     = nrow(GTs_3sp),
                            n_rho      = 100,
                            n_or_draws = 25,
                            rho_w_lim  = list(min=0,max=max_rho),
                            rho_d_lim  = list(min=0,max=max_rho),
                            rho_ld_lim = list(min=0,max=max_rho),
                            alpha_lim  = list(min=1.31,max=4),
                            lmin_lim   = list(min=1,max=30),
                            cores      = 1,
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
#saveRDS(draws_3sp,"draws_3sp.rds")
#draws_wide$draws[OR_size>0,hist(l)]
#tmp <- consistency_score(draws_3s)
table(unique(unlist(draws_3sp$draws$OR)) %in% map_3sp$marker)

# table(tmp$consistency[,SNP] %in% map_3sp$marker)
# tmp$consistency$C[match(map_3sp$marker,tmp$consistency[,SNP])])


map2 <- add_consistency_to_map(map_3sp, consistency_obj = consistency_score(draws_3sp_03_05_05))
map2[,ld_w:=compute_ld_w(ld_struct_3sp_w1000,0.3)]

#map2[,plot(lfmm_F,ld_w)]
#map2[,plot(ld_w)]

#map2[!is.na(Joint_C),plot(Joint_C)]
sign_th = 0.2
#map2[,]
idx <- which(map2[,Joint_C>sign_th])
el  <- get_el(gds_3sp, idx, slide_win_ld = -1,by_chr = TRUE)


ORs_tbl <- detect_or(el,
                     q_vals=map2[,.(Joint_C)],
                     SNP_ids = map_3sp$marker,
                     SNP_chr = map_3sp$Chr,
                     ld_struct=ld_struct_3sp_w1000,
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

##########
decay_sum <- ld_struct_3sp_w1000$by_chr$Chr1$summary


plot(seq(0.85,0.95,length.out=11),d_from_rho(a=decay_sum$robust$a,  d0 = decay_sum$robust$d0, rho = seq(0.85,0.95,length.out=11)))
par(mfcol=c(3,1))
plot(compute_ld_w(ld_struct_3sp_w1000,0.85),main="0.85")
plot(compute_ld_w(ld_struct_3sp_w1000,0.90),main="0.9")
plot(compute_ld_w(ld_struct_3sp_w1000,0.95),main="0.95")
gc()

