library(LDscnR)
data(sim_ex)

geno = sim_ex$GTs
map = sim_ex$map
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
max_rho = 0.99
ld_struct <-  compute_ld_structure(gds,
                                   use         = "robust",
                                   max_rho     = max_rho,
                                   K_target    = 100
                                   )
plot(ld_struct)


# draws_per_method <- ld_rho_draws(
#   ld_struct  = ld_struct,
#   F_vals     = map[,..F_cols],
#   q_vals     = map[,..q_cols],
#   n_inds     = n_inds,
#   SNP_ids    = SNP_ids,
#   n_rho      = n_rho,
#   n_or_draws = n_or_draws,
#   rho_w_lim  = list(min=0.9,max=max_rho),
#   rho_d_lim  = list(min=0.9,max=max_rho),
#   rho_ld_lim = list(min=0.5,max=max_rho),
#   alpha_lim  = list(min=1.31,max=3),
#   lmin_lim   = list(min=1,max=5),
#   cores      = cores,
#   mode       = "per_method"
#
# )

draws_joint <- ld_rho_draws(gds,
                            ld_struct  = ld_struct,
                            F_vals     = map[,..F_cols],
                            q_vals     = map[,..q_cols],
                            n_inds     = nrow(geno),
                            n_rho      = 40,
                            n_or_draws = 25,
                            rho_w_lim  = list(min=0.75,max=max_rho),
                            rho_d_lim  = list(min=0.9,max=max_rho),
                            rho_ld_lim = list(min=0.5,max=max_rho),
                            alpha_lim  = list(min=1.31,max=4),
                            lmin_lim   = list(min=1,max=10),
                            cores      = 8,
                            use        = "robust",
                            mode       = "joint"

)


#draws_joint$draws[,plot(OR_size)]
map2 <- add_consistency_to_map(map, consistency_obj = consistency_score(draws_joint))
#map2 <- add_consistency_to_map(map2, consistency_obj = consistency_score(draws_per_method))
map2[,ld_w:=compute_ld_w(ld_struct,max_rho)]


#map2[,cor.test(ld_w,lfmm_F)]
#map2[,plot(ld_w)]

#map2[!is.na(Joint_C),plot(Joint_C)]
sign_th = 0.05
idx <- which(map2[,Joint_C>sign_th])
el  <- get_el(gds, idx, slide_win_ld = -1,by_chr = TRUE)


ORs_tbl <- detect_or(el,
                     q_vals=map2[,.(Joint_C)],
                     SNP_ids = map$marker,
                     SNP_chr = map$Chr,
                     ld_struct=ld_struct,
                     sign_th=sign_th,
                     sign_if = "greater",
                     rho_d=0.999,
                     rho_ld=0.999,
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
  y_vars = c("lfmm_q","Joint_C","ld_w"),
  y_labels = c("-log10(q) | LFMM","Joint_C","ld_w"),
  thresholds = c(-log10(0.05),0.05, NA),
  or_var = "OR_factor",
  type_var = "type",
  ncol=1
)

