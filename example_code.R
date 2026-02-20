
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
LDscn <- LDscn_pipeline(geno = sim_ex$GTs,
                        map = sim_ex$map,
                        F_cols=c("emx_F","lfmm_F"),
                        q_cols=c("emx_q","lfmm_q"),
                        cores=8,
                        max_rho = max_rho,
                        rho_w_lim = list(min=0.5,max=max_rho),
                        rho_d_lim=list(min=0.5,max=max_rho+1e5),
                        rho_ld_lim=list(min=0.5,max=max_rho+1e5)
                        )
t2 <- Sys.time()
difftime(t2,t1)
plot(LDscn$ld_struct)

max_rho = 0.9999
LDscn$draws <- ld_rho_draws(
  ld_struct  = LDscn$ld_struct,
  F_vals     = map[,..F_cols],
  q_vals     = map[,..q_cols],
  n_inds     = n_inds,
  SNP_ids    = SNP_ids,
  n_rho      = n_rho,
  n_or_draws = n_or_draws,
  rho_w_lim  = list(min=0.75,max=max_rho),
  rho_d_lim  = list(min=0.8,max=max_rho),
  rho_ld_lim = list(min=0.5,max=max_rho),
  alpha_lim  = alpha_lim,
  lmin_lim   = lmin_lim,
  cores      = cores
)
#LDscn$draws$draws[,plot(rho_w, OR_size)]
LDscn$draws$n_or_draws
LDscn$consistency <- consistency_score(LDscn$draws)

map2 <- add_consistency_to_map(map, consistency_obj = LDscn$consistency)
map2[,ld_w:=compute_ld_w(LDscn$ld_struct,max_rho)]

#map2[,plot(-log10(lfmm_q_C))]
#map2[Chr=="Chr17",plot(C_mean)]

ORs_tbl <- detect_or(gds,
                     q_vals=map2[,.(lfmm_F_prime_C)],
                     ld_struct=LDscn$ld_str,
                     sign_th=0.05,
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

#map2[OR_factor == "Chr12_OR10"]

data_manh <-  map2[, .(
  bp = Pos,
  Chr,
  marker,
  C_mean,
  lfmm_q = -log10(lfmm_q),
  ld_w,
  OR_factor,
  type
)]
layout <- prep_manhattan(
  map2[, .(
    bp = Pos,
    Chr,
    marker,
    C_mean,
    lfmm_q = -log10(lfmm_q),
    ld_w,
    OR_factor,
    type
  )],spacer = 0
)


plot_manhattan_gg(
  layout,
  y_vars = c("lfmm_q","C_mean","ld_w"),
  y_labels = c("-log10(q) | LFMM","Joint_C","ld_w"),
  thresholds = c(-log10(0.05),0.05, NA),
  or_var = "OR_factor",
  type_var = "type",
  ncol=1
)

#####
or_dt <- LDscn$draws
consistency_score(LDscn$draws)

#or_dt$draws[method=="lfmm_F_prime",plot(rho_d,OR_size)]

long_dt <- or_dt$draws[
  , .(SNP = unlist(OR)),
  by = .(method)
]

C_dt <- long_dt[, .N, by = .(SNP,method)]

total_draws <- or_dt$n_rho*or_dt$n_or_draws

C_dt[, C := N / total_draws]

C_dt[method=="lfmm_F_prime",plot(C)]

C_obj <- structure(
  list(
    consistency = C_dt,
    total_draws = total_draws,
    n_rho       = or_dt$n_rho,
    n_or        = or_dt$n_or_draws
  ),
  class = "ld_consistency"
)
