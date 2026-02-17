
data(sim_ex)

t1 <- Sys.time()
LDscn <- LDscn_pipeline(geno = sim_ex$GTs,
                        map = sim_ex$map,
                        F_cols=c("emx_F","lfmm_F","lm_F"),
                        q_cols=c("emx_q","lfmm_q","lm_q"),
                        cores=8)
t2 <- Sys.time()
difftime(t2,t1)
plot_ld_decay(LDscn,max_dist = 1e6)
plot(LDscn$decay)
plot(LDscn$ld_struct)

map2 <- add_consistency_to_map(sim_ex$map, consistency_obj = LDscn$consistency)
map2[,ld_w:=compute_ld_w(LDscn$ld_str, LDscn$decay, 0.9)]
map2[is.na(C_mean),C_mean:=0]
ORs_tbl <- detect_or(q_vals=map2[,.(C_mean)],
                     ld_struct=LDscn$ld_str,
                     decay_obj=LDscn$decay,
                     sign_th=0.05,
                     sign_if = "greater",
                     rho_d=0.9,
                     rho_ld=0.9,
                     l_min = 1,
                     ret_table = TRUE)

map2 <- merge(map2,
              ORs_tbl,
              by.x = "marker",
              by.y = "SNP",
              all.x = TRUE)

map2[is.na(OR_id),OR_id:="ns"]
map2[,table(OR_id,method)]

unique_ORs <- unique(na.omit(map2$OR_id))
map2[, OR_factor := factor(OR_id, levels = sample(unique_ORs))]


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
  )]
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

