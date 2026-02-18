
data(sim_ex)


t1 <- Sys.time()
LDscn <- LDscn_pipeline(geno = sim_ex$GTs,
                        map = sim_ex$map,
                        F_cols=c("emx_F","lfmm_F"),
                        q_cols=c("emx_q","lfmm_q"),
                        cores=8,
                        slide_win_ld = 1000
                        )
t2 <- Sys.time()
difftime(t2,t1)

plot_ld_decay
tmp <- LDscn$ld_struct$by_chr$Chr1
tmp$edges <- tmp$edges[sample(.N, 10000)]

LDscn$ld_struct$by_chr <-  NULL
LDscn$ld_struct$by_chr$Chr1 <- tmp

usethis::use_data(LDscn, overwrite = TRUE)

plot_ld_decay(LDscn,max_dist = 5e5)
plot(LDscn$decay)
plot(LDscn$ld_struct)


map[,ld_w:=ld_w_0.9]
plot_ldscn_manhattan(LDscn,SNP_res = sim_ex$map,compute_ld_w = FALSE,y_vars = c("C_mean","ld_w"))

LDscn$decay
saveRDS(ld_w_0.9,"ld_w_0.9.rds")
#saveRDS(LDscn, "LDscn_res.rds")

map2 <- add_consistency_to_map(sim_ex$map, consistency_obj = LDscn$consistency)
map2[,ld_w:=compute_ld_w(LDscn$ld_str, LDscn$decay, 0.9,r2_lower_lim = 0.03)]

map2[,plot(-log10(lfmm_q))]

ORs_tbl <- detect_or(q_vals=map2[,.(C_mean)],
                     ld_struct=LDscn$ld_str,
                     decay_obj=LDscn$decay,
                     sign_th=0.05,
                     sign_if = "greater",
                     rho_d=0.97,
                     rho_ld=0.999,
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

