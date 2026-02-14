library(SNPRelate)
library(dplyr)
library(data.table)



data(sim_example)

map <- sim_example$map
GTs <- sim_example$GTs
n_inds <- nrow(GTs)
map[,emx_q:=p.adjust(pf(emx_F,1,n_inds-2,lower.tail = FALSE),"fdr")]
map[,lfmm_q:=p.adjust(pf(lfmm_F,1,n_inds-2,lower.tail = FALSE),"fdr")]


gds_path = tempfile(fileext = ".gds")
gds <- create_gds_from_geno(geno = sim_example$GTs, map=sim_example$map,gds_path)
decay <- ld_decay(gds)
ld_str <- compute_ld_structure(gds,slide_win_ld = 1000,n_cores = 8)

ld_w  <- compute_ld_w(ld_str, decay, 0.9)

scan <- ld_scan(ld_str,
                decay,
                SNP_ids = sim_example$map$marker,
                F_vals = sim_example$map[,.(emx_F,lfmm_F)],
                rho_w = 0.95,
                n_inds = 125,
                full = TRUE)

plot(scan,method = "lfmm_F")
plot(scan,method = "emx_F")
#x <- copy(scan)
q_vals     = cbind(sim_example$map[,.(emx_q,lfmm_q)],
                   do.call(cbind,lapply(scan$result,function(x) x$q_prime)))

ORs <- detect_or(ld_struct=ld_str,
                 decay_obj=decay,
                 q_vals = q_vals,
                 alpha=0.05,
                 rho_d=0.95,
                 rho_ld=0.95,
                 l_min = 2)


draws <-  ld_rho_draws(ld_struct=ld_str,
                         decay_obj=decay,
                         F_vals=sim_example$map[,.(emx_F,lfmm_F)],
                         q_vals = sim_example$map[,.(emx_q,lfmm_q)],
                         n_inds=n_inds,
                         SNP_ids=sim_example$map$marker,
                         n_rho = 10,
                         rho_w_lim = list(min=0.8,max=0.99),
                         n_or_draws = 25,
                         rho_d_lim=list(min=0.5,max=0.999),
                         rho_ld_lim=list(min=0.9,max=0.999),
                         alpha_lim=list(min=1.31,max=4),
                         lmin_lim=list(min=1,max=10),
                         n_cores = 8,
                         seed = NULL)

C_scores <- consistency_score(rho_draws_obj = draws)

C_obj <- consistency_score(draws)
C_obj <- combine_consistency(C_obj, method = "geometric")

map2 <- add_consistency_to_map(map, C_obj)

layout <- prep_manhattan(
  map2[, .(bp = Pos, Chr, marker, C_log = -log10(1 - C_mean + 1e-8),C=C_mean)],
  spacer = 0,
  chr_cols = c("white", "grey90")
)

plot_manhattan_gg(layout, y = "C_log")

