library(SNPRelate)
library(dplyr)
library(data.table)


sim_example <- list(
  GTs = data$GTs,
  map = data$SNP_res[,.(Chr, Pos, marker, emx_F, lfmm_F, max_LD_with_QTN)]
)

# usethis::use_data(sim_example, overwrite = TRUE)
# devtools::document()
sim_example()

devtools::load_all()
data(sim_example)

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

x <- copy(scan)
plot(scan,method = "lfmm_F")
plot(scan,method = "emx_F")

ORs <- detect_or(scan_obj=scan,
                 ld_struct=ld_str,
                 decay_obj=decay,
                 alpha=0.05,
                 rho_d=0.95,
                 rho_ld=0.95,
                 l_min = 2)

rho_w_lim = list(min=0,max=0.5)

draws <- list()
for(i in 1:10){
  cat("Draw",i,"..")
  scan <- ld_scan(ld_str,
                  decay,
                  SNP_ids = sim_example$map$marker,
                  F_vals = sim_example$map[,.(emx_F,lfmm_F)],
                  rho_w = runif(1, rho_w_lim$min,rho_w_lim$max),
                  n_inds = 125,
                  full = TRUE)

  draws[[i]] <- or_draws(scan_obj=scan,
                    ld_struct,
                    decay_obj,
                    n_draws = 25,
                    rho_d_lim=list(min=0.5,max=0.999),
                    rho_ld_lim=list(min=0.9,max=0.999),
                    alpha_lim=list(min=1.31,max=4),
                    lmin_lim=list(min=1,max=10))
}

return(draws)

draws <-  ld_rho_draws(ld_struct=ld_str,
                         decay_obj=decay,
                         F_vals=sim_example$map[,.(emx_F,lfmm_F)],
                         n_inds=125,
                         SNP_ids=sim_example$map$marker,
                         n_rho = 10,
                         rho_w_lim = list(min=0,max=0.5),
                         n_or_draws = 25,
                         rho_d_lim=list(min=0.5,max=0.999),
                         rho_ld_lim=list(min=0.9,max=0.999),
                         alpha_lim=list(min=1.31,max=4),
                         lmin_lim=list(min=1,max=10),
                         seed = NULL)

C_scores <- consistency_score(rho_draws_obj = draws)
C_scores$consistency


plot(C_scores$consistency$emx_F$consistency,C_scores$consistency$lfmm_F$consistency)

sim_example$map[,C_emx := C_scores$consistency$emx_F$consistency[match(sim_example$map$marker,C_scores$consistency$emx_F$SNP)]]
sim_example$map[,C_lfmm := C_scores$consistency$emx_F$consistency[match(sim_example$map$marker,C_scores$consistency$lfmm_F$SNP)]]



