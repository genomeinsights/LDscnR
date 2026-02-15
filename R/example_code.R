library(SNPRelate)
library(dplyr)
library(data.table)



data(sim_example)

map <- copy(sim_example$map)
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


draws <-  ld_rho_draws(ld_struct=ld_str,
                         decay_obj=decay,
                         F_vals=map[,.(emx_F,lfmm_F)],
                         q_vals=map[,.(emx_q,lfmm_q)],
                         n_inds=n_inds,
                         SNP_ids=map$marker,
                         n_rho = 10,
                         rho_w_lim = list(min=0.8,max=0.99),
                         n_or_draws = 25,
                         rho_d_lim=list(min=0.5,max=0.999),
                         rho_ld_lim=list(min=0.9,max=0.999),
                         alpha_lim=list(min=1.31,max=4),
                         lmin_lim=list(min=1,max=10),
                         n_cores = 8,
                         seed = NULL)

C_obj <- consistency_score(draws,combine=TRUE)
C_obj$combined

map2 <- add_consistency_to_map(map, consistency_obj = C_obj)

ORs_tbl <- detect_or(q_vals=map2[,.(C_mean)],
          ld_struct=ld_str,
          decay_obj=decay,
          sign_th=0.05,
          sign_if = "greater",
          rho_d=0.95,
          rho_ld=0.95,
          l_min = 1,
          ret_table = TRUE)

#map2[marker %in% unlist(ORs_tbl$ORs)]

map2 <- merge(map2,
              ORs_tbl,
              by.x = "marker",
              by.y = "SNP",
              all.x = TRUE)

unique_ORs <- unique(na.omit(map2$OR_id))

map2[, OR_factor := factor(OR_id, levels = unique_ORs)]


layout <- prep_manhattan(
  map2[, .(bp = Pos, Chr, marker, C_log = -log10(1 - C_mean + 1e-8),C=C_mean,OR_factor)],
  spacer = 0,
  chr_cols = c("white", "grey50")
)

don      <- copy(layout$data)
axisdf  <- layout$axis
rect_dt <- layout$rect

## y aesthetic
don[, yval := C]



ggplot() +
  ## chromosome background
  geom_rect(
    data = rect_dt,
    aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
    fill = rect_dt$col,
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  ## points
  geom_point(
    data = don,
    aes(BPcum, yval, color = OR_factor),
    size = 1,
    alpha = 1
  ) +
  #scale_color_manual(values =c("black",rep(col_vector,3)),guide="none")+
  scale_shape_manual(values =c(20,3),guide="none")+
  scale_x_continuous(
    breaks = axisdf$center,
    labels = axisdf$Chr,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(x = NULL, y = y) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 0),
    aspect.ratio = 0.25
  )



