library(SNPRelate)
library(dplyr)
library(data.table)

#sim_ex <- readRDS("../LD-scaling-genome-scans/results_sim/c1_V0.5_rep4.rds")

data(sim_ex)
tmp <- readRDS("../LD-scaling-genome-scans/parsed_data/chr1_V0.5_c1_rep4.rds")
keep_125 <- readRDS("../LD-scaling-genome-scans/data/keep_inds_125.rds")
env <- tmp$env$env[keep_125]
map <- copy(sim_ex$SNP_res[,.(Chr, Pos,marker,type, max_LD_with_QTN,bp_to_focal_QTN,focal_QTN,emx_F,lfmm_F,emx_q,lfmm_q)])
GTs <- sim_ex$GTs

n_inds <- nrow(GTs)
map[,lm_F:=apply(GTs,2,function(x)summary(lm(x~env))$fstatistic[1])]
gif <- map[,median(lm_F)/qf(0.5,1,n_inds-2)]
map[,lm_F:=lm_F/gif]
map[,lm_q:=p.adjust(pf(lm_F,1,n_inds-2,lower.tail = FALSE),"fdr")]
map[,emx_q:=p.adjust(pf(emx_F,1,n_inds-2,lower.tail = FALSE),"fdr")]
map[,lfmm_q:=p.adjust(pf(lfmm_F,1,n_inds-2,lower.tail = FALSE),"fdr")]
map[,plot(-log10(lm_q))]

gds_path = tempfile(fileext = ".gds")
gds <- create_gds_from_geno(geno = sim_example$GTs, map=sim_example$map,gds_path)
decay <- ld_decay(gds)
ld_str <- compute_ld_structure(gds,slide_win_ld = 1000,n_cores = 8)



scan <- ld_scan(ld_str,
                decay,
                SNP_ids = sim_example$map$marker,
                F_vals = sim_example$map[,.(emx_F,lfmm_F)],
                rho_w = 0.95,
                n_inds = 125,
                full = TRUE)

plot(scan,method = "lfmm_F")
plot(scan,method = "emx_F")


LDscn <- LDscn_pipeline(gds = gds,
                        F_vals=map[,.(emx_F,lfmm_F,lm_F)],
                        q_vals=map[,.(emx_q,lfmm_q,lm_q)])



map2 <- add_consistency_to_map(map, consistency_obj = LDscn$consistency)
map2[,ld_w:=compute_ld_w(ld_str, decay, 0.9)]

ORs_tbl <- detect_or(q_vals=map2[,.(C_mean)],
          ld_struct=ld_str,
          decay_obj=decay,
          sign_th=0.05,
          sign_if = "greater",
          rho_d=0.9,
          rho_ld=0.9,
          l_min = 1,
          ret_table = TRUE)


ORs_tbl[!duplicated(SNP),table(OR_id,method)]

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
  map2[, .(bp = Pos, Chr, marker, C_log = -log10(1 - C_mean + 1e-8),C_mean,lfmm_q=-log10(lfmm_q),type,OR_factor,lfmm_q_C,emx_q=-log10(emx_q),lm_q=-log10(lm_q),lfmm_F_prime_C,emx_F_prime_C,lm_F_prime_C,ld_w)],
  spacer = 0,
  chr_cols = c("white", "grey50")
)
don      <- copy(layout$data)
axisdf  <- layout$axis
rect_dt <- layout$rect


to_plot <- list(y=c("lfmm_q","emx_q","lm_q","lfmm_F_prime_C","emx_F_prime_C","lm_F_prime_C","C_mean","ld_w"),
     lab=c("-log10(q) | LFMM","-log10(q) | EMX","-log10(q) | lm","LFMM´_C","EMX´_C","lm´_C","Joint_C","ld_w"),
     th=c(-log10(0.05),-log10(0.05),-log10(0.05),0.05,0.05,0.05,0.05,NULL))

plots <- lapply(1:8,function(x){
  don[, yval := get(to_plot$y[x])]
  don[is.na(yval),yval:=0]

  ggplot()+
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
      data = don[OR_factor=="ns"],
      aes(BPcum, yval),
      size = 1,
      col="grey50",
      alpha = 1
    ) +
    geom_point(
      data = don[OR_factor!="ns"],
      aes(BPcum, yval, color = OR_factor),
      size = 1,
      alpha = 1
    ) +
    geom_point(
      data = don[type=="QTN" & OR_factor!="ns"],
      aes(BPcum, yval, color = OR_factor),
      size = 5,
      shape= 3,
      alpha = 1
    ) +
    geom_hline(yintercept = to_plot$th[x],linetype=2, col="grey30") +
    scale_x_continuous(
      breaks = axisdf$center,
      labels = axisdf$Chr,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(x = NULL, y = to_plot$lab[x]) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 90),
      aspect.ratio = 0.25
    )
})

(plots[[1]] / plots[[2]] / plots[[3]] / plots[[4]]) | (plots[[5]] / plots[[6]] / plots[[7]] / plots[[8]])
