library(data.table)
library(ggplot2)
library(patchwork)
ll_files <- list.files("../LD-scaling-genome-scans/results_sim/",full.names = TRUE)
all_files <- all_files[!grepl("PR",all_files)]
files <- c(all_files[grep("c2_V0.5",all_files)],all_files[grep("c1.5_V1",all_files)],all_files[grep("c1_V2",all_files)])

done <- list.files("./PR_out_robust_fill/")
todo <- do.call(rbind,strsplit(files,"/"))
todo <- todo[,ncol(todo)]
done <- list.files("./PR_out_robust_fill/")
rho_bin_int_data <- rbindlist(lapply(files[which(todo %in% done)],function(file){
  tmp <- strsplit(file,"/")[[1]]
  out_file <- tmp[length(tmp)]
  rho_bin_int_data <- readRDS(paste0("./PR_out_robust_fill/",out_file))
}))

dt_sum <- rho_bin_int_data[,
  .(
    AUC_mean = mean(AUC, na.rm = TRUE),
    AUC_sd   = sd(AUC, na.rm = TRUE),
    n        = .N,
    AUC_se   = sd(AUC, na.rm = TRUE) / sqrt(.N),
    AUC_lwr  = mean(AUC, na.rm = TRUE) - 1.96 * sd(AUC, na.rm = TRUE) / sqrt(.N),
    AUC_upr  = mean(AUC, na.rm = TRUE) + 1.96 * sd(AUC, na.rm = TRUE) / sqrt(.N)
  ),
  by = .(AUC_PR, alpha_min, c, V, ld, rho, bins)
]

base_rhow <- dt_sum[
  AUC_PR == "PR-C" & ld == "rho_w",
  .(base_AUC = mean(AUC_mean, na.rm = TRUE)),
  by = .(alpha_min,c)
]

dt_heat <- merge(
  dt_sum[AUC_PR == "PR-C" & ld == "DPI"],
  base_rhow,
  by = c("alpha_min","c"),
  all.x = TRUE
)

dt_heat[, delta_AUC := AUC_mean - base_AUC]
dt_heat[, alpha_lab := factor(1 / 10^(alpha_min),
                              levels = sort(unique(1 / 10^(alpha_min))),
                              labels = paste0("alpha=", sort(unique(1 / 10^(alpha_min)))))]

rho_default  <- 0.99
bins_default <- 1000

bins_default <- dt_sum[ld == "DPI", sort(unique(bins))][which.min(abs(sort(unique(dt_sum[ld == "DPI", bins])) - bins_default))]

dt_default_dpi <- dt_sum[
  AUC_PR == "PR-C" & ld == "DPI" & rho == rho_default & bins == bins_default
][, method := "DPI"]

dt_default_rhow <- dt_sum[
  AUC_PR == "PR-C" & ld == "rho_w"
][, `:=`(method = "rho_w", rho = rho_default, bins = bins_default)]

dt_default <- rbindlist(list(dt_default_dpi, dt_default_rhow), fill = TRUE)

dt_default[, alpha_lab := factor(1 / 10^(alpha_min),
                                 levels = sort(unique(1 / 10^(alpha_min))),
                                 labels = paste0("alpha=", sort(unique(1 / 10^(alpha_min)))))]

pA <- ggplot(
  dt_heat,
  aes(x = bins, y = rho, fill = delta_AUC)
) +
  geom_tile() +
  facet_grid(c ~ alpha_lab) +
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    name = expression(Delta * "AUC")
  ) +
  geom_point(
    data = dt_heat[rho == rho_default & bins == bins_default],
    shape = 21, size = 3, stroke = 1.1, fill = NA, color = "black"
  ) +
  scale_x_continuous(trans = "log10") +
  labs(
    x = "Number of bins (log scale)",
    y = expression(rho[max])
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 1/1.1,
    panel.grid = element_blank(),
    strip.background = element_blank()
  )
pA

#### B ####
pB <- ggplot(
  dt_default,
  aes(x = alpha_lab, y = AUC_mean, group = method, color = method, fill = method)
) +
  geom_ribbon(
    aes(ymin = AUC_lwr, ymax = AUC_upr),
    alpha = 0.15,
    color = NA
  ) +
  facet_grid(. ~ c) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  scale_color_manual(values = c("rho_w" = "#F8766D", "DPI" = "black")) +
  scale_fill_manual(values = c("rho_w" = "#F8766D", "DPI" = "black")) +
  labs(
    x = "Minimum alpha",
    y = "AUC"
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.title = element_blank(),
    legend.position = "top"
  )

pB
p_main <- pA + pB + plot_layout(widths = c(1.8, 1))
p_main


###### C ######
rho_bin_int_data$PR_dat[[1]]
mean_C <- sapply(rho_bin_int_data$C_scores,function(x) length(x$Joint_C>0))
mean_C <- sapply(rho_bin_int_data$C_scores,function(x) mean(x$Joint_C,na.rm=TRUE)*sd(x$Joint_C,na.rm=TRUE))

plot(mean_C,rho_bin_int_data$AUC)

