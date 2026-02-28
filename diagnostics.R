##### ===== a-stability ===== ######
# gds,
# chr_idx,
# snp_pos=pos_chr,
# b,
n_snps_grid = c(5:10*50,1000,2000)
dist_bins_grid = c(20,40,60,80,100,200,350,500,750,1000,1500)
#overlap = 0.5

library(data.table)
library(ggplot2)
#n_snps <- n_snps_grid[1]
#target_dist_bins <- dist_bins_grid[1]
results <- rbindlist(lapply(n_snps_grid, function(n_snps) {



  # par(mfcol=c(3,4))

  results <- rbindlist(lapply(dist_bins_grid, function(target_dist_bins) {



    decay <- suppressWarnings(
      estimate_decay_chr(
        gds,
        chr_idx,
        snp_pos = pos_chr,
        b,
        window_size = window_size,
        step_size   = step_size,
        q = q,
        target_dist_bins = target_dist_bins,
        n_snps = n_snps,
        cores = cores
      )
    )


    decay_sum <- summarize_decay(decay, prob_robust)

    #decay$decay[,data][[1]][,plot(d_mid,r2_q,main=paste("target bins=",target_dist_bins))]
    # a=median(decay$decay$a, na.rm = TRUE)
    # c=median(decay$decay$c, na.rm = TRUE)
    #decay$decay[,data][[1]][order(d_mid),lines(d_mid, b + ( decay$decay$c[1]  - b) / (1 + decay$decay$a[1]  * d_mid),col="red")]

    data.table(
      n_snps = n_snps,
      target_dist_bins = target_dist_bins,
      median_a = decay_sum$a
    )

  }),fill=TRUE,use.names = TRUE)

}))

# ---- Plot 1: effect of n_snps ----
p1 <- ggplot(results[!duplicated(median_a)],
             aes(x = log(n_snps),
                 y = median_a,
                 color = log(target_dist_bins),
             group=factor(target_dist_bins))) +
  geom_line() +
  geom_point() +
  scale_color_viridis_c(option="turbo")+
  theme_bw(base_size = 14) +
  labs(color = "dist_bins",
       y = "Median a",
       title = "Effect of n_snps on LD decay parameter a")


# ---- Plot 2: effect of dist_bins ----
#results[!duplicated(median_a)]
p2 <- ggplot(results[!duplicated(median_a)],
             aes(x = log(target_dist_bins),
                 y = median_a,
                 color = log(n_snps),
                 group = factor(n_snps))) +
  geom_line() +
  geom_point() +
  scale_color_viridis_c(option="turbo")+
  theme_bw(base_size = 14) +
  labs(color = "n_snps",
       y = "Median a",
       title = "Effect of distance bin resolution on a")

p1 | p2



results <- test_decay_sensitivity
ggplot(results,
       aes(x = n_snps,
           y = target_dist_bins,
           fill = median_a)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw()

dt <- sub



if (k_star > k_max) {
  warning(
    sprintf(
      "Estimated k* = %.0f SNPs exceeds recommended maximum (%d).
Consider thinning SNP density by a factor of ~%.1f
to achieve k* ≈ %d.",
      k_star, k_max, k_star / k_target, k_target
    )
  )
}

thin_to_k_target <- function(pos, d_star, k_target = 1000) {

  L <- max(pos) - min(pos)
  lambda_current <- length(pos) / L
  lambda_target  <- k_target / d_star

  if (lambda_current <= lambda_target) {
    return(seq_along(pos))
  }

  thin_factor <- ceiling(lambda_current / lambda_target)

  keep <- seq(1, length(pos), by = thin_factor)

  keep
}
