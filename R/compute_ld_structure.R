#' Compute Chromosome-wise LD Structure
#'
#' Estimates linkage disequilibrium (LD) structure and LD-decay
#' parameters for each chromosome and (i)ntegrated LD (D)ecay (P)ersistence (I)ndex (iDPI) and optionally compresses and saved per-SNP LD
#' data for downstream LD-based summary statistics.
#'
#' The workflow consists of:
#' \enumerate{
#'   \item Estimation of background LD from inter-chromosomal pairs.
#'   \item Sliding-window estimation of LD-decay parameters per chromosome.
#'   \item Robust summarization of decay parameters.
#'   \item Estimation of (i)ntegrated LD (D)ecay (P)ersistence (I)ndex (iDPI)
#'   \item Optional saving of LD-data to files and histogram compression for downstream estimation of ld_w
#' }
#'
#' LD decay is modeled as:
#' \deqn{r^2(d) = b + \frac{c - b}{1 + a d}}
#'
#' where:
#' \itemize{
#'   \item \eqn{a} is the decay rate,
#'   \item \eqn{b} is background LD,
#'   \item \eqn{c} is short-distance LD,
#'   \item \eqn{d} is physical distance (bp).
#' }
#'
#' @param gds GDS file handle containing genotype data.
#' @param el_data_folder Optional path to folder where edge-list of LD values are saved. If `"NULL"`, data not saved.
#' @param q Quantile used for LD-decay fitting.
#' @param n_sub_bg Number of SNPs used for background LD estimation.
#' @param n_win_decay Number of sliding windows per chromosome for decay fitting.
#' @param overlap Proportion of overlap between consecutive decay windows.
#' @param prob_robust Central proportion of windows retained for robust decay estimation.
#' @param keep_rho_hist Whether to save histograms for estimating ld_w values.
#' @param keep_el Whether to save edge-lists directly. Useful for benchmarking smaller simulated data.
#' @param target_dist_bins_for_decay Number of distance bins for decay fitting.
#' @param k_max Slide window for LD-estimation.
#' @param n_rho_bins_hist Number of distance bins for saved LD histograms.
#' @param rho Maximum distance for pair-wise LD-values relative to decay rate.
#' @param r2_unit Bin width for r² values in histograms.
#' @param cores Number of CPU cores used for parallel computations.
#'
#' @return An object of class `"ld_structure"` containing:
#' \describe{
#'   \item{by_chr}{Per-chromosome LD structure, including histograms and decay parameters.}
#'   \item{decay_sum}{Genome-wide table of robust decay parameters.}
#'   \item{params}{List of parameters used for computation.}
#' }
#'
#' @export
compute_ld_structure <- function(
    gds,
    el_data_folder=NULL,
    ## for LD-decay and bg
    q = 0.95,
    ## for bg
    n_sub_bg = 5000,
    ## for decay
    n_win_decay = 20,
    overlap = 0.5,
    prob_robust = 0.95,
    target_dist_bins_for_decay = 60,
    max_pairs = 5000,
    n_strata = 20,
    keep_el = FALSE,
    ## for ld_int
    n_bins_ld_int = 20,
    ## only for initial filtering of el
    rho = 0.99,
    k_max=1000,
    cores = 1
) {


  if(!is.null(el_data_folder)) if(!dir.exists(el_data_folder)) dir.create(el_data_folder)


  ids  <- .read_gds_ids(gds)
  chrs <- unique(ids$snp_chr)

  b <- estimate_background_ld(gds, n_sub = n_sub_bg, q = q)

  out_by_chr <- vector("list", length(chrs))
  names(out_by_chr) <- chrs

  #ch = "Chr1"
  dist_unit_base = 10000
  for (ch in chrs) {

    chr_idx  <- which(ids$snp_chr == ch)
    pos_chr  <- ids$snp_pos[chr_idx]
    snps_chr <- ids$snp_id[chr_idx]

    cat("Processing ", ch, "-- getting edge list")

    el <- get_el(
      gds = gds,
      idx = chr_idx,
      slide_win_ld = k_max,
      cores = cores,
      by_chr = TRUE,
      symmetric = TRUE,
      edge_symmetry = FALSE
    )

    data.table::setkey(el, SNP)

    window_size <- (max(pos_chr) - min(pos_chr)) / n_win_decay
    step_size   <- window_size * overlap

    # ---- LD-decay
    decay <- suppressWarnings(
      estimate_decay_chr(
        el = el[d<window_size],
        b=b,
        window_size = window_size,
        step_size   = step_size,
        max_pairs=max_pairs,
        n_strata = n_strata,
        q = q,
        cores = cores
      )
    )

    #decay_dt <- decay
    decay[,contrast:=c - b]

    decay[, regime := ifelse(
      contrast < 0.05,
      "weak",
      "structured"
    )]

    decay_sum <- summarize_decay(decay[regime=="structured",])

    decay_sum[,chr_size := max(pos_chr)]

    if (is.null(decay_sum)) next

    decay_sum[,n_w_used := sum(na.omit(decay$regime == "structured"))]
    decay_sum[, Chr := ch]

    d_star <- d_from_rho(decay_sum$a, rho)
    el <- el[d<d_star]

    spacing <- diff(sort(pos_chr))
    q_spacing <- quantile(spacing, 0.75)
    k_star <- ceiling(d_star / q_spacing)

    max_d_available <- max(el$d)
    frac_covered <- signif(max_d_available / d_star,3)
    exceeds_kmax <- k_star > k_max

    cat(" -- frac_covered =", frac_covered)
    cat(" -- k_star =", k_star)

    el[, dist_idx := floor(d / dist_unit_base)]
    #el[,length(unique(dist_idx))]

    shell_dt <- el[, .(
      med_r2 = median(r2)
    ), by = .(SNP, dist_idx)]


    #ld_int_vect <- compute_ld_int(el, snps_chr, decay_sum$a, n_bins=n_bins_ld_int)

    if(!is.null(el_data_folder)) saveRDS(el,paste0(el_data_folder,ch,".rds")); cat(" -- Saving el")

    # ---- collect data
    out_by_chr[[ch]] <- list(
      snp_ids     = snps_chr,
      el          = if(keep_el) el else NULL,
      shell_dt    = shell_dt,
      #ld_int_vect = ld_int_vect,
      decay       = decay,
      decay_sum   = decay_sum
    )

    cat(" -- done\n")
  }

  decay_sum <- rbindlist(lapply(out_by_chr,function(x) x$decay_sum))

  d_mod <- rlm(log(a) ~ log(chr_size), data = decay_sum)

  # decay_sum[,a_pred := exp(predict(d_mod,decay_sum))]
  # decay_sum[,c_pred := median(c)]
  # s_sum <- decay_sum[1,]
  # apply(decay_sum, function(s_sum){
  #
  #   ld_int <- compute_ld_int_adaptive(out_by_chr$Chr1$el,snp_ids = out_by_chr$Chr1$snp_ids,s_sum$a_pred,delta = 5,min_pairs = 0 )
  #   snp_ids <- out_by_chr$Chr1$snp_ids
  #   plot(ld_int)
  #   # bin width based on LD decay scale
  #   shell_dt <- out_by_chr$Chr1$shell_dt
  #   shell_dt[,dist_idx_old := dist_idx]
  #
  #   dist_unit <- 1 / s_sum$a_pred
  #
  #   shell_dt[, dist_idx := floor(d_val / dist_unit)]
  #
  #   #el[,length(unique(dist_idx))]
  #   # shell medians
  #   shell_dt <- shell_dt[, .(
  #     med_r2 = if (.N >= min_pairs) median(med_r2) else NA_real_
  #   ), by = .(SNP, dist_idx)]
  #
  #   shell_dt <- shell_dt[!is.na(med_r2)]
  #
  #   # distance value of each shell
  #   shell_dt[, d_val := dist_idx * dist_unit]
  #
  #   # theoretical decay weight
  #   shell_dt[, weight := (a / (1 + a * d_val)^2) * dist_unit]
  #
  #   # normalize weights per SNP (important for window edges)
  #   shell_dt[, weight := weight / sum(weight), by = SNP]
  #
  #   # weighted median integration
  #   LD_dt <- shell_dt[, .(
  #     LD_int = weighted_median(med_r2, weight)
  #   ), by = SNP]
  #
  #   result <- setNames(rep(NA_real_, length(snp_ids)), snp_ids)
  #   result[LD_dt$SNP] <- LD_dt$LD_int
  #
  #   plot(LD_dt$LD_int,ld_int)
  #
  # })

  out <- list(
    by_chr    = out_by_chr,
    decay_sum = decay_sum,
    params    = list( q = q,
                      n_sub_bg = n_sub_bg,
                      n_win_decay = n_win_decay,
                      overlap = overlap,
                      el_data_folder=el_data_folder,
                      rho = rho,
                      n_strata = n_strata,
                      max_pairs = max_pairs,
                      keep_el   = keep_el,
                      n_bins_ld_int = n_bins_ld_int,
                      prob_robust = prob_robust,
                      target_dist_bins_for_decay = target_dist_bins_for_decay,
                      cores = cores)
  )

  class(out) <- "ld_structure"

  out
}

#' Summarize LD Decay Parameters
#'
#' Computes robust summary estimates of LD-decay parameters across
#' sliding windows within a chromosome.
#'
#' Windows with extreme decay rates are excluded based on symmetric
#' quantile trimming of parameter \eqn{a}.
#'
#' @param decay_dt Data table containing window-wise decay parameters.
#' @param prob_robust Central proportion of windows retained
#'   (e.g. 0.95 keeps the central 95% of \eqn{a} values).
#'
#' @return A \code{data.table} with robust median estimates of
#'   \code{a}, \code{b}, and \code{c}.
#'
#' @export
summarize_decay <- function(decay_dt, prob_robust = 0.95) {

  decay_valid <- decay_dt[!is.na(a) & a > 0 & !is.na(c)]

  if (nrow(decay_valid) < 5)
    return(NULL)

  # compute symmetric trimming bounds
  alpha <- (1 - prob_robust) / 2

  q_lo <- quantile(decay_valid$a, alpha, na.rm = TRUE)
  q_hi <- quantile(decay_valid$a, 1 - alpha, na.rm = TRUE)

  decay_trim <- decay_valid[a >= q_lo & a <= q_hi]

  if (nrow(decay_trim) < 3)
    return(NULL)

  data.table(
    c = median(decay_trim$c, na.rm = TRUE),
    a = median(decay_trim$a, na.rm = TRUE),
    b = decay_trim$b[1]
  )
}

subsample_pairs_for_decay <- function(sub,
                                      max_pairs = 5000,
                                      n_strata = 20) {


  if (nrow(sub) <= max_pairs) max_pairs <- nrow(sub)
  # log-distance strata
  log_d <- log(sub$d)
  breaks <- seq(min(log_d), max(log_d), length.out = n_strata + 1)

  sub[, strata := cut(log_d, breaks = breaks, include.lowest = TRUE)]

  target_per_stratum <- floor(max_pairs / n_strata)

  sub_sampled <- sub[, {
    if (.N > target_per_stratum) {
      .SD[sample(.N, target_per_stratum)]
    } else {
      .SD
    }
  }, by = strata]


  return(sub_sampled)
}

#' Estimate LD Decay for a Chromosome
#'
#' Fits LD-decay curves in sliding genomic windows for a single chromosome.
#'
#' SNP density is adaptively reduced to maintain approximately constant
#' SNP counts per window before decay fitting.
#'
#' @param gds GDS file handle.
#' @param chr_idx SNP indices for the chromosome.
#' @param snp_pos Physical SNP positions (bp).
#' @param b Background LD.
#' @param window_size Physical window size (bp).
#' @param step_size Step size (bp) between windows.
#' @param q Quantile of r² used for decay fitting.
#' @param target_dist_bins_for_decay Number of distance bins used in fitting.
#' @param cores Number of CPU cores.
#' @param n_snps_for_decay Target SNP count per window after thinning.
#'
#' @return A \code{data.table} with window-specific decay parameters.
#'
#' @export
estimate_decay_chr <- function(el,
                               b,
                               window_size,
                               step_size,
                               q = 0.95,
                               max_pairs=5000,
                               n_strata = 20,
                               cores = 1) {



  min_pos <- min(c(el$pos1, el$pos2))
  max_pos <- max(c(el$pos1, el$pos2))

  starts <- seq(min_pos, max_pos - window_size, by = step_size)
  ends   <- starts + window_size
  #i <- 1
  decay <- suppressWarnings(rbindlist(lapply(seq_along(starts), function(i) {

    sub <- el[
      pos1 >= starts[i] & pos1 < ends[i] &
        pos2 >= starts[i] & pos2 < ends[i]
    ]

    if(!nrow(sub)>100) return(NULL)

    sub <- subsample_pairs_for_decay(
      sub,max_pairs = max_pairs, n_strata = n_strata
    )

    coefs <- tryCatch(
      coef_ld_dec(dt_strata=sub, q = q,  b = b),
      error = function(e) NULL
    )

    if (is.null(coefs)) {
      data.table(start = starts[i], end = ends[i])
    } else {
      data.table(start = starts[i], end = ends[i],coefs)
    }
  }),fill=TRUE,use.names = TRUE))

  return(decay)
}


estimate_background_ld <- function(gds,
                                   idx=NULL,
                                   n_sub = 5000,
                                   q = 0.95) {


  ids <- .read_gds_ids(gds)

  if (is.null(idx)) idx <- seq_along(ids$snp_id)

  n_snps <- length(idx)
  if (n_snps < 2L) {
    stop("Not enough SNPs to estimate background LD.")
  }

  chr_vec <- ids$snp_chr[idx]
  chr_levels <- unique(chr_vec)

  if (length(chr_levels) < 2L) {
    stop("Need at least two chromosomes to estimate background LD.")
  }


  ## proportional sampling across chromosomes
  snp_pool <- unlist(lapply(chr_levels, function(ch) {

    ix <- which(chr_vec == ch)
    n_ch <- length(ix)

    sample(ix,size = min(n_ch, ceiling(n_sub * n_ch / n_snps)),replace = FALSE)

  }))

  if (length(unique(chr_vec[snp_pool])) < 2L) {
    stop("Sampled SNPs fall on a single chromosome; increase n_sub.")
  }

  ld <- snpgdsLDMat(
    gds,
    snp.id = ids$snp_id[idx][snp_pool],
    method = "r",
    slide = -1,
    verbose = FALSE
  )

  r2 <- ld$LD^2
  chr_sub <- chr_vec[snp_pool]

  ## upper triangle only
  inter_idx <- which(outer(chr_sub, chr_sub, FUN = "!=") & upper.tri(r2))
  r2_inter <- r2[inter_idx]

  if (length(r2_inter) == 0L) {
    stop("No inter-chromosomal SNP pairs found.")
  }

  b <- quantile(r2_inter, probs = q, na.rm = TRUE)

  message("Background LD: ", round(b, 3))
  return(b)
}



#' Fit LD Decay Model
#'
#' Fits the nonlinear decay model
#' \deqn{r^2(d) = b + \frac{c - b}{1 + a(d - d_0)}}
#' to binned LD quantiles.
#'
#' @param dt_strata LD edge data table with strata.
#' @param q Quantile for LD binning.
#' @param b Background LD.
#'
#' @return Named vector with parameters \code{c}, \code{a}, \code{d0}.
#'
#' @export
coef_ld_dec <- function(dt_strata,
                        q = 0.95,
                        b) {

  if (nrow(dt_strata) < 100) return(NULL)
  #dt_strata <- sub
  agg <- dt_strata[, .(
    d_mid = exp(mean(log(d))),   # geometric midpoint
    r2_q  = quantile(r2, q, na.rm = TRUE),
    n     = .N
  ), by = strata]


  if (nrow(agg) < 10) return(NULL)

  # Initial parameter guesses
  c_start <- max(agg$r2_q, na.rm = TRUE)
  a_start <- 1 / median(agg$d_mid)

  fit <- tryCatch(
    nls(
      r2_q ~ b + (c - b)/(1 + a * d_mid),
      data = agg,
      start = list(c = min(max(agg$r2_q), 1), a = a_start),
      algorithm = "port",
      lower = c(c = b, a = 0),
      upper = c(c = 1, a = Inf),
      weights = sqrt(n),
      control = nls.control(warnOnly = TRUE)
    ),
    error = function(e) NULL
  )


  if (is.null(fit)) return(NULL)

  coefs <- coef(fit)

  return(data.table(c=coefs["c"],a=coefs["a"],b,agg=list(agg),raw=list(dt_strata)))

}



#' Convert Relative LD Threshold to Distance
#'
#' Computes the physical distance corresponding to a relative LD
#' threshold \eqn{\rho}.
#'
#' @param a Decay rate.
#' @param rho Relative LD threshold.
#' @param d0 Offset parameter.
#'
#' @return Distance in base pairs.
#'
#' @export
d_from_rho <- function(a, rho){
  rho / (a * (1 - rho))
}

ld_from_rho <- function(b, c = 1, rho){
  b + (c - b) * (1 - rho)
}


parallel_apply <- function(X, FUN, cores = 1) {

  if (cores > 1 && .Platform$OS.type != "windows") {
    parallel::mclapply(X, FUN, mc.cores = cores)
  } else {
    lapply(X, FUN)
  }
}


thin_to_k_target <- function(pos, d_star, k_target = 1000) {

  pos <- sort(pos)

  # current effective density
  lambda_current <- 1 / median(diff(pos))

  # target density
  lambda_target <- k_target / d_star

  # if already acceptable
  if (lambda_current <= lambda_target) {
    return(seq_along(pos))
  }

  # target spacing
  target_spacing <- 1 / lambda_target

  keep <- integer()
  last_kept <- -Inf

  for (i in seq_along(pos)) {
    if (pos[i] - last_kept >= target_spacing) {
      keep <- c(keep, i)
      last_kept <- pos[i]
    }
  }

  keep
}

#' Print Method for ld_structure Objects
#'
#' Displays a concise summary of the LD structure object.
#'
#' @param x An object of class `"ld_structure"`.
#' @param ... Additional arguments (unused).
#'
#' @export
print.ld_structure <- function(x, ...) {

  if (!inherits(x, "ld_structure"))
    stop("Object must be of class 'ld_structure'.")

  cat("LD Structure Object\n")
  cat("-------------------\n")

  # Chromosomes
  chrs <- names(x$by_chr)
  n_chr <- length(chrs)
  cat("Chromosomes:", n_chr, "\n")

  # Total SNPs
  n_snps <- sum(
    vapply(x$by_chr, function(chr)
      length(chr$snp_ids), numeric(1)),
    na.rm = TRUE
  )
  cat("Total SNPs:", n_snps, "\n\n")

  # Decay summary
  if (!is.null(x$decay_sum) && nrow(x$decay_sum) > 0) {

    cat("Robust LD-decay parameters (median across chromosomes):\n")

    if ("a" %in% names(x$decay_sum))
      cat("  a =", signif(median(x$decay_sum$a, na.rm = TRUE), 3), "\n")

    if ("c" %in% names(x$decay_sum))
      cat("  c =", signif(median(x$decay_sum$c, na.rm = TRUE), 3), "\n")

    if ("b" %in% names(x$decay_sum))
      cat("  b =", signif(median(x$decay_sum$b, na.rm = TRUE), 3), "\n")

    if ("k_star" %in% names(x$decay_sum))
      cat("  k* (median) =",
          signif(median(x$decay_sum$k_star, na.rm = TRUE), 3), "\n")
  }

  cat("\nParameters used:\n")

  if (!is.null(x$params)) {
    for (nm in names(x$params)) {
      cat(" ", nm, "=", x$params[[nm]], "\n")
    }
  }

  cat("\nLD histograms stored per chromosome (not printed).\n")

  invisible(x)
}


#' Plot Method for ld_structure Objects
#'
#' Visualizes LD-decay curves or decay-rate distributions.
#'
#' @param x An object of class `"ld_structure"`.
#' @param type Plot type: `"decay"` or `"a_dist"`.
#' @param use Decay summary to display: `"robust"` or `"median"`.
#' @param alpha Transparency for window-specific curves.
#' @param ... Additional graphical parameters.
#'
#' @export
plot.ld_structure <- function(x,
                              type = c("decay", "a_dist"),
                              use = c("robust", "median"),
                              alpha = 0.4,
                              ...) {

  if(!is.null(dev.list())) dev.off()
  if (!inherits(x, "ld_structure"))
    stop("Object must be of class 'ld_structure'.")

  type <- match.arg(type)
  #use  <- match.arg(use)

  chrs <- names(x$by_chr)


  if (type == "decay") {

    n_chr <- length(chrs)
    ncol = ceiling(sqrt(n_chr))
    nrow_plot <- ceiling(n_chr / ncol)

    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))

    par(mfrow = c(nrow_plot, ncol))
    par(mar = c(2, 2, 2, 1))

    for (ch in chrs) {

      z     <- x$by_chr[[ch]]
      decay <- z$decay

      if (is.null(decay) || nrow(decay) < 5) {
        plot.new()
        next
      }

      b <- as.numeric(na.omit(decay$b)[1])

      # Distance grid
      max_d <- max(decay$end, na.rm = TRUE)
      d <- seq(0, max_d, length.out = 200)

      plot(NULL,
           xlim = c(0, max_d),
           ylim = c(0, 1),
           xlab = "",
           ylab = "",
           main = ch)

      abline(h = b, lty = 2, col = "grey40")

      # Plot window-specific fits
      for (i in seq_len(nrow(decay))) {

        a_i  <- decay$a[i]
        c_i  <- decay$c[i]
        d0_i <- if ("d0" %in% names(decay)) decay$d0[i] else 0

        if (is.na(a_i) || a_i <= 0) next

        r2 <- b + (c_i - b) /
          (1 + a_i * pmax(d - d0_i, 0))

        lines(d, r2,
              lwd = 1.5,
              col = scales::alpha("steelblue", alpha))
      }

      # Representative curve
      #y <- use[1]
      lapply(use,function(y){
        summary_use <- z$summary[[y]]

        if (!is.null(summary_use)) {

          a_rep  <- summary_use$a
          c_rep  <- summary_use$c
          d0_rep <- summary_use$d0

          r2_rep <- b + (c_rep - b) /
            (1 + a_rep * pmax(d - d0_rep, 0))

          lines(d, r2_rep,
                lwd = 3,
                col = c("firebrick","salmon")[use==y])
        }
      })

    }
  }

  if (type == "a_dist") {

    a_vals <- rbindlist(lapply(x$by_chr, function(z) {
      z$decay[!is.na(a) & a > 0, .(a)]
    }))

    hist(a_vals$a,
         breaks = 40,
         col = "grey70",
         border = "white",
         main = "Distribution of LD-decay rate (a)",
         xlab = "a")

    abline(v = median(a_vals$a, na.rm = TRUE),
           col = "firebrick",
           lwd = 2)
  }

  invisible(x)
}

compute_ld_w <- function(
    ld_structure,
    rho = 0.95,
    cores = 1
) {

  ld_w <- unlist(parallel_apply(ld_structure$by_chr, function(chr_obj) {

    a <- chr_obj$decay_sum$a
    b <- chr_obj$decay_sum$b

    d_window <- d_from_rho(a, rho)

    if(is.null(chr_obj$el)) stop("No edge list present")

    ld_w <- chr_obj$el[d<d_window,.(r2_median=median(r2)),by=SNP]

    ld_w[match(chr_obj$snp_ids,ld_w$SNP),r2_median]

  }, cores = cores))

  return(ld_w)
}


weighted_median <- function(x, w) {
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  w <- w / sum(w)
  cs <- cumsum(w)
  x[which(cs >= 0.5)[1]]
}


compute_ld_int <- function(el, snp_ids, a, n_bins) {

  # restrict to integration radius
  if (nrow(el) == 0)
    return(setNames(rep(NA_real_, length(snp_ids)), snp_ids))

  # distance grid
  max_d     <- max(el$d)
  dist_unit <- signif(max_d / n_bins, 2)

  el[, dist_idx := floor(d / dist_unit)]

  # shell medians per SNP
  shell_dt <- el[, .(
    med_r2 = median(r2)
  ), by = .(SNP, dist_idx)]

  shell_dt[, d_val := dist_idx * dist_unit]

  shell_dt[, weight := (a / (1 + a * d_val)^2) * dist_unit]

  LD_dt <- shell_dt[, .(
    LD_int = weighted_median(med_r2, weight)
  ), by = SNP]



  #plot(LD_dt[,LD_int],LD_dt_el[,LD_int])


  result <- setNames(rep(NA_real_, length(snp_ids)), snp_ids)
  result[LD_dt$SNP] <- LD_dt$LD_int

  result
}

# n_bins_target <- 30
#a=0.000001
compute_ld_int_adaptive <- function(el, snp_ids, a, ld_resolution=2, min_shell_n = 3) {

  if (nrow(el) == 0)
    return(setNames(rep(NA_real_, length(snp_ids)), snp_ids))

  # bin_frac_target <- max(el$d) * a / n_bins_target
  #
  # dist_unit <- (1 / a) * bin_frac
  # el[, dist_idx := floor(d / dist_unit)]

  max_d <- max(el$d)
  n_bins <- ceiling(a * max_d / ld_resolution)

  #el[,length(unique(dist_idx))]

  shell_dt <- el[, .(
    med_r2 = if (.N >= min_shell_n) median(r2) else NA_real_
  ), by = .(SNP, dist_idx)]

  shell_dt <- shell_dt[is.finite(med_r2)]
  shell_dt[, d_val := dist_idx * dist_unit]
  shell_dt[, weight := (a / (1 + a * d_val)^2) * dist_unit]

  shell_dt[, weight := weight / sum(weight), by = SNP]
  #shell_dt[SNP=="Chr1:10559176"]

  LD_dt <- shell_dt[, .(
    LD_int = weighted_median(med_r2, weight)
  ), by = SNP]

  #LD_dt[,plot(LD_int)]

  #LD_dt_01 <- copy(LD_dt)
  #LD_dt_1 <- copy(LD_dt)

  #plot(LD_dt$LD_int,LD_dt_24$LD_int)

  #plot(map[match(LD_dt$SNP,marker),max_LD_with_QTN],LD_dt$LD_int)

  result <- setNames(rep(NA_real_, length(snp_ids)), snp_ids)
  result[LD_dt$SNP] <- LD_dt$LD_int
  result
}

compute_ld_int_adaptive <- function(el, snp_ids, a, delta = 1, min_pairs = 3) {

  if (nrow(el) == 0)
    return(setNames(rep(NA_real_, length(snp_ids)), snp_ids))

  # bin width based on LD decay scale
  dist_unit <- delta / a

  # assign linear bins
  el[, dist_idx := floor(d / dist_unit)]
  #el[,length(unique(dist_idx))]
  # shell medians
  shell_dt <- el[, .(
    med_r2 = if (.N >= min_pairs) median(r2) else NA_real_
  ), by = .(SNP, dist_idx)]

  shell_dt <- shell_dt[!is.na(med_r2)]

  # distance value of each shell
  shell_dt[, d_val := dist_idx * dist_unit]

  # theoretical decay weight
  shell_dt[, weight := (a / (1 + a * d_val)^2) * dist_unit]

  # normalize weights per SNP (important for window edges)
  shell_dt[, weight := weight / sum(weight), by = SNP]

  # weighted median integration
  LD_dt <- shell_dt[, .(
    LD_int = weighted_median(med_r2, weight)
  ), by = SNP]

  result <- setNames(rep(NA_real_, length(snp_ids)), snp_ids)
  result[LD_dt$SNP] <- LD_dt$LD_int

  result
}
