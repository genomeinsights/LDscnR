#' Compute Chromosome-wise LD Structure
#'
#' Estimates linkage disequilibrium (LD) structure and LD-decay
#' parameters for each chromosome, and constructs per-SNP LD
#' histograms for downstream LD-based summary statistics.
#'
#' The workflow consists of:
#' \enumerate{
#'   \item Estimation of background LD from inter-chromosomal pairs.
#'   \item Sliding-window estimation of LD-decay parameters per chromosome.
#'   \item Robust summarization of decay parameters.
#'   \item Construction of per-SNP LD histograms within an adaptive
#'         LD radius derived from the decay model.
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
#' @param q Quantile used for LD-decay fitting.
#' @param n_sub_bg Number of SNPs used for background LD estimation.
#' @param n_win_decay Number of sliding windows per chromosome for decay fitting.
#' @param overlap Proportion of overlap between consecutive decay windows.
#' @param prob_robust Central proportion of windows retained for robust decay estimation.
#' @param target_dist_bins_for_decay Number of distance bins for decay fitting.
#' @param n_snps_for_decay Target number of SNPs per decay window after thinning.
#' @param n_dist_target_for_hist Number of distance bins for LD histograms.
#' @param eps Relative LD tail tolerance used to derive integration radius.
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
    ## for LD-decay and bg
    q = 0.95,
    ## for bg
    n_sub_bg = 5000,
    ## for decay
    n_win_decay = 20,
    overlap = 0.5,
    prob_robust = 0.95,
    target_dist_bins_for_decay = 200,
    n_snps_for_decay = 500,
    ## for histograms
    n_dist_target_for_hist = 100,
    eps = 0.005,
    r2_unit = 0.001,
    compression = c("hist", "median_only"),
    cores = 1,
    k_max=1000
) {
  #compression = "median_only"
  ids  <- .read_gds_ids(gds)
  chrs <- unique(ids$snp_chr)

  b <- estimate_background_ld(gds, n_sub = n_sub_bg, q = q)

  out_by_chr <- vector("list", length(chrs))
  names(out_by_chr) <- chrs
  #ch = "Chr3"
  for (ch in chrs) {

    cat("Processing ", ch, "-- estimating LD-decay")

    chr_idx  <- which(ids$snp_chr == ch)
    pos_chr  <- ids$snp_pos[chr_idx]
    snps_chr <- ids$snp_id[chr_idx]


    window_size <- (max(pos_chr) - min(pos_chr)) / n_win_decay
    step_size   <- window_size * overlap


    # ---- LD-decay
    decay <- suppressWarnings(
      estimate_decay_chr(
        gds=gds,
        chr_idx=chr_idx,
        snp_pos = pos_chr,
        b=b,
        window_size = window_size,
        step_size   = step_size,
        q = q,
        target_dist_bins_for_decay = target_dist_bins_for_decay,
        n_snps_for_decay = n_snps_for_decay,
        cores = cores
      )
    )


    decay[,contrast:=c - b]
    decay[,regime := ifelse(contrast < 0.02,"weak","structured")]

    decay_sum <- summarize_decay(decay[regime=="structured",], prob_robust)
    decay_sum[,n_w_used := sum(decay$regime == "structured")]

    if (is.null(decay_sum)) next

    decay_sum[, Chr := ch]

    # ---  Histograms

    d_star <- derive_ld_radius(decay_sum$a, eps)

    cat(" -- compressing data to histograms")
    histograms <- build_ld_histograms(
      gds=gds,
      chr_idx   = chr_idx,
      snps_chr  = snps_chr,
      pos_chr   = pos_chr,
      d_star    = d_star,
      r2_unit   = r2_unit,
      n_dist_target = n_dist_target_for_hist,
      k_max     = k_max,
      compression = compression,
      cores     = cores
    )

    # ---- collect data
    out_by_chr[[ch]] <- list(
      snp_ids    = snps_chr,
      hist_obj   = histograms,
      decay      = decay,
      decay_sum  = decay_sum
    )

    cat(" -- done\n")
  }

  out <- list(
    by_chr    = out_by_chr,
    decay_sum = rbindlist(lapply(out_by_chr,function(x) x$decay_sum)),
    params    = list( q = q,
                      n_sub_bg = n_sub_bg,
                      n_win_decay = n_win_decay,
                      overlap = overlap,
                      n_dist_target_for_hist = n_dist_target_for_hist,
                      eps = eps,
                      r2_unit = r2_unit,
                      prob_robust = prob_robust,
                      target_dist_bins_for_decay = target_dist_bins_for_decay,
                      n_snps_for_decay = n_snps_for_decay,
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
estimate_decay_chr <- function(gds,
                               chr_idx,
                               snp_pos,
                               b,
                               window_size,
                               step_size,
                               overlap = 0.5,
                               q = 0.95,
                               target_dist_bins_for_decay = 60,
                               cores = 1,
                               n_snps_for_decay = 500,
                               k_max=1000) {




  thin <- adaptive_thinning(snp_pos, W_LD = window_size, n_snps_for_decay = n_snps_for_decay)
  reduction <- 1-length(thin$keep)/length(chr_idx)

  el <- get_el(
    gds,
    idx = chr_idx[thin$keep],
    slide_win_ld = thin$k_star,
    cores = cores
  )

  # short_idx <- el[,order(d)[1:floor(length(d)*0.1)]]
  # max_short_range <- el[,quantile(r2[short_idx], 0.95)]
  #
  # el[,r2_std := (r2 - b) / (max_short_range - b)]
  # el[,plot(d,r2)]
  #
  min_pos <- min(c(el$pos1, el$pos2))
  max_pos <- max(c(el$pos1, el$pos2))

  starts <- seq(min_pos, max_pos - window_size, by = step_size)
  ends   <- starts + window_size

  # i <- 1
  decay <- suppressWarnings(rbindlist(lapply(seq_along(starts), function(i) {

    sub <- el[
      pos1 >= starts[i] & pos1 < ends[i] &
        pos2 >= starts[i] & pos2 < ends[i]
    ]

    # d_small <- quantile(sub$d, 0.1)
    # max_short_range <- sub[,quantile(r2[d <= d_small], 0.95)]
    #
    # short_idx <- sub[,order(d)[1:floor(length(d)*0.1)]]
    # max_short_range <- sub[,quantile(r2[short_idx], 0.95)]
    #
    # sub[,r2_std := (r2 - b) / (max_short_range - b)]
    # sub[,plot(d,r2_std)]

    coefs <- tryCatch(
      coef_ld_dec(sub, q = q,  b = b,n_bins = target_dist_bins_for_decay),
      error = function(e) NULL
    )

    if (is.null(coefs)) {
      data.table(start = starts[i], end = ends[i])
    } else {
      data.table(start = starts[i], end = ends[i],coefs)
    }
  }),fill=TRUE,use.names = TRUE))

  decay[, k_star := thin$k_star]
  decay[, reduction := signif(reduction,2)]

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


##' Adaptive SNP Thinning
#'
#' Reduces SNP density to maintain approximately constant SNP counts
#' within LD-decay windows.
#'
#' @param pos SNP physical positions (bp).
#' @param W_LD Target LD window size (bp).
#' @param n_snps_for_decay Desired number of SNPs per window.
#'
#' @return A list containing:
#' \describe{
#'   \item{keep}{Indices of retained SNPs.}
#'   \item{k_star}{Effective sliding window size (in SNP count).}
#' }
#'
#' @export
adaptive_thinning <- function(pos,
                              W_LD,
                              n_snps_for_decay = 1000) {

  pos <- sort(pos)
  spacing <- diff(pos)
  median_spacing <- median(spacing)

  target_spacing <- W_LD / n_snps_for_decay

  if (median_spacing >= target_spacing) {

    keep <- seq_along(pos)
    eff_spacing <- median_spacing

  } else {

    thin_factor <- ceiling(target_spacing / median_spacing)
    keep <- seq(1, length(pos), by = thin_factor)
    eff_spacing <- median_spacing * thin_factor
  }

  k_star <- floor(W_LD / eff_spacing)

  list(
    keep = keep,
    k_star = k_star
  )
}

#' Fit LD Decay Model
#'
#' Fits the nonlinear decay model
#' \deqn{r^2(d) = b + \frac{c - b}{1 + a(d - d_0)}}
#' to binned LD quantiles.
#'
#' @param el_ld LD edge data table.
#' @param q Quantile for LD binning.
#' @param dist_unit Distance bin size.
#' @param b Background LD.
#'
#' @return Named vector with parameters \code{c}, \code{a}, \code{d0}.
#'
#' @export
coef_ld_dec <- function(dt,
                        q = 0.95,
                        b,
                        n_bins = 40,
                        min_dist_quantile = 0.02,
                        return_full=FALSE) {

  if (nrow(dt) < 100) return(NULL)

  d <- dt$d
  r2 <- dt$r2

  # Remove zero or extremely small distances
  # d_min <- quantile(d[d > 0], min_dist_quantile, na.rm = TRUE)
  # keep <- d > d_min
  #
  # d <- d[keep]
  # r2 <- r2[keep]

  if (length(d) < 100) return(NULL)

  # Log-distance bins (more stable near zero)
  d <- d[d > 0]
  log_d <- log(d)

  breaks <- seq(min(log_d), max(log_d), length.out = n_bins + 1)

  bin_id <- cut(log_d, breaks = breaks, include.lowest = TRUE)

  agg <- data.table::data.table(d = d, r2 = r2, bin = bin_id)[,
                                                              .(
                                                                d_mid = exp(mean(log(d))),
                                                                r2_q = quantile(r2, q, na.rm = TRUE),
                                                                n = .N
                                                              ),
                                                              by = bin
  ]

  agg <- agg[!is.na(r2_q) & n > 10]

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

  return(data.table(c=coefs["c"],a=coefs["a"],b,data=list(agg)))

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


parallel_apply <- function(X, FUN, cores = 1) {

  if (cores > 1 && .Platform$OS.type != "windows") {
    parallel::mclapply(X, FUN, mc.cores = cores)
  } else {
    lapply(X, FUN)
  }
}

build_hist_from_dt <- function(dt, dist_unit, r2_unit) {

  if (nrow(dt) == 0) return(NULL)

  # --- symmetry correction ---
  dt[, side := data.table::fifelse(pos_other > pos, "R", "L")]

  maxL <- if (any(dt$side == "L")) max(dt$d[dt$side == "L"]) else 0
  maxR <- if (any(dt$side == "R")) max(dt$d[dt$side == "R"]) else 0
  S <- min(maxL, maxR)

  if (S > 0 && maxL != maxR) {
    long_side <- if (maxL > maxR) "L" else "R"
    tail_dt <- dt[side == long_side & d > S]
    if (nrow(tail_dt) > 0)
      dt <- data.table::rbindlist(list(dt, tail_dt))
  }

  # --- integer bins ---
  dt[, dist_idx := floor(d / dist_unit)]
  dt[, r2_idx   := floor(r2 / r2_unit)]

  hist2d <- dt[, .N, by = .(dist_idx, r2_idx)]

  mat <- xtabs(N ~ dist_idx + r2_idx, data = hist2d)
  mat <- as.matrix(mat)

  # convert indices to physical scale
  dist_vals <- as.numeric(rownames(mat)) * dist_unit
  r2_vals   <- as.numeric(colnames(mat)) * r2_unit

  rownames(mat) <- signif(dist_vals, 6)
  colnames(mat) <- signif(r2_vals, 6)
  return(mat)

}

summarize_hist_shell <- function(hist_row,
                                 b = NULL,
                                 type = c("mean", "median", "excess")) {

  type <- match.arg(type)

  total <- sum(hist_row)
  if (total == 0) return(NA_real_)

  r2_vals <- as.numeric(names(hist_row))

  if (type == "median") {
    cs <- cumsum(hist_row)
    idx <- which(cs >= total / 2)[1]
    return(r2_vals[idx])
  }

  if (type == "mean") {
    return(sum(r2_vals * hist_row) / total)
  }

  if (type == "excess") {
    excess <- pmax(r2_vals - b, 0)
    return(sum(excess * hist_row) / total)
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

derive_ld_radius <- function(a, eps) {
  (1 / a) * (1/eps - 1)
}

#' Build Per-SNP LD Histograms
#'
#' Constructs two-dimensional histograms of LD (r²) as a function
#' of physical distance for each SNP within an adaptive LD radius.
#'
#' @param gds GDS file handle.
#' @param chr_idx SNP indices for chromosome.
#' @param snps_chr SNP IDs.
#' @param pos_chr Physical SNP positions.
#' @param d_star Maximum LD integration radius (bp).
#' @param r2_unit Bin width for r² values.
#' @param n_dist_target Target number of distance bins.
#' @param cores Number of CPU cores.
#'
#' @return A list of matrices, one per SNP. Each matrix contains
#'   counts of LD pairs binned by distance (rows) and r² (columns).
#'
#' @export
build_ld_histograms <- function(
    gds,
    chr_idx,
    snps_chr,
    pos_chr,
    d_star,
    r2_unit = 0.001,
    n_dist_target = 100,
    k_max = 1000,
    cores = 1,
    compression = c("hist", "median_only")
) {

  compression <- match.arg(compression)

  spacing <- diff(sort(pos_chr))
  q_spacing <- quantile(spacing, 0.75)
  k_star <- ceiling(d_star / q_spacing)

  el <- get_el(
    gds = gds,
    idx = chr_idx,
    slide_win_ld = min(k_max, k_star),
    cores = cores
  )

  el <- data.table::rbindlist(list(
    el[, .(SNP = SNP1, pos = pos1, pos_other = pos2, r2, d)],
    el[, .(SNP = SNP2, pos = pos2, pos_other = pos1, r2, d)]
  ))

  data.table::setkey(el, SNP)

  if (compression == "hist") {
    # i <- 1
    hist_list <- parallel_apply(seq_along(snps_chr), function(i) {

      dt <- el[J(snps_chr[i]), nomatch = 0]
      dt[, r2 := pmin(r2, 1)]

      max_d <- max(dt$d)
      dist_unit <- signif(max_d / n_dist_target, 2)

      build_hist_from_dt(dt, dist_unit, r2_unit)

    }, cores = cores)

    return(hist_list)
  }

  if (compression == "median_only") {

    shell_list <- parallel_apply(seq_along(snps_chr), function(i) {

      dt <- el[J(snps_chr[i]), nomatch = 0]
      dt[, r2 := pmin(r2, 1)]

      compute_shell_medians(
        dt = dt,
        d_star = d_star,
        n_dist_target = n_dist_target
      )

    }, cores = cores)

    return(shell_list)
  }
}

integrate_ld_kernel <- function(hist_mat,
                                a,
                                b,
                                d_star,
                                shell_type = "excess") {

  if (is.null(hist_mat)) return(NA_real_)

  dist_vals <- as.numeric(rownames(hist_mat))
  keep <- dist_vals <= d_star

  if (!any(keep)) return(NA_real_)

  h_sub <- hist_mat[keep, , drop = FALSE]
  dist_vals_sub <- as.numeric(rownames(h_sub))

  # Kernel weights
  w_d <- a / (1 + a * dist_vals_sub)^2

  if (length(dist_vals_sub) > 1) {
    delta_d <- c(diff(dist_vals_sub),
                 tail(diff(dist_vals_sub), 1))
  } else {
    delta_d <- 1
  }

  w_use <- w_d * delta_d
  if (sum(w_use) == 0) return(NA_real_)
  w_use <- w_use / sum(w_use)

  # Shell summaries
  shell_vals <- apply(h_sub, 1, summarize_hist_shell,
                      b = b,
                      type = shell_type)

  sum(shell_vals * w_use)
}

#' Compute LD-Based Summary Statistic
#'
#' Projects stored LD histograms into scalar LD summary statistics
#' using either kernel-weighted integration or hard truncation.
#'
#' @param ld_structure An object of class `"ld_structure"`.
#' @param method One of:
#'   \describe{
#'     \item{"ld_int"}{Kernel-weighted integration of shell summaries.}
#'     \item{"rho_w"}{Hard distance truncation equivalent to rho thresholding.}
#'     \item{"ld_int_median"}{Kernel-weighted integration using weighted median across shells.}
#'   }
#' @param eps Tail tolerance for LD integration (used for `"ld_int"`).
#' @param d_window Physical distance cutoff (used for `"rho_w"`).
#' @param shell_type Shell summary statistic:
#'   `"mean"`, `"median"`, or `"excess"`.
#' @param cores Number of CPU cores.
#'
#' @return Numeric vector of LD summary values (one per SNP).
#'
#' @export
compute_ld_summary <- function(
    ld_structure,
    method = c("ld_int", "rho_w", "ld_int_median"),
    eps = 0.005,
    d_window = NULL,
    shell_type = "median",
    cores = 1
) {

  method <- match.arg(method)

  results <- parallel_apply(ld_structure$by_chr, function(chr_obj) {

    hist_obj <- chr_obj$hist_obj
    a <- chr_obj$decay_sum$a
    b <- chr_obj$decay_sum$b

    d_star <- derive_ld_radius(a, eps)

    # -----------------------------
    # NEW: median_shell backend
    # -----------------------------
    if (is.list(hist_obj) &&
        length(hist_obj) > 0 &&
        all(vapply(hist_obj,
                   function(x) is.null(x) || is.data.frame(x),
                   logical(1)))) {

      return(sapply(hist_obj, function(shell_dt) {

        if (is.null(shell_dt)) return(NA_real_)

        integrate_ld_kernel_median_shell(
          shell_dt = shell_dt,
          a = a,
          d_star = d_star
        )
      }))
    }

    # -----------------------------
    # OLD: histogram backend
    # -----------------------------
    if (method == "ld_int" || method == "ld_int_median") {

      return(sapply(hist_obj, function(h)
        integrate_ld_kernel_median(
          hist_mat = h,
          a = a,
          b = b,
          d_star = d_star,
          shell_type = shell_type
        )
      ))
    }

    if (method == "rho_w") {

      if (is.null(d_window))
        stop("d_window must be supplied for rho_w.")

      return(sapply(hist_obj, function(h)
        compute_ld_rhow_from_hist(
          hist_mat = h,
          d_window = d_window,
          shell_type = shell_type,
          b = b
        )
      ))
    }

  }, cores = cores)

  unlist(results)
}
weighted_median <- function(x, w) {

  if (length(x) == 0) return(NA_real_)

  ord <- order(x)
  x <- x[ord]
  w <- w[ord]

  w <- w / sum(w)

  cw <- cumsum(w)

  x[which(cw >= 0.5)[1]]
}

integrate_ld_kernel_median <- function(hist_mat,
                                       a,
                                       b,
                                       d_star,
                                       shell_type = "median") {

  if (is.null(hist_mat)) return(NA_real_)

  dist_vals <- as.numeric(rownames(hist_mat))
  keep <- dist_vals <= d_star

  if (!any(keep)) return(NA_real_)

  h_sub <- hist_mat[keep, , drop = FALSE]
  dist_vals_sub <- as.numeric(rownames(h_sub))

  # decay weights
  w_d <- a / (1 + a * dist_vals_sub)^2

  if (length(dist_vals_sub) > 1) {
    delta_d <- c(diff(dist_vals_sub),
                 tail(diff(dist_vals_sub), 1))
  } else {
    delta_d <- 1
  }

  w_use <- w_d * delta_d
  if (sum(w_use) == 0) return(NA_real_)

  # shell summaries
  shell_vals <- apply(h_sub, 1, summarize_hist_shell,
                      b = b,
                      type = shell_type)

  weighted_median(shell_vals, w_use)
}

#d_window
compute_ld_rhow_from_hist <- function(hist_mat,
                                      d_window,
                                      shell_type = "median",
                                      b = NULL) {

  if (is.null(hist_mat)) return(NA_real_)

  dist_vals <- as.numeric(rownames(hist_mat))
  keep <- dist_vals <= d_window

  if (!any(keep)) return(NA_real_)

  h_sub <- hist_mat[keep, , drop = FALSE]

  # Collapse all rows into one distribution
  combined <- colSums(h_sub)

  summarize_hist_shell(combined,
                       b = b,
                       type = shell_type)
}


uniform_thin_to_distance <- function(pos,
                                     d_star,
                                     k_max,
                                     tol = 1e-8) {

  pos <- sort(pos)

  # required spacing
  target_spacing <- d_star / k_max

  # current median spacing
  med_spacing <- median(diff(pos))

  # if already sufficient, no thinning needed
  if (k_max * med_spacing >= d_star) {
    return(seq_along(pos))
  }

  keep <- integer()
  last_kept_pos <- -Inf

  for (i in seq_along(pos)) {
    if ((pos[i] - last_kept_pos) >= (target_spacing - tol)) {
      keep <- c(keep, i)
      last_kept_pos <- pos[i]
    }
  }

  keep
}


compute_shell_medians <- function(dt,
                                  d_star,
                                  n_dist_target) {

  if (nrow(dt) == 0) return(NULL)

  dt <- dt[d <= d_star]
  if (nrow(dt) == 0) return(NULL)

  dist_unit <- d_star / n_dist_target
  dt[, dist_idx := floor(d / dist_unit)]

  shell_dt <- dt[, .(
    median_r2 = median(r2)
  ), by = dist_idx]

  shell_dt
}

integrate_ld_kernel_median_shell <- function(shell_dt,
                                             a,
                                             d_star) {

  if (is.null(shell_dt) || nrow(shell_dt) == 0)
    return(NA_real_)

  dist_unit <- d_star / max(shell_dt$dist_idx)

  dist_vals <- shell_dt$dist_idx * dist_unit

  w_d <- a / (1 + a * dist_vals)^2

  if (length(dist_vals) > 1) {
    delta_d <- c(diff(dist_vals), tail(diff(dist_vals), 1))
  } else {
    delta_d <- 1
  }

  w_use <- w_d * delta_d
  if (sum(w_use) == 0) return(NA_real_)
  w_use <- w_use / sum(w_use)

  weighted_median(shell_dt$median_r2, w_use)
}

compute_ld_int_kernel_smooth <- function(dt, a, d_star, n_bins = 20) {

  dt <- dt[d <= d_star]
  if (nrow(dt) == 0) return(NA_real_)

  # sort by distance
  dt <- dt[order(d)]

  # equal-width bins
  breaks <- seq(0, d_star, length.out = n_bins + 1)
  dt[, bin := cut(d, breaks = breaks, include.lowest = TRUE, labels = FALSE)]

  shell_dt <- dt[, .(
    median_r2 = median(r2)
  ), by = bin]

  # representative distance per bin
  dist_vals <- breaks[shell_dt$bin]

  # kernel weights
  w_d <- a / (1 + a * dist_vals)^2

  delta_d <- diff(breaks)
  w_use <- w_d * delta_d[shell_dt$bin]
  w_use <- w_use / sum(w_use)

  weighted_median(shell_dt$median_r2, w_use)
}

compute_ldw_from_hist <- function(hist_mat, d_window) {

  if (is.null(hist_mat)) return(NA_real_)

  dist_vals <- as.numeric(rownames(hist_mat))
  keep <- dist_vals <= d_window

  if (!any(keep)) return(NA_real_)

  h_sub <- hist_mat[keep, , drop = FALSE]

  # Collapse into 1D r² distribution
  r2_counts <- colSums(h_sub)

  total <- sum(r2_counts)
  if (total == 0) return(NA_real_)

  r2_vals <- as.numeric(colnames(hist_mat))

  # weighted median
  cs <- cumsum(r2_counts)
  r2_vals[which(cs >= total/2)[1]]
}
