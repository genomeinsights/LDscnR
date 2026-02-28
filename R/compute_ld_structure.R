#' Compute LD Structure and Decay Parameters
#'
#' Computes chromosome-wise linkage disequilibrium (LD) structure,
#' including LD-decay estimation and filtered edge lists for LD-weighted
#' statistics.
#'
#' The function performs the following steps:
#' \enumerate{
#'   \item Estimates background LD from genome-wide or inter-chromosomal pairs.
#'   \item Fits LD-decay curves in sliding windows per chromosome.
#'   \item Computes both median and robust decay summaries.
#'   \item Derives a maximum LD distance (\eqn{W_{\max}}) from the decay model.
#'   \item Constructs reduced LD edge lists for downstream LD-weighted analysis.
#' }
#'
#' LD decay is modeled as:
#' \deqn{f(d) = b + \frac{c - b}{1 + a d}}
#'
#' where:
#' \itemize{
#'   \item \eqn{b} is background LD,
#'   \item \eqn{a} is the decay rate,
#'   \item \eqn{c} is short-distance LD,
#'   \item \eqn{d} is physical distance (bp).
#' }
#'
#' @param gds GDS file handle containing genotype data.
#' @param q Quantile used for LD-decay fitting (default 0.95).
#' @param n_sub Number of SNPs used for background LD estimation.
#' @param window_size Physical window size (bp) for decay fitting.
#' @param step_size Step size (bp) for sliding decay windows.
#' @param cores Number of CPU cores for LD computation.
#' @param dist_unit Distance bin size (bp) used for decay fitting.
#' @param K_target Target number of pairwise LD comparisons per SNP.
#' @param n_win Target number of SNPs per decay window.
#' @param prob_robust Quantile threshold for excluding slowest-decaying
#'   windows in robust estimation.
#'
#' @return An object of class `"ld_structure"` containing per-chromosome:
#' \describe{
#'   \item{edges}{Filtered LD edge lists for LD-weighted statistics.}
#'   \item{decay}{Window-wise decay fits.}
#'   \item{summary}{Median and robust decay summaries.}
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
    target_dist_bins_for_decay = 40,
    n_snps_for_decay = 500,
    ## for histograms
    n_dist_target_for_hist = 100,
    eps = 0.005,
    r2_unit = 0.001,
    ## cores
    cores = 1
) {

  ids  <- .read_gds_ids(gds)
  chrs <- unique(ids$snp_chr)

  b <- estimate_background_ld(gds, n_sub = n_sub_bg, q = q)

  out_by_chr <- vector("list", length(chrs))
  names(out_by_chr) <- chrs

  for (ch in chrs) {

    message("Processing ", ch)

    chr_idx  <- which(ids$snp_chr == ch)
    pos_chr  <- ids$snp_pos[chr_idx]
    snps_chr <- ids$snp_id[chr_idx]


    window_size <- (max(pos_chr) - min(pos_chr)) / n_win_decay
    step_size   <- window_size * overlap


    # ---- decay ----
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

    decay_sum <- summarize_decay(decay, prob_robust)
    decay_sum[,k_star:=decay[1,k_star]]

    if (is.null(decay_sum)) next

    decay_sum[, Chr := ch]

    # --- Layer 2: Histograms ---

    d_star <- derive_ld_radius(decay_sum$a, eps)

    histograms <- build_ld_histograms(
      gds=gds,
      chr_idx   = chr_idx,
      snps_chr  = snps_chr,
      pos_chr   = pos_chr,
      d_star    = d_star,
      r2_unit   = r2_unit,
      n_dist_target = n_dist_target_for_hist,
      cores = cores
    )

    # ---- collect data ----
    out_by_chr[[ch]] <- list(
      snp_ids    = snps_chr,
      hist_obj   = histograms,
      decay      = decay,
      decay_sum  = decay_sum
    )
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


# target_snps_per_window = 500
# mean_spacing <- median(diff(sort(snp_pos)))
# window_size  <- target_snps_per_window * mean_spacing
# step_size    <- window_size * 0.5

#' Summarize LD Decay Parameters
#'
#' Computes median and robust summaries of LD-decay parameters
#' across sliding windows.
#'
#' Robust estimates exclude the slowest-decaying windows based on
#' a quantile threshold of parameter \eqn{a}.
#'
#' @param decay_dt Data table containing window-wise decay parameters.
#' @param prob_robust Quantile threshold for robust filtering
#'   (default 0.5 keeps upper 50% of \eqn{a} values).
#'
#' @return A list with elements:
#' \describe{
#'   \item{median}{Median decay parameters.}
#'   \item{robust}{Robust decay parameters.}
#' }
#'
#' @keywords internal
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
#' Decay is estimated using nonlinear least squares applied to binned
#' quantiles of \eqn{r^2}.
#'
#' @param gds GDS file handle.
#' @param chr_idx SNP indices for the chromosome.
#' @param snp_pos SNP physical positions (bp).
#' @param b Background LD.
#' @param window_size Window size (bp).
#' @param step_size Step size (bp).
#' @param q Quantile for LD summary (default 0.95).
#' @param dist_unit Distance bin size (bp).
#' @param cores Number of CPU cores.
#' @param n_win Target number of SNPs per window.
#'
#' @return A \code{data.table} with window-wise decay parameters.
#'
#' @keywords internal
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
                               n_snps_for_decay = 500) {




  thin <- adaptive_thinning(snp_pos, W_LD = window_size, n_snps_for_decay = n_snps_for_decay)
  reduction <- 1-length(thin$keep)/length(chr_idx)

  el <- get_el(
    gds,
    idx = chr_idx[thin$keep],
    slide_win_ld = thin$k_star,
    cores = cores
  )

  min_pos <- min(c(el$pos1, el$pos2))
  max_pos <- max(c(el$pos1, el$pos2))

  starts <- seq(min_pos, max_pos - window_size, by = step_size)
  ends   <- starts + window_size
  # i <- 1
  decay <- rbindlist(lapply(seq_along(starts), function(i) {

    sub <- el[
      pos1 >= starts[i] & pos1 < ends[i] &
        pos2 >= starts[i] & pos2 < ends[i]
    ]

    coefs <- tryCatch(
      coef_ld_dec(sub, q = q,  b = b,n_bins = target_dist_bins_for_decay),
      error = function(e) NULL
    )

    if (is.null(coefs)) {
      data.table(start = starts[i], end = ends[i])
    } else {
      data.table(start = starts[i], end = ends[i],coefs)
    }
  }),fill=TRUE,use.names = TRUE)

  decay[, k_star := k_star]
  decay[, reduction := reduction]

  return(decay)
}

#' Estimate LD Decay for a Chromosome
#'
#' Fits LD-decay curves in sliding genomic windows for a single chromosome.
#'
#' Decay is estimated using nonlinear least squares applied to binned
#' quantiles of \eqn{r^2}.
#'
#' @param gds GDS file handle.
#' @param chr_idx SNP indices for the chromosome.
#' @param snp_pos SNP physical positions (bp).
#' @param b Background LD.
#' @param window_size Window size (bp).
#' @param step_size Step size (bp).
#' @param q Quantile for LD summary (default 0.95).
#' @param dist_unit Distance bin size (bp).
#' @param cores Number of CPU cores.
#' @param n_win Target number of SNPs per window.
#'
#' @return A \code{data.table} with window-wise decay parameters.
#'
#' @keywords internal
estimate_background_ld <- function(gds,
                                   n_sub = 5000,
                                   q = 0.95) {

  b <- get_bg_ld(gds, idx = NULL, n_sub = n_sub, q = q)

  message("Background LD: ", round(b, 3))
  b
}


#' Adaptive SNP Thinning for Decay Estimation
#'
#' Reduces SNP density to maintain approximately constant SNP counts
#' within LD-decay windows.
#'
#' @param pos SNP physical positions.
#' @param W_LD LD window size.
#' @param n_snps_for_decay Target SNP count per window.
#'
#' @return List with retained SNP indices and slide window size.
#'
#' @keywords internal
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
#' @keywords internal
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
  d_min <- quantile(d[d > 0], min_dist_quantile, na.rm = TRUE)
  keep <- d > d_min

  d <- d[keep]
  r2 <- r2[keep]


  if (length(d) < 100) return(NULL)

  # Log-distance bins (more stable near zero)
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

  #agg[,plot(d_mid,r2_q)]

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
      weights = n,
      control = nls.control(warnOnly = TRUE)
    ),
    error = function(e) NULL
  )


  if (is.null(fit)) return(NULL)
  coefs <- coef(fit)
  #coefs
  #agg[,points(d_mid,r2_q,col="red")]
  #agg[,b:=b]
  #agg[order(d_mid),lines(d_mid, b + (coefs["c"]  - b) / (1 + coefs["a"]  * d_mid),col="red")]

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
#' @keywords internal
#'
d_from_rho <- function(a, rho){
  rho / (a * (1 - rho))
}

ld_from_rho <- function(b, c = 1, rho){
  b + (c - b) * (1 - rho)
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
#' @keywords internal
get_bg_ld <- function(gds,
                      idx = NULL,
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

  b <- stats::quantile(r2_inter, probs = q, na.rm = TRUE)

  return(b)
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
  #x <- out
  if (!inherits(x, "ld_structure"))
    stop("Object must be of class 'ld_structure'.")

  chrs <- names(x$histograms)
  n_chr <- length(chrs)

  n_snps <-length(x$SNP_ids)

  cat("LD Structure Object\n")
  cat("-------------------\n")
  cat("Chromosomes:", n_chr, "\n")
  cat("Total SNPs:", n_snps, "\n\n")

  if (!is.null(x$decay_summary) && nrow(x$decay_summary) > 0) {

    cat("Robust LD-decay parameters (median across chromosomes):\n")

    cat("  a  =", signif(median(x$decay_summary$a, na.rm = TRUE), 3), "\n")
    cat("  c  =", signif(median(x$decay_summary$c, na.rm = TRUE), 3), "\n")

    if ("b" %in% names(x$decay_summary))
      cat("  b  =", signif(median(x$decay_summary$b, na.rm = TRUE), 3), "\n")
  }

  cat("\nParameters used:\n")

  for (nm in names(x$params)) {
    cat(" ", nm, "=", x$params[[nm]], "\n")
  }

  cat("\nHistograms stored (not printed).\n")

  invisible(x)
}


#' Summary Method for ld_structure Objects
#'
#' Returns chromosome-wise summary statistics of LD structure.
#'
#' @param object An object of class `"ld_structure"`.
#' @param ... Additional arguments (unused).
#'
#' @return A \code{data.table} summarizing decay parameters and edge counts.
#'
#' @export
summary.ld_structure <- function(object, ...) {

  if (!inherits(object, "ld_structure"))
    stop("Object must be of class 'ld_structure'.")
  object$summary

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

ld_exc_from_hist <- function(cum_mat, row_index, b, r2_unit) {

  if (is.null(cum_mat) || !(row_index %in% rownames(cum_mat)))
    return(NA_real_)

  counts <- cum_mat[row_index, , drop = TRUE]
  total <- sum(counts)

  if (total == 0) return(NA_real_)

  r2_idx <- as.numeric(colnames(cum_mat))
  r2_vals <- r2_idx * r2_unit

  excess <- pmax(r2_vals - b, 0)

  sum(excess * counts) / total
}

median_from_hist <- function(cum_mat, row_index) {

  row_counts <- cum_mat[row_index, ]
  total <- sum(row_counts)

  if (total == 0) return(NA_real_)

  cs <- cumsum(row_counts)
  idx <- which(cs >= total / 2)[1]

  as.numeric(colnames(cum_mat)[idx])
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

build_ld_histograms <- function(
    gds,
    chr_idx,
    snps_chr,
    pos_chr,
    d_star,
    r2_unit,
    n_dist_target,
    cores = 1
) {

  library(data.table)

  # Convert LD radius to SNP window
  spacing <- diff(sort(pos_chr))
  q_spacing <- quantile(spacing, 0.75)
  k_star <- ceiling(d_star / q_spacing)

  el <- get_el(
    gds = gds,
    idx = chr_idx,
    slide_win_ld = k_star,
    cores = cores
  )

  el <- rbindlist(list(
    el[, .(SNP = SNP1, pos = pos1, pos_other = pos2, r2, d)],
    el[, .(SNP = SNP2, pos = pos2, pos_other = pos1, r2, d)]
  ))

  setkey(el, SNP)

  hist_list <- vector("list", length(snps_chr))

  for (i in seq_along(snps_chr)) {

    dt <- el[J(snps_chr[i]), nomatch = 0]

    if (nrow(dt) == 0) {
      hist_list[[i]] <- NULL
      next
    }

    dt[, r2 := pmin(r2, 1)]

    max_d <- max(dt$d)
    dist_unit <- signif(max_d / n_dist_target, 2)

    hist_list[[i]] <- build_hist_from_dt(dt, dist_unit, r2_unit)
  }

  return(hist_list)

}

integrate_ld_kernel <- function(
    hist_mat,
    a,
    b,
    d_star
) {

  if (is.null(hist_mat)) return(NA_real_)

  dist_vals <- as.numeric(rownames(hist_mat))
  keep <- dist_vals <= d_star

  if (!any(keep)) return(NA_real_)

  h_sub <- hist_mat[keep, , drop = FALSE]
  dist_vals_sub <- as.numeric(rownames(h_sub))

  # kernel weights
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

  # excess LD
  r2_vals <- as.numeric(colnames(h_sub))
  excess  <- pmax(r2_vals - b, 0)

  E_shell <- as.numeric(h_sub %*% excess)
  N_shell <- rowSums(h_sub)

  LD_shell <- E_shell / pmax(N_shell, 1)

  sum(LD_shell * w_use)
}


compute_ld_integral <- function(
    hist_obj,
    eps = 0.001
) {

  hist_list <- hist_obj$hist_list
  k_star    <- hist_obj$k_star
  a         <- hist_obj$decay_sum$a
  b         <- hist_obj$decay_sum$b

  # derive LD radius dynamically
  d_star <- derive_ld_radius(a, eps)

  LD_int <- sapply(hist_list, function(h)
    integrate_ld_kernel(
      hist_mat = h,
      a = a,
      b = b,
      d_star = d_star
    )
  )

  list(
    LD_int = LD_int,
    d_star = d_star,
    k_star = k_star
  )
}

summarize_hist_shell <- function(hist_row, b = NULL) {

  if (is.null(hist_row))
    return(list(mean_excess = NA_real_, median = NA_real_))

  total <- sum(hist_row)
  if (total == 0)
    return(list(mean_excess = NA_real_, median = NA_real_))

  r2_vals <- as.numeric(names(hist_row))

  # median
  cs <- cumsum(hist_row)
  idx <- which(cs >= total / 2)[1]
  median_r2 <- r2_vals[idx]

  # excess mean
  if (!is.null(b)) {
    excess <- pmax(r2_vals - b, 0)
    mean_excess <- sum(excess * hist_row) / total
  } else {
    mean_excess <- NA_real_
  }

  list(
    mean_excess = mean_excess,
    median = median_r2
  )
}

# plot(unlist(decay[,sapply(data,nrow)]),as.vector(na.omit(decay$a)))
#
# par(mfcol=c(4,5))
# for(i in which(!is.na(decay$a))){
#   decay[i,plot(data[[1]]$d_mid,data[[1]]$r2_q)]
#   a = decay[i,a]
#   c = decay[i,c]
#   decay[i,data[[1]]][order(d_mid),lines(d_mid, b + (c  - b) / (1 + a  * d_mid),col="red")]
# }
#
# par(mfcol=c(1,1))
# a = decay[1,a]
# c = decay[1,c]
# decay[1,data[[1]]][order(d_mid),plot(d_mid, b + (c  - b) / (1 + a  * d_mid),type="l",col="red",ylim=c(0,1))]
# for(i in which(!is.na(decay$a))[-1]){
#   decay[i,data[[1]]][order(d_mid),lines(d_mid, b + (c  - b) / (1 + a  * d_mid),type="l",ylim=c(0,1),col="red")]
# }
#
