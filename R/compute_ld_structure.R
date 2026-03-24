#' Compute Chromosome-wise LD Structure and Decay
#'
#' Estimates linkage disequilibrium (LD) structure and LD-decay parameters
#' across chromosomes from genotype data stored in a GDS file. The function
#' combines background LD estimation, chromosome-wise decay fitting, and
#' genome-wide summarization into a unified workflow.
#'
#' The procedure consists of:
#' \enumerate{
#'   \item Estimation of background LD from inter-chromosomal SNP pairs.
#'   \item Sliding-window estimation of LD-decay parameters per chromosome.
#'   \item Robust aggregation of decay parameters across windows.
#'   \item Prediction of decay parameters based on chromosome size.
#'   \item Derivation of recommended LD window sizes for user-defined LD thresholds (\eqn{\rho}).
#' }
#'
#' LD decay is modeled as:
#' \deqn{r^2(d) = b + \frac{c - b}{1 + a d}}
#'
#' where:
#' \itemize{
#'   \item \eqn{a} controls the rate of decay,
#'   \item \eqn{b} is background LD (long-distance baseline),
#'   \item \eqn{c} is short-range LD,
#'   \item \eqn{d} is physical distance (bp).
#' }
#'
#' @param gds A GDS file handle containing genotype data.
#' @param el_data_folder Optional path to store LD edge lists (currently not used internally).
#' @param q Quantile of \eqn{r^2} used for decay fitting (default = 0.95).
#' @param n_sub_bg Number of SNPs used to estimate background LD.
#' @param n_win_decay Number of sliding windows per chromosome.
#' @param overlap Proportion of overlap between consecutive windows (0–1).
#' @param max_SNPs_decay Maximum number of SNPs sampled per chromosome for decay estimation.
#' @param prob_robust Central proportion of windows retained for robust summarization.
#' @param max_pairs Maximum number of SNP pairs per window used in decay fitting.
#' @param n_strata Number of distance strata used when subsampling SNP pairs.
#' @param keep_el Logical; whether to store full LD edge lists per chromosome.
#' @param slide Sliding window size in number of SNPs used for LD estimation.
#' @param rho_targets Numeric vector of target LD thresholds used to derive recommended window sizes.
#' @param cores Number of CPU cores for parallel computation.
#'
#' @return An object of class \code{"ld_decay"} containing:
#' \describe{
#'   \item{by_chr}{List of per-chromosome results including decay fits and optional LD edge lists.}
#'   \item{decay_sum}{Data table of chromosome-wise decay parameters and derived quantities.}
#'   \item{decay_model}{Robust regression model linking decay rate \eqn{a} to chromosome size.}
#'   \item{recommendation}{Suggested LD window sizes for specified \eqn{\rho} thresholds.}
#'   \item{params}{List of parameters used in the computation.}
#' }
#'
#' @details
#' The function estimates LD decay independently per chromosome and then
#' fits a genome-wide model relating decay rate (\eqn{a}) to chromosome size.
#' This allows extrapolation of LD behavior and derivation of consistent
#' window sizes across heterogeneous genomic architectures.
#'
#' Recommended window sizes are provided in SNP units to match downstream
#' functions that operate on marker indices rather than physical distance.
#'
#' @export
compute_LD_decay <- function(
    gds,
    el_data_folder = NULL,
    q = 0.95,
    n_sub_bg = 5000,
    n_win_decay = 20,
    overlap = 0.5,
    max_SNPs_decay = Inf,
    prob_robust = 0.95,
    max_pairs = 5000,
    n_strata = 20,
    keep_el = FALSE,
    slide = 1000,
    rho_targets = c(0.90, 0.95, 0.99),
    cores = 1
) {

  if (!is.null(el_data_folder)) {
    if (!dir.exists(el_data_folder)) dir.create(el_data_folder, recursive = TRUE)
  }

  ids  <- .read_gds_ids(gds)
  chrs <- unique(ids$snp_chr)

  b <- estimate_background_ld(gds, n_sub = n_sub_bg, q = q)

  out_by_chr <- vector("list", length(chrs))
  names(out_by_chr) <- chrs


  message("Estimating LD-decay")
  pb <- txtProgressBar(min = 0, max = length(chrs)-1, style = 3)
  for (ch in chrs) {

    chr_idx  <- which(ids$snp_chr == ch)
    pos_chr  <- ids$snp_pos[chr_idx]
    snps_chr <- ids$snp_id[chr_idx]




    setTxtProgressBar(pb, which(chrs==ch))


    ## robust chromosome-specific bp per SNP scale
    chr_size_bp <- max(pos_chr) - min(pos_chr)
    n_snps_chr  <- length(chr_idx)

    bp_per_snp <- if (n_snps_chr > 1) {
      chr_size_bp / (n_snps_chr - 1)
    } else {
      NA_real_
    }

    ## current slide converted to bp
    slide_bp <- slide * bp_per_snp

    sample_idx <- sort(sample(
      chr_idx,
      size = min(max_SNPs_decay, length(chr_idx)),
      replace = FALSE
    ))

    el <- get_el(
      gds = gds,
      idx = sample_idx,
      slide_win_ld = slide,   # still in SNPs for get_el
      cores = cores,
      by_chr = TRUE,
      symmetric = TRUE,
      edge_symmetry = FALSE
    )

    data.table::setkey(el, SNP)

    window_size <- chr_size_bp / n_win_decay
    step_size   <- window_size * overlap

    ## ---- LD-decay
    decay <- suppressWarnings(
      estimate_decay_chr(
        el = el[d < window_size],
        b = b,
        window_size = window_size,
        step_size   = step_size,
        max_pairs   = max_pairs,
        n_strata    = n_strata,
        q = q,
        cores = cores
      )
    )

    decay[, contrast := c - b]

    decay[, regime := ifelse(
      contrast < 0.05,
      "weak",
      "structured"
    )]

    decay_sum_chr <- summarize_decay(decay[regime == "structured", ])

    if (is.null(decay_sum_chr)) next

    decay_sum_chr[, chr_size := max(pos_chr)]
    decay_sum_chr[, n_snp_chr := n_snps_chr]
    decay_sum_chr[, bp_per_snp := bp_per_snp]
    decay_sum_chr[, slide_snp := slide]
    decay_sum_chr[, slide_bp := slide_bp]
    decay_sum_chr[, n_w_used := sum(na.omit(decay$regime == "structured"))]
    decay_sum_chr[, Chr := ch]

    ## rho actually covered by the user-supplied slide, using raw chromosome a
    if (is.finite(decay_sum_chr$a) && decay_sum_chr$a > 0 &&
        is.finite(slide_bp) && slide_bp > 0) {
      decay_sum_chr[, rho_slide_raw := (a * slide_bp) / (1 + a * slide_bp)]
    } else {
      decay_sum_chr[, rho_slide_raw := NA_real_]
    }

    out_by_chr[[ch]] <- list(
      snp_ids   = snps_chr,
      el        = if (keep_el) el else NULL,
      decay     = decay,
      decay_sum = decay_sum_chr
    )
  }
  close(pb)

  message("Predicting a from chromosome size")
  decay_sum <- data.table::rbindlist(
    lapply(out_by_chr, function(x) x$decay_sum),
    fill = TRUE
  )


    ## robust background decay model
    if(length(unique(decay_sum$chr_size))<3){
      message("Not enough unique chromosome sizes to predict, using median")
      d_mod <- NA
      decay_sum[, a_pred := median(a, na.rm = TRUE)]
      decay_sum[, c_pred := median(c, na.rm = TRUE)]
    }else{
      d_mod <- MASS::rlm(log(a) ~ log(chr_size), data = decay_sum)

      decay_sum[, a_pred := exp(predict(d_mod, decay_sum))]
      decay_sum[, c_pred := median(c, na.rm = TRUE)]

    }


  ## rho covered by the current slide using predicted background a
  decay_sum[, rho_slide_pred := (a_pred * slide_bp) / (1 + a_pred * slide_bp)]

  ## recommended slide for each requested rho target
  rho_targets <- sort(unique(rho_targets))
  rho_targets <- rho_targets[is.finite(rho_targets) & rho_targets > 0 & rho_targets < 1]

  for (rho_t in rho_targets) {
    nm_bp  <- paste0("slide_bp_rho_",  formatC(rho_t, format = "f", digits = 2))
    nm_snp <- paste0("slide_snp_rho_", formatC(rho_t, format = "f", digits = 2))

    ## required physical window
    decay_sum[, (nm_bp) := rho_t / (a_pred * (1 - rho_t))]

    ## convert required bp window back to SNP window
    decay_sum[, (nm_snp) := get(nm_bp) / bp_per_snp]
  }

  ## optional compact summary for user-facing reporting
  rec_cols_snp <- grep("^slide_snp_rho_", names(decay_sum), value = TRUE)

  recommendation <- list(
    slide_input_snp = slide,
    rho_targets = rho_targets,

    rho_covered = decay_sum[, .(
      Chr,
      chr_size,
      n_snp_chr,
      bp_per_snp,
      a_pred,
      slide_snp,
      slide_bp,
      rho_slide_pred
    )],

    suggested_slide_by_chr = decay_sum[, c(
      "Chr", "chr_size", "n_snp_chr", "bp_per_snp", "a_pred", rec_cols_snp
    ), with = FALSE],

    suggested_slide_summary = if (length(rec_cols_snp) > 0) {
      tmp <- lapply(rec_cols_snp, function(cc) {
        vals <- decay_sum[[cc]]
        data.table::data.table(
          target = sub("^slide_snp_rho_", "", cc),
          median = stats::median(vals, na.rm = TRUE),
          p90    = stats::quantile(vals, probs = 0.90, na.rm = TRUE),
          max    = max(vals, na.rm = TRUE)
        )
      })
      data.table::rbindlist(tmp)
    } else {
      NULL
    }
  )

  out <- list(
    by_chr         = out_by_chr,
    decay_sum      = decay_sum,
    decay_model    = d_mod,
    recommendation = recommendation,
    params         = list(
      q = q,
      n_sub_bg = n_sub_bg,
      n_win_decay = n_win_decay,
      overlap = overlap,
      el_data_folder = el_data_folder,
      n_strata = n_strata,
      max_pairs = max_pairs,
      keep_el = keep_el,
      prob_robust = prob_robust,
      slide = slide,             # in SNPs
      rho_targets = rho_targets,
      cores = cores
    )
  )

  class(out) <- "ld_decay"

  if (!is.null(recommendation$suggested_slide_summary)) {
    message("Current slide covers the following rho range (using predicted background a):")
    print(decay_sum[, .(
      min_rho = min(rho_slide_pred, na.rm = TRUE),
      median_rho = median(rho_slide_pred, na.rm = TRUE),
      max_rho = max(rho_slide_pred, na.rm = TRUE)
    )])

    message("Suggested slide windows in SNPs for target rho:")
    print(recommendation$suggested_slide_summary)
  }

  out
}


#' Summarize LD Decay Parameters Across Windows
#'
#' Computes robust summary statistics of LD-decay parameters estimated
#' across sliding windows within a chromosome.
#'
#' Extreme windows are removed using symmetric quantile trimming of the
#' decay parameter \eqn{a}, and median values are returned.
#'
#' @param decay_dt Data table containing window-wise decay parameters
#'   (columns \code{a}, \code{b}, \code{c}).
#' @param prob_robust Central proportion of windows retained (default = 0.95).
#'
#' @return A \code{data.table} with median estimates of:
#' \describe{
#'   \item{a}{Decay rate.}
#'   \item{b}{Background LD.}
#'   \item{c}{Short-range LD.}
#' }
#'
#' @details
#' Trimming is based on the distribution of \eqn{a}, which is typically
#' the most variable parameter across genomic windows.
#'
#' Returns \code{NULL} if insufficient valid windows are available.
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

#' Estimate LD Decay Within a Chromosome
#'
#' Fits LD-decay models in sliding genomic windows for a single chromosome
#' using pairwise LD data (edge list format).
#'
#' For each window:
#' \enumerate{
#'   \item SNP pairs are restricted to the window.
#'   \item Pairs are subsampled across distance strata.
#'   \item LD decay is fitted using nonlinear regression.
#' }
#'
#' @param el Data table of SNP pairs with LD values and positions.
#' @param b Background LD.
#' @param window_size Window size in base pairs.
#' @param step_size Step size between windows in base pairs.
#' @param q Quantile of \eqn{r^2} used for fitting.
#' @param max_pairs Maximum number of SNP pairs per window.
#' @param n_strata Number of distance strata for subsampling.
#' @param cores Number of CPU cores.
#'
#' @return A \code{data.table} with window coordinates and decay parameters.
#'
#' @details
#' Windows with insufficient data or failed fits return partial rows with
#' missing parameter estimates.
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

#' Estimate Background LD from Inter-chromosomal SNP Pairs
#'
#' Computes background LD as a high quantile of \eqn{r^2} between SNPs
#' located on different chromosomes.
#'
#' @param gds GDS file handle.
#' @param idx Optional vector of SNP indices.
#' @param n_sub Number of SNPs sampled for estimation.
#' @param q Quantile used to define background LD.
#'
#' @return Numeric scalar representing background LD (\eqn{b}).
#'
#' @details
#' Inter-chromosomal LD provides an empirical estimate of baseline
#' correlation unrelated to physical linkage. This value anchors the
#' LD-decay model and stabilizes parameter estimation.
#'
#' Sampling is performed proportionally across chromosomes.
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



#' Fit LD Decay Model to Stratified LD Data
#'
#' Fits a nonlinear LD-decay model to binned LD values derived from
#' distance-stratified SNP pairs.
#'
#' The model is:
#' \deqn{r^2(d) = b + \frac{c - b}{1 + a d}}
#'
#' @param dt_strata Data table of SNP pairs with distance strata.
#' @param q Quantile used to summarize LD within strata.
#' @param b Background LD.
#'
#' @return A \code{data.table} containing:
#' \describe{
#'   \item{a}{Decay rate.}
#'   \item{c}{Short-range LD.}
#'   \item{b}{Background LD.}
#'   \item{agg}{Binned data used for fitting.}
#'   \item{raw}{Original stratified data.}
#' }
#'
#' @details
#' The model is fitted using weighted nonlinear least squares,
#' where weights correspond to the number of SNP pairs per stratum.
#'
#' Returns \code{NULL} if fitting fails or insufficient data are available.
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



#' Convert Relative LD Threshold to Physical Distance
#'
#' Computes the physical distance corresponding to a relative LD threshold \eqn{\rho}.
#'
#' @param a Decay rate.
#' @param rho Relative LD threshold (0 < rho < 1).
#'
#' @return Distance in base pairs.
#'
#' @export
d_from_rho <- function(a, rho){
  rho / (a * (1 - rho))
}

#' Convert Relative LD Threshold to Expected LD Value
#'
#' Computes the expected \eqn{r^2} corresponding to a relative LD threshold \eqn{\rho}.
#'
#' @param b Background LD.
#' @param c Short-range LD (default = 1).
#' @param rho Relative LD threshold.
#'
#' @return Expected \eqn{r^2}.
#'
#' @export
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

#' Compute Local LD Support (ld_w)
#'
#' Computes per-SNP local LD support within a distance defined by a
#' relative LD threshold \eqn{\rho}.
#'
#' For each SNP, LD support is defined as the median \eqn{r^2} with
#' neighboring SNPs within the corresponding distance window.
#'
#' @param ld_decay Object of class \code{"ld_decay"}.
#' @param rho Relative LD threshold used to define the window.
#' @param cores Number of CPU cores.
#'
#' @return Numeric vector of LD support values (one per SNP).
#'
#' @details
#' The physical window is derived using:
#' \deqn{d = \frac{\rho}{a(1 - \rho)}}
#'
#' where \eqn{a} is the chromosome-specific decay rate.
#'
#' Requires \code{keep_el = TRUE} when calling \code{compute_LD_decay()}.
#'
#' @export
print.ld_decay <- function(x, digits = 3) {

  cat("<ld_decay>\n")

  ## ---- parameters
  if (!is.null(x$params)) {

    cat("\nRun parameters:\n")

    cat(
      "  slide window:",
      format(x$params$slide, big.mark = ","),
      "SNPs\n"
    )

    cat(
      "  background LD quantile q:",
      signif(x$params$q, digits),
      "\n"
    )

    cat(
      "  chromosomes analysed:",
      length(x$by_chr),
      "\n"
    )

    cat(
      "  keep_el:",
      x$params$keep_el,
      "\n"
    )
  }

  ds <- x$decay_sum

  if (!is.null(ds) && nrow(ds) > 0) {

    cat("\nDecay parameter summary:\n")

    if ("a_pred" %in% names(ds)) {

      cat(
        "  predicted a:",
        "median =", signif(stats::median(ds$a_pred, na.rm = TRUE), digits),
        " range = [",
        signif(min(ds$a_pred, na.rm = TRUE), digits), ", ",
        signif(max(ds$a_pred, na.rm = TRUE), digits), "]\n",
        sep = ""
      )
    }

    if ("c_pred" %in% names(ds)) {

      cat(
        "  predicted c:",
        signif(stats::median(ds$c_pred, na.rm = TRUE), digits),
        "\n"
      )
    }

    if ("n_w_used" %in% names(ds)) {

      cat(
        "  structured decay windows:",
        "median =", stats::median(ds$n_w_used, na.rm = TRUE),
        " range = [",
        min(ds$n_w_used, na.rm = TRUE), ", ",
        max(ds$n_w_used, na.rm = TRUE), "]\n",
        sep = ""
      )
    }
  }

  ## ---- rho coverage
  rec <- x$recommendation

  if (!is.null(rec)) {

    rho_cov <- rec$rho_covered

    if (!is.null(rho_cov) && "rho_slide_pred" %in% names(rho_cov)) {

      cat("\nCurrent slide window coverage:\n")

      cat(
        "  rho covered:",
        "median =", signif(stats::median(rho_cov$rho_slide_pred, na.rm = TRUE), digits),
        " range = [",
        signif(min(rho_cov$rho_slide_pred, na.rm = TRUE), digits), ", ",
        signif(max(rho_cov$rho_slide_pred, na.rm = TRUE), digits), "]\n",
        sep = ""
      )

      cat(
        "  (slide =", format(x$params$slide, big.mark = ","), "SNPs)\n"
      )
    }

    ## ---- suggested slide windows
    ss <- rec$suggested_slide_summary

    if (!is.null(ss) && nrow(ss) > 0) {

      cat("\nSuggested slide windows for target rho:\n")

      for (i in seq_len(nrow(ss))) {

        cat(
          "  rho =", ss$target[i],
          ": median =", format(round(ss$median[i]), big.mark = ","), "SNPs",
          ", p90 =", format(round(ss$p90[i]), big.mark = ","), "SNPs",
          ", max =", format(round(ss$max[i]), big.mark = ","), "SNPs",
          "\n"
        )
      }

      cat(
        "  (median = typical chromosome, p90 = covers most chromosomes)\n"
      )
    }
  }

  ## ---- stored components
  cat("\nStored components:\n")

  cat(
    "  by_chr, decay_sum, decay_model, recommendation, params\n"
  )

  invisible(x)
}


#' Compute Local LD Support (ld_w)
#'
#' Computes per-SNP local LD support within a distance defined by a
#' relative LD threshold \eqn{\rho}.
#'
#' For each SNP, LD support is defined as the median \eqn{r^2} with
#' neighboring SNPs within the corresponding distance window.
#'
#' @param ld_decay Object of class \code{"ld_decay"}.
#' @param rho Relative LD threshold used to define the window.
#' @param cores Number of CPU cores.
#'
#' @return Numeric vector of LD support values (one per SNP).
#'
#' @details
#' The physical window is derived using:
#' \deqn{d = \frac{\rho}{a(1 - \rho)}}
#'
#' where \eqn{a} is the chromosome-specific decay rate.
#'
#' Requires \code{keep_el = TRUE} when calling \code{compute_LD_decay()}.
#'
#' @export
compute_ld_w <- function(
    ld_decay,
    rho = 0.95,
    cores = 1
) {

  ld_w <- unlist(parallel_apply(ld_decay$by_chr, function(chr_obj) {

    a <- ld_decay$decay_sum[Chr==chr_obj$decay_sum$Chr,a_pred]
    b <- ld_decay$decay_sum[Chr==chr_obj$decay_sum$Chr,b]

    d_window <- d_from_rho(a, rho)

    if(is.null(chr_obj$el)) stop("No edge list present")

    ld_w <- chr_obj$el[d<d_window,.(r2_median=median(r2)),by=SNP]

    ld_w[match(chr_obj$snp_ids,ld_w$SNP),r2_median]

  }, cores = cores))

  return(ld_w)
}

#' Plot LD-decay results
#'
#' S3 plotting method for objects of class \code{"ld_decay"} produced by
#' \code{compute_LD_decay()}.
#'
#' @param x Object of class \code{"ld_decay"}.
#' @param type Type of plot. One of \code{"summary"}, \code{"chr"},
#'   or \code{"recommendation"}.
#' @param chr Optional chromosome name for \code{type = "chr"}.
#'   Defaults to the first chromosome in \code{x$by_chr}.
#' @param rho Optional target rho value for highlighting in
#'   \code{type = "recommendation"}.
#' @param ask Logical; if \code{TRUE}, ask before advancing multi-panel plots.
#' @param ... Further graphical arguments passed to low-level plotting functions.
#'
#' @details
#' \describe{
#'   \item{\code{type = "summary"}}{
#'     Shows a multi-panel overview of chromosome-wise LD decay, including
#'     observed and predicted decay rates, slide-window rho coverage, and the
#'     number of informative windows retained per chromosome.
#'   }
#'   \item{\code{type = "chr"}}{
#'     Shows chromosome-specific decay estimates across windows, including
#'     window-wise parameter estimates and fitted LD-decay curves from the
#'     stratified data stored in \code{agg}.
#'   }
#'   \item{\code{type = "recommendation"}}{
#'     Visualizes suggested slide-window sizes in SNP units for the requested
#'     target rho values.freco
#'   }
#' }
#'
#' @return Invisibly returns \code{x}.
#'
#' @method plot ld_decay
#' @export
plot.ld_decay <- function(
    x,
    type = c("summary", "chr", "recommendation"),
    chr = NULL,
    rho = NULL,
    ...
){

  type <- match.arg(type)

  if (!inherits(x, "ld_decay")) {
    stop("`x` must be an object of class 'ld_decay'.")
  }

  ds <- x$decay_sum

  if (is.null(ds) || !nrow(ds)) {
    stop("No `decay_sum` found in `x`.")
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)



  if (type == "summary") {

    par(
      mfrow = c(2, 2),
      mar = c(3, 3, 2, 1),   # tighter margins
      mgp = c(2, 0.7, 0),    # axis label spacing
      tcl = -0.3             # shorter ticks
    )

    ## ------------------------------------------------------------
    ## 1. Observed vs predicted decay rate by chromosome size
    ## ------------------------------------------------------------
    ok <- is.finite(ds$chr_size) & ds$chr_size > 0 &
      is.finite(ds$a) & ds$a > 0

    plot(
      log(ds$chr_size[ok]),
      log(ds$a[ok]),
      xlab = "log chromosome size (bp)",
      ylab = "log observed decay rate (a)",
      main = "Decay rate vs chromosome size",
      pch = 19,
      ...
    )

    if ("a_pred" %in% names(ds)) {
      ord <- order(ds$chr_size)
      lines(
        log(ds$chr_size[ord]),
        log(ds$a_pred[ord]),
        lwd = 2
      )
    }

    if ("Chr" %in% names(ds)) {
      text(
        log(ds$chr_size[ok]),
        log(ds$a[ok]),
        labels = ds$Chr[ok],
        pos = 3,
        cex = 0.8
      )
    }

    ## ------------------------------------------------------------
    ## 2. Rho covered by current slide
    ## ------------------------------------------------------------
    if ("rho_slide_pred" %in% names(ds)) {
      plot(
        seq_len(nrow(ds)),
        ds$rho_slide_pred,
        xaxt = "n",
        xlab = "Chromosome",
        ylab = expression(rho),
        main = "Rho covered by current slide",
        pch = 19,
        ...
      )
      #axis(1, at = seq_len(nrow(ds)), labels = ds$Chr, las = 2)
      axis(1, at = seq_len(nrow(ds)), labels = ds$Chr, las = 2, cex.axis = 0.7)
      abline(h = median(ds$rho_slide_pred, na.rm = TRUE), lty = 2)
    } else {
      plot.new()
      title("Rho covered by current slide")
      text(0.5, 0.5, "No rho_slide_pred available")
    }

    ## ------------------------------------------------------------
    ## 3. Structured windows used per chromosome
    ## ------------------------------------------------------------
    if ("n_w_used" %in% names(ds)) {
      barplot(
        height = ds$n_w_used,
        names.arg = ds$Chr,
        las = 2,
        ylab = "Number of structured windows",
        main = "Informative windows per chromosome",
        ...
      )
    } else {
      plot.new()
      title("Informative windows per chromosome")
      text(0.5, 0.5, "No n_w_used available")
    }

    ## ------------------------------------------------------------
    ## 4. Suggested slide windows across rho targets
    ## ------------------------------------------------------------
    rec_cols <- grep("^slide_snp_rho_", names(ds), value = TRUE)

    if (length(rec_cols) > 0) {
      mat <- as.matrix(ds[, ..rec_cols])
      colnames(mat) <- sub("^slide_snp_rho_", "", rec_cols)

      matplot(
        x = seq_len(nrow(mat)),
        y = mat,
        type = "b",
        pch = 19,
        lty = 1,
        xaxt = "n",
        xlab = "Chromosome",
        ylab = "Suggested slide (SNPs)",
        main = "Recommended slide sizes"
      )
      axis(1, at = seq_len(nrow(ds)), labels = ds$Chr, las = 2)
      legend(
        "topleft",
        legend = paste0("rho=", colnames(mat)),
        lty = 1,
        pch = 19,
        bty = "n"
      )
    } else {
      plot.new()
      title("Recommended slide sizes")
      text(0.5, 0.5, "No slide recommendations available")
    }
  }

  if (type == "recommendation") {

    rec_cols <- grep("^slide_snp_rho_", names(ds), value = TRUE)

    if (!length(rec_cols)) {
      stop("No recommendation columns found in `x$decay_sum`.")
    }

    if (!is.null(rho)) {
      rho_lab <- formatC(rho, format = "f", digits = 2)
      rec_cols <- rec_cols[sub("^slide_snp_rho_", "", rec_cols) == rho_lab]

      if (!length(rec_cols)) {
        stop("Requested `rho` not found among recommendation columns.")
      }
    }

    mat <- as.matrix(ds[, ..rec_cols])
    colnames(mat) <- sub("^slide_snp_rho_", "", rec_cols)

    if (ncol(mat) == 1) {
      plot(
        seq_len(nrow(ds)),
        mat[, 1],
        type = "b",
        pch = 19,
        xaxt = "n",
        xlab = "Chromosome",
        ylab = "Suggested slide (SNPs)",
        main = paste0("Recommended slide size (rho=", colnames(mat)[1], ")"),
        ...
      )
      axis(1, at = seq_len(nrow(ds)), labels = ds$Chr, las = 2)
      abline(h = median(mat[, 1], na.rm = TRUE), lty = 2)
    } else {

      cols <- c("salmon","steelblue","black","firebrick","darkorange")[seq_len(ncol(mat))]

      matplot(
        x = seq_len(nrow(ds)),
        y = mat,
        type = "b",
        pch = 19,
        col = cols,
        lty = 1,
        xaxt = "n",
        xlab = "Chromosome",
        ylab = "Suggested slide (SNPs)",
        main = "Recommended slide sizes",
        ...
      )
      axis(1, at = seq_len(nrow(ds)), labels = ds$Chr, las = 2)
      legend(
        "topleft",
        legend = paste0("rho=", colnames(mat)),
        lty = 1,
        pch = 19,
        col = cols,
        bty = "n"
      )
    }
  }

  if (type == "chr") {

    if (is.null(chr)) {
      chr <- names(x$by_chr)[1]
    }

    if (!chr %in% names(x$by_chr)) {
      stop("`chr` not found in `x$by_chr`.")
    }

    chr_obj <- x$by_chr[[chr]]
    decay   <- chr_obj$decay

    if (is.null(decay) || !nrow(decay)) {
      stop("No chromosome-specific decay data found for `chr`.")
    }


    par(
      mfrow = c(2, 2),
      mar = c(3, 3, 2, 1),   # tighter margins
      mgp = c(2, 0.7, 0),    # axis label spacing
      tcl = -0.3             # shorter ticks
    )

    ## ------------------------------------------------------------
    ## 1. Window-wise a
    ## ------------------------------------------------------------
    if ("a" %in% names(decay)) {
      mid <- rowMeans(decay[, .(start, end)], na.rm = TRUE)

      plot(
        mid,
        decay$a,
        xlab = "Genomic position (bp)",
        ylab = "Decay rate (a)",
        main = paste0(chr, ": window-wise decay rate"),
        pch = 19,
        ...
      )

      chr_sum <- chr_obj$decay_sum
      if (!is.null(chr_sum) && "a" %in% names(chr_sum)) {
        abline(h = chr_sum$a[1], lty = 2, lwd = 2)
      }
    } else {
      plot.new()
      title(paste0(chr, ": window-wise decay rate"))
      text(0.5, 0.5, "No window-wise a estimates")
    }

    ## ------------------------------------------------------------
    ## 2. Window-wise c
    ## ------------------------------------------------------------
    if ("c" %in% names(decay)) {
      mid <- rowMeans(decay[, .(start, end)], na.rm = TRUE)

      plot(
        mid,
        decay$c,
        xlab = "Genomic position (bp)",
        ylab = "Short-range LD (c)",
        main = paste0(chr, ": window-wise short-range LD"),
        pch = 19,
        ...
      )

      chr_sum <- chr_obj$decay_sum
      if (!is.null(chr_sum) && "c" %in% names(chr_sum)) {
        abline(h = chr_sum$c[1], lty = 2, lwd = 2)
      }
    } else {
      plot.new()
      title(paste0(chr, ": window-wise short-range LD"))
      text(0.5, 0.5, "No window-wise c estimates")
    }

    ## ------------------------------------------------------------
    ## 3. Contrast by window
    ## ------------------------------------------------------------
    if (all(c("c", "b") %in% names(decay))) {
      mid <- rowMeans(decay[, .(start, end)], na.rm = TRUE)
      contrast <- decay$c - decay$b

      plot(
        mid,
        contrast,
        xlab = "Genomic position (bp)",
        ylab = expression(c - b),
        main = paste0(chr, ": LD contrast across windows"),
        pch = 19,
        ...
      )
      abline(h = 0.05, lty = 2)
    } else {
      plot.new()
      title(paste0(chr, ": LD contrast across windows"))
      text(0.5, 0.5, "No c/b estimates available")
    }

    ## ------------------------------------------------------------
    ## 4. Stratified decay curves from stored agg data
    ## ------------------------------------------------------------
    agg_list <- decay$agg[!vapply(decay$agg, is.null, logical(1))]
    agg_list <- agg_list[lengths(agg_list) > 0]

    if (length(agg_list) > 0) {

      xlim <- range(unlist(lapply(agg_list, function(z) z$d_mid)), na.rm = TRUE)
      ylim <- range(unlist(lapply(agg_list, function(z) z$r2_q)), na.rm = TRUE)

      plot(
        NA,
        xlim = xlim,
        ylim = ylim,
        xlab = "Distance (bp)",
        ylab = expression(r^2),
        main = paste0(chr, ": fitted decay curves"),
        ...
      )

      for (i in seq_along(agg_list)) {
        agg <- agg_list[[i]]
        points(agg$d_mid, agg$r2_q, pch = 16, cex = 0.6)

        if (all(c("a", "c", "b") %in% names(decay))) {
          ai <- decay$a[match(i, which(!vapply(decay$agg, is.null, logical(1))))]
          ci <- decay$c[match(i, which(!vapply(decay$agg, is.null, logical(1))))]
          bi <- decay$b[match(i, which(!vapply(decay$agg, is.null, logical(1))))]

          if (all(is.finite(c(ai, ci, bi)))) {
            dseq <- seq(min(agg$d_mid), max(agg$d_mid), length.out = 200)
            yhat <- bi + (ci - bi) / (1 + ai * dseq)
            lines(dseq, yhat)
          }
        }
      }

      chr_sum <- chr_obj$decay_sum
      if (!is.null(chr_sum) && all(c("a", "c", "b") %in% names(chr_sum))) {
        dseq <- seq(xlim[1], xlim[2], length.out = 300)
        yhat <- chr_sum$b[1] + (chr_sum$c[1] - chr_sum$b[1]) / (1 + chr_sum$a[1] * dseq)
        lines(dseq, yhat, lwd = 3,col="salmon")
      }

    } else {
      plot.new()
      title(paste0(chr, ": fitted decay curves"))
      text(0.5, 0.5, "No stored stratified fits available")
    }
  }

  invisible(x)
}

