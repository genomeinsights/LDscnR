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
compute_LD_decay <- function(
    gds,
    el_data_folder = NULL,
    ## for LD-decay and bg
    q = 0.95,
    ## for bg
    n_sub_bg = 5000,
    ## for decay
    n_win_decay = 20,
    overlap = 0.5,
    max_SNPs_decay = Inf,
    prob_robust = 0.95,
    max_pairs = 5000,
    n_strata = 20,
    keep_el = FALSE,
    ## for ld_int / LD support
    slide = 1000,                  # in SNPs
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
  for (ch in chrs) {

    chr_idx  <- which(ids$snp_chr == ch)
    pos_chr  <- ids$snp_pos[chr_idx]
    snps_chr <- ids$snp_id[chr_idx]

    cat(ch, " - ")

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

  message("Predicting a from chromosome size")
  decay_sum <- data.table::rbindlist(
    lapply(out_by_chr, function(x) x$decay_sum),
    fill = TRUE
  )

  ## robust background decay model
  d_mod <- MASS::rlm(log(a) ~ log(chr_size), data = decay_sum)

  decay_sum[, a_pred := exp(predict(d_mod, decay_sum))]
  decay_sum[, c_pred := median(c, na.rm = TRUE)]

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



#' Print Method for ld_structure Objects
#'
#' Displays a concise summary of the LD structure object.
#'
#' @param x An object of class `"ld_structure"`.
#' @param ... Additional arguments (unused).
#'
#' @export
print.ld_decay <- function(x, digits = 3, ...) {

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

summary.ld_decay <- function(object, ...) {
  object$recommendation$suggested_slide_summary
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

compute_ld_int_fill <- function(
    el,
    snp_ids,
    a,
    n_bins = 20,
    rho_max = 0.90,
    edge_mode = c("none", "trim", "fill")
) {

  edge_mode <- match.arg(edge_mode)

  if (nrow(el) == 0L) {
    return(setNames(rep(NA_real_, length(snp_ids)), snp_ids))
  }

  if (!is.finite(a) || a <= 0) {
    stop("`a` must be a positive finite number.")
  }

  if (!is.numeric(rho_max) || length(rho_max) != 1L || rho_max <= 0 || rho_max >= 1) {
    stop("`rho_max` must be a single number in (0, 1).")
  }

  n_bins <- as.integer(n_bins)
  if (!is.finite(n_bins) || n_bins < 1L) {
    stop("`n_bins` must be a positive integer.")
  }

  if (!all(c("SNP", "d", "r2", "pos1", "pos2") %in% names(el))) {
    stop("`el` must contain columns: SNP, d, r2, pos1, pos2")
  }

  # truncate to rho support
  d_cap <- rho_max / (a * (1 - rho_max))
  el_use <- el[d <= d_cap, .(SNP, d, r2, pos1, pos2)]

  if (nrow(el_use) == 0L) {
    return(setNames(rep(NA_real_, length(snp_ids)), snp_ids))
  }

  # left/right leg
  el_use[, side := data.table::fifelse(pos2 > pos1, "R", "L")]

  # optional symmetric trimming to common support
  if (edge_mode == "trim") {
    support_dt <- el_use[, .(
      maxL = if (any(side == "L")) max(d[side == "L"]) else 0,
      maxR = if (any(side == "R")) max(d[side == "R"]) else 0
    ), by = SNP]

    support_dt[, S := pmin(maxL, maxR)]

    el_use <- support_dt[el_use, on = "SNP"]
    el_use <- el_use[d <= S, .(SNP, d, r2, pos1, pos2, side)]
  }

  if (nrow(el_use) == 0L) {
    return(setNames(rep(NA_real_, length(snp_ids)), snp_ids))
  }

  # equal-rho bins via direct mapping
  rho <- (a * el_use$d) / (1 + a * el_use$d)
  dist_idx <- as.integer(floor(n_bins * rho / rho_max) + 1L)
  dist_idx[dist_idx < 1L] <- 1L
  dist_idx[dist_idx > n_bins] <- n_bins
  el_use[, dist_idx := dist_idx]

  # side-specific bin medians
  shell_side <- el_use[, .(
    n_pairs = .N,
    med_r2 = median(r2,na.rm = TRUE)
  ), by = .(SNP, side, dist_idx)]

  shell_side <- shell_side[is.finite(med_r2)]

  if (nrow(shell_side) == 0L) {
    return(setNames(rep(NA_real_, length(snp_ids)), snp_ids))
  }

  # wide table: one row per SNP x bin, columns for left and right medians
  shell_wide <- data.table::dcast(
    shell_side,
    SNP + dist_idx ~ side,
    value.var = "med_r2"
  )

  # bin-level summary depending on edge handling
  if (edge_mode == "none" || edge_mode == "trim") {
    # use whatever is observed in that bin
    shell_wide[, bin_val := data.table::fifelse(
      !is.na(L) & !is.na(R), (L + R) / 2,
      data.table::fifelse(!is.na(L), L, R)
    )]
  }

  if (edge_mode == "fill") {
    # if one side is missing, fill it from the other side
    shell_wide[, bin_val := data.table::fifelse(
      !is.na(L) & !is.na(R), (L + R) / 2,
      data.table::fifelse(!is.na(L), L, R)
    )]
    # algebraically same bin value, but conceptually this is a symmetric fill
    # because missing leg bins are completed by the observed leg before combining
  }

  shell_wide <- shell_wide[is.finite(bin_val)]

  if (nrow(shell_wide) == 0L) {
    return(setNames(rep(NA_real_, length(snp_ids)), snp_ids))
  }

  # final DPI: median across occupied rho-bins
  LD_dt <- shell_wide[, .(
    LD_int = median(bin_val)
  ), by = SNP]

  result <- setNames(rep(NA_real_, length(snp_ids)), snp_ids)
  result[LD_dt$SNP] <- LD_dt$LD_int

  result
}

compute_ld_int <- function(
    el,
    snp_ids,
    a,
    b = NULL,
    c = NULL,
    n_bins = 20,
    rho_max = 0.90,
    min_pairs = 3,
    edge_mode = c("none", "trim", "fill"),
    scale = c("raw", "excess", "contrast"),
    clip_contrast = TRUE,
    return_support = FALSE
) {

  edge_mode <- match.arg(edge_mode)
  scale <- match.arg(scale)

  if (nrow(el) == 0L) {
    result <- setNames(rep(NA_real_, length(snp_ids)), snp_ids)
    if (return_support) {
      return(list(
        LD_int = result,
        support = data.table::data.table()
      ))
    }
    return(result)
  }

  if (!is.finite(a) || a <= 0) {
    stop("`a` must be a positive finite number.")
  }

  if (scale %in% c("excess", "contrast") && !is.finite(b)) {
    stop("For `scale = 'excess'` or `scale = 'contrast'`, `b` must be provided.")
  }

  if (scale == "contrast") {
    if (!is.finite(c)) stop("For `scale = 'contrast'`, `c` must be provided.")
    if ((c - b) <= 0) stop("For `scale = 'contrast'`, `c - b` must be > 0.")
  }

  if (!is.numeric(rho_max) || length(rho_max) != 1L || rho_max <= 0 || rho_max >= 1) {
    stop("`rho_max` must be a single number in (0, 1).")
  }

  n_bins <- as.integer(n_bins)
  if (!is.finite(n_bins) || n_bins < 1L) {
    stop("`n_bins` must be a positive integer.")
  }

  req_cols <- c("SNP", "d", "r2", "pos1", "pos2")
  if (!all(req_cols %in% names(el))) {
    stop("`el` must contain columns: ", paste(req_cols, collapse = ", "))
  }

  ## ---- truncate to rho support
  d_cap <- rho_max / (a * (1 - rho_max))
  el_use <- el[d <= d_cap, .(SNP, d, r2, pos1, pos2)]

  if (nrow(el_use) == 0L) {
    result <- setNames(rep(NA_real_, length(snp_ids)), snp_ids)
    if (return_support) {
      return(list(
        LD_int = result,
        support = data.table::data.table()
      ))
    }
    return(result)
  }

  ## ---- left/right leg
  el_use[, side := data.table::fifelse(pos2 > pos1, "R", "L")]

  ## ---- optional symmetric trimming
  if (edge_mode == "trim") {
    support_trim <- el_use[, .(
      maxL = if (any(side == "L")) max(d[side == "L"]) else 0,
      maxR = if (any(side == "R")) max(d[side == "R"]) else 0
    ), by = SNP]

    support_trim[, S := pmin(maxL, maxR)]

    el_use <- support_trim[el_use, on = "SNP"]
    el_use <- el_use[d <= S, .(SNP, d, r2, pos1, pos2, side)]
  }

  if (nrow(el_use) == 0L) {
    result <- setNames(rep(NA_real_, length(snp_ids)), snp_ids)
    if (return_support) {
      return(list(
        LD_int = result,
        support = data.table::data.table()
      ))
    }
    return(result)
  }

  ## ---- equal-rho bins via direct mapping
  rho <- (a * el_use$d) / (1 + a * el_use$d)
  dist_idx <- as.integer(floor(n_bins * rho / rho_max) + 1L)
  dist_idx[dist_idx < 1L] <- 1L
  dist_idx[dist_idx > n_bins] <- n_bins
  el_use[, dist_idx := dist_idx]

  ## ---- side-specific bin medians
  shell_side <- el_use[, .(
    n_pairs = .N,
    med_r2 = if (.N >= min_pairs) median(r2) else NA_real_
  ), by = .(SNP, side, dist_idx)]

  shell_side <- shell_side[is.finite(med_r2)]

  if (nrow(shell_side) == 0L) {
    result <- setNames(rep(NA_real_, length(snp_ids)), snp_ids)
    if (return_support) {
      return(list(
        LD_int = result,
        support = data.table::data.table()
      ))
    }
    return(result)
  }

  ## ---- wide table: one row per SNP x bin, columns L and R
  shell_wide <- data.table::dcast(
    shell_side,
    SNP + dist_idx ~ side,
    value.var = "med_r2"
  )

  ## make sure missing side columns exist
  if (!"L" %in% names(shell_wide)) shell_wide[, L := NA_real_]
  if (!"R" %in% names(shell_wide)) shell_wide[, R := NA_real_]

  ## ---- combine sides within bins
  shell_wide[, bin_val := data.table::fifelse(
    !is.na(L) & !is.na(R), (L + R) / 2,
    data.table::fifelse(!is.na(L), L, R)
  )]

  shell_wide <- shell_wide[is.finite(bin_val)]

  if (nrow(shell_wide) == 0L) {
    result <- setNames(rep(NA_real_, length(snp_ids)), snp_ids)
    if (return_support) {
      return(list(
        LD_int = result,
        support = data.table::data.table()
      ))
    }
    return(result)
  }

  ## ---- transform bin values if requested
  if (scale == "raw") {
    shell_wide[, bin_val_scaled := bin_val]
  }

  if (scale == "excess") {
    shell_wide[, bin_val_scaled := pmax(0, bin_val - b)]
  }

  if (scale == "contrast") {
    shell_wide[, bin_val_scaled := (bin_val - b) / (c - b)]
    if (clip_contrast) {
      shell_wide[, bin_val_scaled := pmin(1, pmax(0, bin_val_scaled))]
    }
  }

  ## ---- fill mode: duplicate missing leg at the bin-summary level
  if (edge_mode == "fill") {
    shell_long <- shell_wide[, .(
      SNP,
      dist_idx,
      val1 = data.table::fifelse(!is.na(L), bin_val_scaled, bin_val_scaled),
      val2 = data.table::fifelse(!is.na(R), bin_val_scaled, bin_val_scaled),
      L_present = !is.na(L),
      R_present = !is.na(R)
    )]

    ## if both sides present, keep two equal representatives
    ## if only one side present, fill the missing side with the observed bin value
    shell_long <- data.table::melt(
      shell_long,
      id.vars = c("SNP", "dist_idx", "L_present", "R_present"),
      measure.vars = c("val1", "val2"),
      value.name = "bin_val_scaled"
    )

    LD_dt <- shell_long[, .(
      LD_int = median(bin_val_scaled, na.rm = TRUE)
    ), by = SNP]

  } else {
    LD_dt <- shell_wide[, .(
      LD_int = median(bin_val_scaled, na.rm = TRUE)
    ), by = SNP]
  }

  result <- setNames(rep(NA_real_, length(snp_ids)), snp_ids)
  result[LD_dt$SNP] <- LD_dt$LD_int

  if (!return_support) {
    return(result)
  }

  support_dt <- shell_wide[, .(
    n_bins_obs   = .N,
    n_bins_both  = sum(!is.na(L) & !is.na(R)),
    n_bins_Lonly = sum(!is.na(L) &  is.na(R)),
    n_bins_Ronly = sum( is.na(L) & !is.na(R))
  ), by = SNP]

  list(
    LD_int = result,
    support = support_dt
  )
}


#' @export
compute_DPI <- function(ld_decay,bins=1000,rho_max=0.99,edge_mode="fill",scale=c("raw", "excess", "contrast"),cores=1){

  dpi <- unlist(parallel_apply(ld_decay$by_chr, function(chr_obj) {

    a <- ld_decay$decay_sum[Chr==chr_obj$decay_sum$Chr,a]

    compute_ld_int(chr_obj$el, snp_ids=chr_obj$snp_ids, a,n_bins = bins,rho_max = rho_max,edge_mode = edge_mode)

  }, cores = cores))

  return(dpi)
}
