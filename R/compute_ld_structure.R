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
compute_ld_structure <- function(gds,
                                 q = 0.95,
                                 n_sub = 5000,
                                 window_size = 1e6,
                                 step_size = 5e5,
                                 cores = 1,
                                 dist_unit = 5000,
                                 K_target = 1000,
                                 n_win = 1000,
                                 prob_robust = 0.5,
                                 use="robust",
                                 compress=TRUE,
                                 slide_win = 1000,
                                 max_rho=0.99) {

  ids  <- .read_gds_ids(gds)
  chrs <- unique(ids$snp_chr)

  b <- estimate_background_ld(gds, n_sub = n_sub, q = q)

  out <- list(by_chr = list())
  #ch <- "Chr1"
  for (ch in chrs) {

    message("Processing ", ch)

    chr_idx <- which(ids$snp_chr == ch)
    pos_chr <- ids$snp_pos[chr_idx]

    # ---- 1. Estimate decay
    decay <- estimate_decay_chr(
      gds, chr_idx, pos_chr, b,
      window_size = window_size,
      step_size = step_size,
      q = q,
      dist_unit = dist_unit,
      cores = cores,
      n_win = n_win
    )

    # ---- 2. Summarize decay
    decay_sum <- summarize_decay(decay, prob_robust)
    # add chromosome
    decay_sum <- lapply(decay_sum,function(x) cbind(Chr=ch,x))

    if (is.null(decay_sum))
      next

    # ---- 3. Derive LD scale
    decay_sum$robust[,W_max:=d_from_rho(a=decay_sum$robust$a,  d0 = decay_sum$robust$d0, rho = max_rho)]
    decay_sum$median[,W_max:=d_from_rho(a=decay_sum$median$a, d0 = decay_sum$median$d0, rho = max_rho)]

    edges <- get_el(gds = gds,
                    idx = chr_idx,
                    slide_win_ld = slide_win,
                    cores = cores)
    edges <- edges[d < decay_sum[[use]]$W_max]

    edges[,r2:=r2-min(r2)]

    if(compress){
      message("Compressing edge list \n")
      #focal <- ids$snp_id[chr_idx][1]
      cum_histograms2d <- parallel::mclapply(ids$snp_id[chr_idx],function(focal) {

        build_2d_ld_hist(

            el = edges,
            focal = focal,
            dist_unit = 1e6,
            r2_unit = 0.001

          )},mc.cores=cores)
    }

    out$by_chr[[ch]] <- list(
      snp_ids = ids$snp_id[chr_idx],
      edges   = if(!compress) edges else NULL,
      cum_histograms2d   = if(compress) cum_histograms2d else NULL,
      decay   = decay,
      summary = decay_sum
    )

  }
  out$summary <- list()
  out$summary$median <- rbindlist(lapply(out$by_chr,function(x)x$summary$median))
  out$summary$robust <- rbindlist(lapply(out$by_chr,function(x)x$summary$robust))

  out$params <- list(q = q,
                     n_sub = n_sub,
                     window_size = window_size,
                     step_size = step_size,
                     dist_unit = dist_unit,
                     K_target = K_target,
                     n_win = n_win,
                     prob_robust = prob_robust,
                     compress = compress,
                     use=use)

  class(out) <- "ld_structure"
  gc()
  out
}
#' Build LD Edge List for a Chromosome
#'
#' Constructs a filtered LD edge list for LD-weighted statistics
#' within a chromosome.
#'
#' Edges are retained if:
#' \itemize{
#'   \item \eqn{r^2 > b}
#'   \item \eqn{|d| \le W_{\max}}
#' }
#'
#' Adaptive thinning ensures approximately \code{K_target}
#' pairwise comparisons per SNP.
#'
#' @param gds GDS file handle.
#' @param chr_idx Indices of SNPs on the chromosome.
#' @param snp_pos SNP physical positions (bp).
#' @param b Background LD.
#' @param W_max Maximum LD distance derived from decay.
#' @param K_target Target number of pairwise LD comparisons per SNP.
#' @param cores Number of CPU cores.
#'
#' @return A \code{data.table} of LD edges.
#'
#' @keywords internal
build_ld_edges_chr <- function(gds,
                               chr_idx,
                               snp_pos,
                               b,
                               W_max,
                               K_target = 1000,
                               cores = 1) {

  thin_res <- adaptive_thin_pairs(
    snp_pos,
    W_max = W_max,
    K_target = K_target
  )

  edge_list <- vector("list", thin_res$thin_factor)

  for (o in seq_len(thin_res$thin_factor)) {

    keep_local <- thin_res$keep[[o]]

    edge_list[[o]] <- get_el(
      gds,
      idx = chr_idx[keep_local],
      slide_win_ld = -1,
      cores = cores
    )[r2 > b & abs(pos2 - pos1) <= W_max]
  }

  edge_list <- rbindlist(edge_list)
  edge_list[,r2:=r2-b]
  return(edge_list)
}
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
summarize_decay <- function(decay_dt, prob_robust = 0.5) {

  decay_valid <- decay_dt[!is.na(a) & a > 0 & !is.na(c)]

  if (nrow(decay_valid) < 5)
    return(NULL)

  # median summary
  median_summary <- decay_valid[, .(
    c  = median(c,  na.rm = TRUE),
    a  = median(a,  na.rm = TRUE),
    d0 = median(d0, na.rm = TRUE),
    b  = b[1]
  )]

  # robust summary (exclude slow decay windows)
  a_cut <- quantile(decay_valid$a,
                    probs = prob_robust,
                    na.rm = TRUE)

  robust <- decay_valid[a >= a_cut]

  robust_summary <- robust[, .(
    c  = median(c,  na.rm = TRUE),
    a  = median(a,  na.rm = TRUE),
    d0 = median(d0, na.rm = TRUE),
    b  = b[1]
  )]

  list(
    median = median_summary,
    robust = robust_summary
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
                               window_size = 1e6,
                               step_size = 5e5,
                               q = 0.95,
                               dist_unit = 5000,
                               cores = 1,
                               n_win = 1000) {

  thin <- adaptive_thinning(snp_pos, W_LD = window_size, n_win = n_win)

  el <- get_el(
    gds,
    idx = chr_idx[thin$keep],
    slide_win_ld = thin$slide_win_ld,
    cores = cores
  )

  min_pos <- min(c(el$pos1, el$pos2))
  max_pos <- max(c(el$pos1, el$pos2))

  starts <- seq(min_pos, max_pos - window_size, by = step_size)
  ends   <- starts + window_size

  decay <- rbindlist(lapply(seq_along(starts), function(i) {

    sub <- el[
      pos1 >= starts[i] & pos1 < ends[i] &
        pos2 >= starts[i] & pos2 < ends[i]
    ]

    coefs <- tryCatch(
      coef_ld_dec(sub, q = q, dist_unit = dist_unit, b = b),
      error = function(e) NULL
    )

    if (is.null(coefs)) {
      data.table(start = starts[i], end = ends[i],
                 c = NA, a = NA, d0 = NA)
    } else {
      data.table(start = starts[i], end = ends[i],
                 t(coefs))
    }
  }))

  decay[, b := b]

  decay
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
#' @param n_win Target SNP count per window.
#'
#' @return List with retained SNP indices and slide window size.
#'
#' @keywords internal
adaptive_thinning <- function(pos,
                              W_LD = 1e7,
                              n_win = 1000) {

  target_spacing <- W_LD / n_win
  mean_spacing <- mean(diff(pos))

  if (mean_spacing >= target_spacing) {
    keep <- seq_along(pos)
  } else {
    thin_factor <- ceiling(target_spacing / mean_spacing)
    keep <- seq(1, length(pos), by = thin_factor)
  }

  list(
    keep = keep,
    slide_win_ld = floor(W_LD / mean_spacing)
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
coef_ld_dec <- function(el_ld, q = 0.95, dist_unit = 5000, b = 0.05) {

  el_ld[, dist := abs(pos2 - pos1)]
  el_ld[, dist_bin := floor(dist / dist_unit)]

  bin_dt <- el_ld[, .(
    r2_q = quantile(r2, na.rm = TRUE, prob = q),
    w = .N,
    dist_mean = mean(dist)
  ), by = dist_bin][order(dist_bin)]

  if (nrow(bin_dt) < 5)
    return(c(c = NA, a = NA, d0 = NA))

  # ensure first bin above background
  c_start <- max(b + 1e-6, min(1, bin_dt$r2_q[1]))

  starts <- list(
    c  = c_start,
    a  = 1 / (median(bin_dt$dist_mean[bin_dt$dist_mean > 0], na.rm = TRUE) + 1e-6),
    d0 = 0
  )

  fit <- tryCatch({

    nls(
      r2_q ~ b + (c - b) / (1 + a * pmax(dist_mean - d0, 0)),
      data = bin_dt,
      start = starts,
      weights = w,
      algorithm = "port",
      lower = c(b + 1e-6, 0, 0),
      upper = c(1, Inf, max(bin_dt$dist_mean)),
      control = nls.control(warnOnly = TRUE)
    )

  }, error = function(e) NULL)

  if (is.null(fit))
    return(c(c = NA, a = NA, d0 = NA))

  coef(fit)
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
d_from_rho <- function(a, rho, d0 = 0) {
  d0 + (1 / a) * (1 / (1 - rho) - 1)
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

  if (!inherits(x, "ld_structure"))
    stop("Object must be of class 'ld_structure'.")

  chrs <- names(x$by_chr)

  n_chr <- length(chrs)

  n_snps <- sum(sapply(x$by_chr, function(z)
    length(z$snp_ids)))

  n_edges <- sum(sapply(x$by_chr, function(z)
    if (!is.null(z$edges)) nrow(z$edges) else 0))

  cat("LD Structure Object\n")
  cat("-------------------\n")
  cat("Chromosomes:", n_chr, "\n")
  cat("Total SNPs:", n_snps, "\n")
  cat("Total LD edges:", n_edges, "\n\n")

  # robust decay summary
  robust <- rbindlist(lapply(x$by_chr, function(z)
    z$summary$robust), fill = TRUE)

  if (!is.null(robust) && nrow(robust) > 0) {

    cat("Robust LD-decay parameters (median across chromosomes):\n")

    cat("  a  =", signif(median(robust$a,  na.rm = TRUE), 3), "\n")
    cat("  c  =", signif(median(robust$c,  na.rm = TRUE), 3), "\n")
    cat("  d0 =", signif(median(robust$d0, na.rm = TRUE), 3), "\n")
    cat("  b  =", signif(median(robust$b,  na.rm = TRUE), 3), "\n")
  }
  cat("\n")
  cat("Parameters used :\n")

  cat("  q           =", x$params$q, "\n")
  cat("  n_sub.      =", x$params$n_sub, "\n")
  cat("  window_size =", x$params$window_size, "\n")
  cat("  step_size   =", x$params$step_size, "\n")
  cat("  dist_unit   =", x$params$dist_unit, "\n")
  cat("  K_target    =", x$params$K_target, "\n")
  cat("  n_win       =", x$params$n_win, "\n")
  cat("  prob_robust =", x$params$prob_robust, "\n")
  cat("  Compressed  =", x$params$compress, "\n")
  cat("  use         =", x$params$use, "\n")


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



#el = edges
build_2d_ld_hist <- function(el, focal,
                             dist_unit = 1e6,
                             r2_unit   = 0.001) {


  # 1️⃣ Symmetric neighbor extraction
  dt <- rbind(
    el[SNP1 == focal,.(pos = pos1, pos_other = pos2, r2)],
    el[SNP2 == focal,.(pos = pos2, pos_other = pos1, r2)]
  )

  if (nrow(dt) == 0)
    return(NULL)

  # 2️⃣ Distance + side
  dt[, d := abs(pos_other - pos)]
  dt[, side := ifelse(pos_other > pos, "R", "L")]

  # 3️⃣ Compute maxL / maxR
  maxL <- if (any(dt$side == "L")) max(dt$d[dt$side == "L"]) else 0
  maxR <- if (any(dt$side == "R")) max(dt$d[dt$side == "R"]) else 0

  S <- min(maxL, maxR)

  # 4️⃣ Identify long side
  if (maxL > maxR) {
    long_side <- "L"
  } else if (maxR > maxL) {
    long_side <- "R"
  } else {
    long_side <- NA_character_
  }

  # 5️⃣ Reweight: duplicate tail beyond S on long side
  if (!is.na(long_side) && S > 0) {

    tail_dt <- dt[side == long_side & d > S]

    if (nrow(tail_dt) > 0)
      dt <- rbind(dt, tail_dt)
  }
  #dist_unit <- 1000
  # 6️⃣ Bin distance + r2
  dt[, dist_bin := floor(d / dist_unit) * dist_unit]
  dt[, r2_bin   := round(r2 / r2_unit) * r2_unit]


  # 7️⃣ Build 2D histogram
  hist2d <- dt[, .N, by = .(dist_bin, r2_bin)]

  hist_mat <- dcast(hist2d,
                    dist_bin ~ r2_bin,
                    value.var = "N",
                    fill = 0)

  mat <- as.matrix(hist_mat[, -1])
  rownames(mat) <- hist_mat$dist_bin

  # Ensure sorted rows
  ord <- order(as.numeric(rownames(mat)))
  mat <- mat[ord, , drop = FALSE]

  # 8️⃣ Convert to cumulative along distance


  cum_mat <- apply(mat, 2, cumsum)

  rownames(cum_mat) <- rownames(mat)
  colnames(cum_mat) <- colnames(mat)

  return(cum_mat)

}

#b

el <- edges
compute_ld_auc_norm <- function(el) {
  print(focal)
  # 1️⃣ Symmetric neighbor extraction
  dt <- rbind(
    el[,.(SNP=SNP1, pos = pos1, pos_other = pos2, r2,d)],
    el[,.(SNP=SNP2, pos = pos2, pos_other = pos1, r2,d)]
  )

  if (nrow(dt) == 0)
    return(NULL)

  # 2️⃣ Distance + side
  dt[, side := ifelse(pos_other > pos, "R", "L")]

  # 3️⃣ Compute maxL / maxR
  maxL <- if (any(dt$side == "L")) max(dt$d[dt$side == "L"]) else 0
  maxR <- if (any(dt$side == "R")) max(dt$d[dt$side == "R"]) else 0

  S <- min(maxL, maxR)

  # 4️⃣ Identify long side
  if (maxL > maxR) {
    long_side <- "L"
  } else if (maxR > maxL) {
    long_side <- "R"
  } else {
    long_side <- NA_character_
  }

  # 5️⃣ Reweight: duplicate tail beyond S on long side
  if (!is.na(long_side) && S > 0) {

    tail_dt <- dt[side == long_side & d > S]

    if (nrow(tail_dt) > 0)
      dt <- rbind(dt, tail_dt)
  }

  AUC_LD <- dt[, .(
    LD_AUC = mean(pmax(r2 - b, 0))
  ), by = SNP]


  res <- data.table::data.table(
    SNP = map[Chr==ch,marker]
  )

  AUC_LD <- AUC_LD[res, on="SNP"]


  # rho_w <- 0.95
  # rho_Ws <- c()
  # a_chr <-decay_sum$robust[Chr == ch, a]
  # d0_chr <- decay_sum$robust[Chr == ch, d0]
  # d_th  <- d_from_rho(a_chr, rho_w,d0 = d0_chr)
  #
  # ldw <- dt[d<d_th, .(
  #   ld_w = median(r2,na.rm = TRUE)
  # ), by = SNP]
  #
  #
  # res <- data.table::data.table(
  #   SNP = map[Chr==ch,marker]
  # )
  #
  # ldw <- ldw[res, on="SNP"]
  #
  # map[Chr==ch,LD_AUC := AUC_LD[,LD_AUC]]
  # map[Chr==ch,ld_w := ldw[,ld_w]]
  #
  # map[Chr==ch,plot(ld_w,max_LD_with_QTN)]
  #
  #
  # map[Chr==ch,summary(lm(max_LD_with_QTN~LD_AUC))]
  # map[Chr==ch,summary(lm(max_LD_with_QTN~LD_AUC))]
  #
  # plot(AUC_LD[,LD_AUC],map[Chr==ch,max_LD_with_QTN])
  # plot(AUC_LD[,LD_AUC],map[Chr==ch,max_LD_with_QTN])
  # abline(v=0.025)
  #AUC_LD[,table(LD_AUC>0.025)]
  #plot(AUC_LD[,ld_w],map[Chr==ch,max_LD_with_QTN])

  return(AUC_LD)
}

sapply(SNP)
AUC_LD <-

AUC_LD <-  rbindlist(parallel::mclapply(ids$snp_id[chr_idx],compute_ld_auc_norm,el,mc.cores=8))





lapply(compute_ld_auc_norm(dt))
compute_ld_auc <- function(hist_mat, b = 0) {

  r2_vals <- as.numeric(colnames(hist_mat))

  # subtract background
  r2_excess <- pmax(r2_vals - b, 0)

  # weighted sum
  ld_auc <- sum(hist_mat %*% r2_excess)



  return(as.numeric(ld_auc))
}

focal <- ids$snp_id[chr_idx][10000]
r2_counts <- cum_mat[1,]
ids$snp_id[c]

compress_cum_mat <- function(cum_mat) {
  # drop zero columns
  keep_cols <- colSums(cum_mat) > 0
  cum_mat <- cum_mat[, keep_cols, drop = FALSE]

  # drop redundant rows
  keep_rows <- c(TRUE, rowSums(abs(diff(cum_mat))) > 0)
  cum_mat <- cum_mat[keep_rows, , drop = FALSE]
}
cum_mat <- compress_cum_mat(cum_mat)
dim(cum_mat)
plot(apply(cum_mat,1,function(r2_counts){
  median_pos <- ceiling(sum(r2_counts) / 2)
  cum_r2 <- cumsum(r2_counts)
  j <- which(cum_r2 >= median_pos)[1]
  as.numeric(colnames(cum_mat))[j]
}),type="l")

medians <- apply(cum_mat,1,function(r2_counts){
  median_pos <- ceiling(sum(r2_counts) / 2)
  cum_r2 <- cumsum(r2_counts)
  j <- which(cum_r2 >= median_pos)[1]
  as.numeric(colnames(cum_mat))[j]
})
keep_rows <- c(TRUE, abs(diff(medians)) > 0)

plot(names(medians[keep_rows]),medians[keep_rows],type="l")
