#' LD-decay model with short-range plateau
#'
#' To accommodate situations where LD remains approximately constant
#' over short physical distances (e.g., due to limited recombination or
#' dense marker spacing), a plateau parameter \eqn{d_0} is
#' included in the decay model. The LD-decay curve is then defined as
#'
#' \deqn{
#' r^2(d) = b + \frac{c - b}{1 + a \max(d - d_0, 0)},
#' }
#'
#' where \eqn{b} is the background LD level, \eqn{c} is the short-range
#' LD plateau (bounded by 1), \eqn{a} controls the decay rate, and
#' \eqn{d_0} represents the distance up to which LD remains approximately
#' constant.
#'
#' The use of \code{pmax(d - d_0, 0)} ensures a continuous transition from
#' the plateau region (\eqn{d \le d_0}) to the decay regime
#' (\eqn{d > d_0}) without introducing discontinuities that may hinder
#' numerical optimization.
#'
#' When SNP density is high, \eqn{d_0} captures the short-distance LD
#' plateau explicitly and stabilizes estimation of \eqn{c}. In sparse
#' marker datasets, \eqn{d_0} typically converges to zero, reducing the
#' model to the standard hyperbolic decay formulation.
#'
#' @param gds An open GDS object.
#' @param idx Optional vector of SNP indices.
#' @param ld_struct Per chromosome edge list of LD-values produced e.g. by \code{compute_ld_structure}
#' @param b Optional background LD value. If missing, estimated automatically.
#' @param q Quantile used for LD summarization (default 0.95).
#' @param n_sub Number of SNPs used to estimate background LD.
#' @param slide_win_ld Window size passed to `snpgdsLDMat`.
#' @param window_size Genomic window size (bp).
#' @param step_size Step size (bp) between windows.
#' @param n_cores_ld Number of cores for LD computation.
#' @param dist_unit Distance bin width (bp).
#'
#' @return An object of class `"ld_decay"` containing:
#' \describe{
#'   \item{data}{Window-wise LD-decay parameters.}
#'   \item{summary}{Per-chromosome median parameters.}
#'   \item{b}{Background LD estimate.}
#' }
#'
#' @examples
#' \dontrun{
#' decay <- ld_decay(gds,ld_struct)
#' plot(decay)
#' }
#'
#' @export
ld_decay <- function(gds,
                     idx,
                     ld_struct,
                     q = 0.95,
                     n_sub = 5000,
                     slide_win_ld = 1000,
                     window_size = 1e6,
                     step_size = 5e5,
                     cores = 1,
                     dist_unit = 5000) {

  t1 <- Sys.time()

  ids <- .read_gds_ids(gds)

  if (missing(idx) || is.null(idx)) idx <- seq_along(ids$snp_id)

  message("Estimating background LD")
  b <- get_bg_ld(gds,idx,n_sub = n_sub, q = q)

  message("Background LD (b) =", sprintf("%.4f", b))

  decay_data <- ld_decay_by_chr_win(
    ld_struct=ld_struct,
    q = q,
    b = b,
    window_size = window_size,
    step_size = step_size,
    cores = cores
  )

  decay_summary <- decay_data[, .(
    Chr_size = max(end),
    c = median(c, na.rm=TRUE),
    a = median(a, na.rm=TRUE),
    d0 = median(d0, na.rm=TRUE),
    b = b
  ), by=Chr]

  t2 <- Sys.time()
  print(difftime(t2,t1))
  cat("\n")

  out <- list(
    data = decay_data,
    summary = decay_summary,
    b = b,
    q = q,
    window_size = window_size,
    step_size = step_size,
    call = match.call()
  )

  class(out) <- "ld_decay"
  out
}


ld_decay_by_chr_win <- function(ld_struct, q = 0.95, dist_unit = 5000, window_size=1e7,
                                step_size =5e+05,b = 0.05,cores=1) {

  rbindlist(parallel::mclapply(names(ld_struct$by_chr),function(ch){
    cat(ch,"..")

    el <- ld_struct$by_chr[[ch]]$edges

    # get windows for ld-decay analyses
    min_pos <- min(el$pos1, el$pos2)
    max_pos <- max(el$pos1, el$pos2)
    starts <- seq(min_pos, max_pos - window_size, by = step_size)
    ends <- starts + window_size

    #i <- 1
    # el_ld <- sub
    rbindlist(lapply(seq_along(starts),function(i) {

      sub <- el[pos1 >= starts[i] & pos1 < ends[i] & pos2 >= starts[i] & pos2 < ends[i]]

      coefs <- tryCatch({
        coef_ld_dec(sub,q = 0.95, dist_unit = 5000, b = b)
      }, error = function(e) NULL)

      if(is.null(coefs)){
        data.table(Chr = ch, start=starts[i],end=ends[i],c = NA, a = NA, d0=NA)
      }else{
        data.table(Chr = ch, start=starts[i],end=ends[i],t(coefs))
      }
    }))

  },mc.cores=cores), use.names = TRUE, fill = TRUE)

}

# el_ld <- sub

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

  ## ---- starting values ----
  starts <- list(
    c  = min(1, bin_dt$r2_q[1]),  # avoid starting >1
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
      lower = c(b, 0, 0),
      upper = c(1, Inf, max(bin_dt$dist_mean)),
      control = nls.control(warnOnly = TRUE)
    )

  }, error = function(e) NULL)

  if (is.null(fit))
    return(c(c = NA, a = NA, d0 = NA))

  coef(fit)
}



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


d_from_rho <- function(a, rho, d0 = 0) {
  d0 + (1 / a) * (1 / (1 - rho) - 1)
}

ld_from_rho <- function(b, c = 1, rho){
  b + (c - b) * (1 - rho)
}

# ---- S3 methods ---- #
#' @export
print.ld_decay <- function(x, ...) {

  cat("\nLD-decay object\n")
  cat("---------------\n")

  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }

  cat("Background LD (b): ", sprintf("%.4f", x$b), "\n")
  cat("Quantile (q):      ", x$q, "\n")
  cat("Window size:       ", format(x$window_size, scientific = FALSE), "bp\n")
  cat("Step size:         ", format(x$step_size, scientific = FALSE), "bp\n\n")

  n_chr <- length(unique(x$summary$Chr))
  n_win <- nrow(x$data)

  cat("Chromosomes: ", n_chr, "\n")
  cat("Total windows:", n_win, "\n\n")

  cat("Genome-wide medians:\n")
  cat("  median(c): ", signif(median(x$data$c, na.rm=TRUE), 4), "\n")
  cat("  median(d0): ", signif(median(x$data$d0, na.rm=TRUE), 4), "\n")
  cat("  median(a): ", signif(median(x$data$a, na.rm=TRUE), 4), "\n\n")

  invisible(x)
}

#' Plot LD-decay curves
#'
#' Visualizes fitted LD-decay curves based on per-chromosome
#' median parameters.
#'
#' @param x An object of class `"ld_decay"`.
#' @param chr Optional chromosome to plot.
#' @param max_dist Maximum distance (bp) shown on x-axis.
#' @param n_points Number of points used to draw curve.
#' @param ... Additional arguments.
#'
#' @export

plot.ld_decay <- function(x,
                          chr = NULL,
                          max_dist = NULL,
                          n_points = 500,
                          ...) {

  summ <- x$summary

  if (!is.null(chr)) {
    summ <- summ[Chr == chr]
    if (nrow(summ) == 0)
      stop("Chromosome not found.")
  }

  if (is.null(max_dist)) {
    max_dist <- x$window_size
  }

  d <- seq(0, max_dist, length.out = n_points)

  plot(NULL,
       xlim = c(0, max_dist),
       ylim = c(0, 1),
       xlab = "Distance (bp)",
       ylab = expression(r^2),
       main = if (is.null(chr))
         "LD-decay curve"
       else
         paste("LD-decay curve -", chr))

  ## Background LD line
  abline(h = x$b, lty = 2, col = "grey40")

  for (i in seq_len(nrow(summ))) {

    a      <- summ$a[i]
    c_par  <- summ$c[i]
    b      <- summ$b[i]
    d0     <- if ("d0" %in% names(summ)) summ$d0[i] else 0

    ## Shifted plateau model
    r2 <- b + (c_par - b) / (1 + a * pmax(d - d0, 0))

    lines(d, r2,
          lwd = 2,
          col = scales::alpha("steelblue", 0.5))
  }

  invisible(x)
}

