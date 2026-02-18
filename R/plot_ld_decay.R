#' Plot empirical LD structure and fitted LD-decay curve
#'
#' Visualises pairwise LD values together with the fitted LD-decay model
#' for a selected chromosome. This diagnostic plot overlays empirical
#' \eqn{r^2} values with the estimated decay curve:
#'
#' \deqn{
#' r^2(d) = b + \frac{c - b}{1 + a \max(d - d_0, 0)},
#' }
#'
#' where \eqn{b} is background LD, \eqn{c} is the short-range plateau,
#' \eqn{a} controls the decay rate, and \eqn{d_0} represents the plateau
#' distance (if estimated).
#'
#' The function subsamples LD pairs when necessary to reduce overplotting.
#'
#' @param x An object returned by \code{\link{ld_pipeline}} containing
#'   components \code{ld_struct} and \code{decay}.
#' @param chr Character string. Chromosome to visualise. If \code{NULL},
#'   the first chromosome in \code{x$ld_struct} is used.
#' @param n_points Integer. Maximum number of LD pairs to display.
#'   If the number of available LD pairs exceeds this value, a random
#'   subset is drawn for plotting.
#' @param n_curve Integer. Number of points used to draw the fitted
#'   decay curve.
#' @param alpha Numeric in (0,1). Transparency level for plotted LD points.
#'
#' @details
#' This function is intended as a diagnostic tool to evaluate goodness-of-fit
#' of the LD-decay model. When a plateau parameter \eqn{d_0} is present,
#' LD remains approximately constant up to distance \eqn{d_0}, after which
#' decay follows a hyperbolic form.
#'
#' For large datasets, LD structures can contain millions of pairs.
#' The \code{n_points} argument prevents excessive overplotting and
#' reduces rendering time.
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \dontrun{
#' pipe <- ld_pipeline(geno, map, F_cols = c("emx_F", "lfmm_F"))
#' plot_ld_decay(pipe, chr = "Chr1")
#' }
#'
#' @export
plot_ld_decay <- function(x,
                       chr = NULL,
                       n_points = 5000,
                       n_curve = 500,
                       max_dist = NULL,
                       alpha = 0.75) {

  if (is.null(chr))
    chr <- names(x$ld_struct$by_chr)[1]

  if (!chr %in% names(x$ld_struct$by_chr))
    stop("Chromosome not found in LD structure.")

  ## Extract LD edges
  el <- x$ld_struct$by_chr[[chr]]$edges

  ## Distance grid for curve
  if(is.null(max_dist)){
    max_dist <- max(el$d, na.rm = TRUE)
  }

  el <- el[d<max_dist]

  if (nrow(el) == 0)
    stop("No LD edges available for this chromosome.")

  ## Subsample for plotting
  if (nrow(el) > n_points) {
    set.seed(1)
    el <- el[sample(.N, n_points)]
  }

  ## Extract decay parameters
  summ <- x$decay$summary[x$decay$summary$Chr == chr, ]
  if (nrow(summ) == 0)
    stop("Chromosome not found in decay summary.")

  a  <- summ$a
  c  <- summ$c
  b  <- summ$b
  d0 <- if ("d0" %in% names(summ)) summ$d0 else 0

  d_seq <- seq(0, max_dist, length.out = n_curve)

  ## Decay curve with plateau
  r2_curve <- b + (c - b) / (1 + a * pmax(d_seq - d0, 0))

  ## ---- Plot ----
  plot(el$d,
       el$r2,
       pch = 16,
       cex = 0.5,
       col = grDevices::adjustcolor("grey40", alpha.f = alpha),
       xlab = "Distance (bp)",
       ylab = expression(r^2),
       main = paste("LD structure + fitted decay —", chr))

  ## Decay curve
  lines(d_seq, r2_curve,
        col = "steelblue",
        lwd = 3)

  ## Background LD
  abline(h = b,
         col = "firebrick",
         lty = 2,
         lwd = 2)

  ## Plateau marker
  if (!is.na(d0) && d0 > 0) {
    abline(v = d0,
           col = "darkgreen",
           lty = 3,
           lwd = 2)
  }

  legend(x=max_dist*0.6,
         y=1,
         legend = c("LD pairs",
                    "Fitted decay",
                    "Background LD",
                    if (d0 > 0) "Plateau distance d0"),
         col = c("grey40",
                 "steelblue",
                 "firebrick",
                 if (d0 > 0) "darkgreen"),
         lty = c(NA, 1, 2, if (d0 > 0) 3),
         pch = c(16, NA, NA, NA),
         bty = "n")

  invisible(x)
}
