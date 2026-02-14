#' LD-scaled genome scan
#'
#' Computes LD-scaled F statistics (F') for a given LD-decay object
#' and LD-structure, using a decay quantile \eqn{\rho_w}.
#'
#' @param ld_struct An object of class `"ld_structure"`.
#' @param decay_obj An object of class `"ld_decay"`.
#' @param F_vals Numeric vector or matrix of F-statistics (n_snp x n_method).
#' @param rho_w Numeric in (0,1). LD-decay quantile defining window size.
#' @param n_inds Number of individuals (for F null df2 = n_inds - 2).
#' @param n_rep Number of circular-shift permutations.
#' @param full Logical; if TRUE return QQ diagnostics.
#'
#' @return An object of class `"ld_scan"`.
#' @export
ld_scan <- function(ld_struct,
                    decay_obj,
                    SNP_ids,
                    F_vals,
                    rho_w,
                    n_inds,
                    n_rep = 10,
                    full = TRUE) {

  if (!inherits(ld_struct, "ld_structure"))
    stop("`ld_struct` must be of class 'ld_structure'.")

  if (!inherits(decay_obj, "ld_decay"))
    stop("`decay_obj` must be of class 'ld_decay'.")

  if (!is.finite(rho_w) || rho_w <= 0 || rho_w >= 1)
    stop("`rho_w` must lie in (0,1).")

  ## 1️⃣ Compute ld_w from structure
  ld_w <- compute_ld_w(ld_struct, decay_obj, rho_w)

  ## 2️⃣ Compute F' using existing logic
  scan_res <- .compute_Fprime(
    F_vals = F_vals,
    ld_w   = ld_w,
    n_rep  = n_rep,
    n_inds = n_inds,
    full   = full
  )

  ## 3️⃣ Wrap into S3 object
  out <- list(
    rho_w   = rho_w,
    F_vals  = F_vals,
    n_inds  = n_inds,
    SNP_ids = SNP_ids,
    ld_w    = ld_w,
    result  = scan_res,
    call    = match.call()
  )

  class(out) <- "ld_scan"
  out
}

.compute_Fprime <- function(F_vals, ld_w, n_rep, n_inds, full) {

  F_mat <- as.matrix(F_vals)
  n <- nrow(F_mat)

  if (length(ld_w) != n)
    stop("Length of `ld_w` must match nrow(F_vals).")

  df2 <- as.integer(n_inds - 2)

  ld_w <- as.numeric(ld_w)

  take_sorted_quantiles <- function(sorted_x, n_out) {
    m <- length(sorted_x)
    if (m == 0L) return(rep(NA_real_, n_out))
    pos <- round(seq(1, m, length.out = n_out))
    sorted_x[pos]
  }
  #j = 1
  res <- lapply(seq_len(ncol(F_mat)), function(j) {

    Fval <- F_mat[, j]
    F_mean <- mean(Fval, na.rm = TRUE)

    F_prime_obs <- Fval * ld_w
    denom <- mean(F_prime_obs, na.rm = TRUE)

    if (!is.finite(denom) || denom == 0) {
      F_prime_obs <- rep(F_mean, length(Fval))
    } else {
      F_prime_obs <- F_prime_obs / denom * F_mean
    }

    ord <- order(F_prime_obs)
    F_prime_obs_sorted <- F_prime_obs[ord]

    null_all <- numeric(n * n_rep)
    #k=1
    for (k in seq_len(n_rep)) {
      start <- sample.int(n, 1)
      shifted <- c(ld_w[start:n], ld_w[1:(start - 1)])

      null_k <- Fval * shifted
      denom_k <- mean(null_k, na.rm = TRUE)

      if (!is.finite(denom_k) || denom_k == 0) {
        null_k <- rep(F_mean, n)
      } else {
        null_k <- null_k / denom_k * F_mean
      }

      null_all[((k - 1) * n + 1):(k * n)] <- null_k
    }

    null_all <- null_all[is.finite(null_all)]
    null_perm_sorted <- sort(null_all)
    null_perm <- take_sorted_quantiles(null_perm_sorted, n_out = n)

    null_true <- stats::qf(
      rev(stats::ppoints(n)),
      df1 = 1,
      df2 = df2,
      lower.tail = FALSE
    )

    q_F_prime_cor <- F_prime_obs_sorted - (null_perm - null_true)
    q_F_prime_cor_adj <- cummax(q_F_prime_cor)

    F_prime <- numeric(n)
    F_prime[ord] <- q_F_prime_cor_adj

    p_prime <- stats::pf(F_prime, df1 = 1, df2 = df2, lower.tail = FALSE)
    q_prime <- stats::p.adjust(p_prime, "fdr")

    if (!full) {
      return(list(
        F_prime = F_prime,
        p_prime = p_prime,
        q_prime = q_prime
      ))
    }

    qq_data <- data.table::data.table(
      null_true         = null_true,
      F_prime_obs       = F_prime_obs_sorted,
      F_obs             = sort(Fval),
      F_prime_cor       = q_F_prime_cor,
      F_prime           = q_F_prime_cor_adj,
      null_perm         = null_perm
    )

    list(
      qq_data = qq_data,
      F_prime = F_prime,
      p_prime = p_prime,
      q_prime = q_prime
    )
  })

  names(res) <- colnames(F_mat)
  res
}
#plot(res$emx_F$q_prime,res$lfmm_F$q_prime)

#' @export
compute_ld_structure <- function(gds,
                                 slide_win_ld = 1000,
                                 n_cores = 1) {

  ids <- .read_gds_ids(gds)
  chrs <- unique(ids$snp_chr)

  out <- vector("list", length(chrs))
  names(out) <- chrs

  for (ch in chrs) {
    chr_idx <- which(ids$snp_chr == ch)

    el <- get_el(gds,
                 idx = chr_idx,
                 slide_win_ld = slide_win_ld,
                 n_cores = n_cores)

    data.table::setorder(el, d)  # sort by distance once

    out[[ch]] <- list(
      snp_ids = ids$snp_id[chr_idx],
      edges   = el
    )
  }

  structure(
    list(
      by_chr = out
    ),
    class = "ld_structure"
  )
}

.unbiased_ld_w <- function(sub) {

  sub <- data.table::rbindlist(list(
    sub[, .(SNP = SNP1, OTHER = SNP2, pos = pos1, pos_other = pos2, r2)],
    sub[, .(SNP = SNP2, OTHER = SNP1, pos = pos2, pos_other = pos1, r2)]
  ))

  sub[, d := abs(pos_other - pos)]
  sub[, side := ifelse(pos_other > pos, "R", "L")]

  sub[, `:=`(
    maxL = if (any(side == "L")) max(d[side == "L"]) else 0,
    maxR = if (any(side == "R")) max(d[side == "R"]) else 0
  ), by = SNP]

  sub[, S := pmin(maxL, maxR)]

  sub[, long_side := data.table::fifelse(
    maxL > maxR, "L",
    data.table::fifelse(maxR > maxL, "R", NA_character_)
  )]

  sub[, tail := (side == long_side & d > S)]
  sub[is.na(tail), tail := FALSE]

  use <- data.table::rbindlist(list(
    sub[, .(SNP, r2)],
    sub[tail == TRUE, .(SNP, r2)]
  ))

  ld_w2 <- use[, .(
    ld_w = median(r2, na.rm = TRUE),
    ld_w_mean = mean(r2, na.rm = TRUE),
    N = .N
  ), by = SNP]

  return(ld_w2)
}


compute_ld_w <- function(ld_struct,
                         decay_obj,
                         rho_w) {

  result <- vector("list", length(ld_struct$by_chr))
  names(result) <- names(ld_struct$by_chr)

  for (ch in names(ld_struct$by_chr)) {

    a_chr <- decay_obj$summary[Chr == ch, a]
    d_th  <- d_from_rho(a_chr, rho_w)

    el <- ld_struct$by_chr[[ch]]$edges

    # Filter by distance
    sub <- el[d < d_th]

    if (nrow(sub) == 0) {
      result[[ch]] <- rep(NA_real_,
                          length(ld_struct$by_chr[[ch]]$snp_ids))
      next
    }

    ld_dt <- .unbiased_ld_w(sub)

    res <- data.table::data.table(
      SNP = ld_struct$by_chr[[ch]]$snp_ids
    )

    res <- ld_dt[res, on = "SNP"]

    result[[ch]] <- res$ld_w
  }

  unlist(result)
}


d_from_rho <- function(a, rho) {
  (1 / a) * (1 / (1 - rho) - 1)
}

#' Plot LD-scaled quantile transformation
#'
#' Visualizes the quantile transformation underlying the LD-scaled
#' statistic, decomposing LD-induced distortion and
#' selection-consistent excess.
#'
#' @param x An object of class `"ld_scan"`.
#' @param method Optional method name if multiple statistics were supplied.
#' @param ... Additional arguments.
#'
#' @return A ggplot object.
#' @export
#' #' Plot LD-scan results
#'
#' Produces a two-panel figure:
#' 1) Quantile decomposition (LD distortion + selection excess)
#' 2) Comparison of original, LD-weighted, and LD-scaled statistics
#'
#' @param x An object of class "ld_scan".
#' @param ... Unused.
#'
#' @return A patchwork ggplot object.
#' @export
#' @export

# plot.ld_scan <- function(x, ...) {
#
#   if (!inherits(x, "ld_scan"))
#     stop("Object must be of class 'ld_scan'.")
#
#   if (is.null(x$result[[1]]$qq_data))
#     stop("QQ data not stored (run ld_scan with full = TRUE).")
#
#   method_names <- names(x$result)
#   if (is.null(method_names))
#     method_names <- paste0("Method_", seq_along(x$result))
#
#   ############################################################
#   ## Stack QQ data across methods
#   ############################################################
#
#   dt_all <- data.table::rbindlist(
#     lapply(seq_along(x$result), function(i) {
#       dt <- data.table::copy(x$result[[i]]$qq_data)
#       dt[, method := method_names[i]]
#       dt
#     }),
#     use.names = TRUE
#   )
#
#   ############################################################
#   ## PANEL 1 — LD distortion + signal
#   ############################################################
#
#   dt_ribbon <- data.table::rbindlist(list(
#     dt_all[, .(
#       q = null_true,
#       ymin = null_true,
#       ymax = null_perm,
#       type = "LD-induced distortion",
#       method
#     )],
#     dt_all[, .(
#       q = null_true,
#       ymin = null_perm,
#       ymax = F_prime_obs,
#       type = "Selection-consistent excess",
#       method
#     )]
#   ))
#
#   dt_lines <- data.table::rbindlist(list(
#     dt_all[, .(q = null_true, value = null_true,
#                curve = "Theoretical null", method)],
#     dt_all[, .(q = null_true, value = null_perm,
#                curve = "Permutation null", method)],
#     dt_all[, .(q = null_true, value = F_prime_obs,
#                curve = "Raw LD-weighted statistic", method)]
#   ))
#
#   p1 <- ggplot2::ggplot() +
#
#     ggplot2::geom_ribbon(
#       data = dt_ribbon,
#       ggplot2::aes(q, ymin = ymin, ymax = ymax,
#                    fill = type, group = interaction(type, method)),
#       alpha = 0.25
#     ) +
#
#     ggplot2::geom_line(
#       data = dt_lines,
#       ggplot2::aes(q, value,
#                    colour = method,
#                    linetype = curve),
#       linewidth = 1
#     ) +
#
#     ggplot2::labs(
#       title = paste("Raw LD-weighted statistic vs permuted null:",
#                     paste(method_names, collapse = ", ")),
#       x = "Expected quantiles under null",
#       y = "Observed quantiles"
#     ) +
#     ggplot2::coord_equal() +
#     ggplot2::theme_bw(base_size = 11)
#
#   ############################################################
#   ## PANEL 2 — LD-scaled vs original
#   ############################################################
#
#   dt2 <- data.table::rbindlist(list(
#     dt_all[, .(q = null_true, value = F_obs,
#                curve = "Original F-value", method)],
#     dt_all[, .(q = null_true, value = F_prime,
#                curve = "LD-scaled statistic (F\u2032)", method)],
#     dt_all[, .(q = null_true, value = null_true,
#                curve = "Theoretical null", method)]
#   ))
#
#   dt_excess2 <- dt_all[
#     F_prime > null_true,
#     .(q = null_true,
#       ymin = null_true,
#       ymax = F_prime,
#       method)
#   ]
#
#   p2 <- ggplot2::ggplot() +
#
#     ggplot2::geom_ribbon(
#       data = dt_excess2,
#       ggplot2::aes(q, ymin = ymin, ymax = ymax,
#                    fill = method),
#       alpha = 0.15
#     ) +
#
#     ggplot2::geom_line(
#       data = dt2,
#       ggplot2::aes(q, value,
#                    colour = method,
#                    linetype = curve),
#       linewidth = 1
#     ) +
#
#     ggplot2::labs(
#       title = paste("LD-scaled statistic (F\u2032) vs original:",
#                     paste(method_names, collapse = ", ")),
#       x = "Expected quantiles under null",
#       y = "Observed quantiles"
#     ) +
#
#     ggplot2::coord_equal() +
#     ggplot2::theme_bw(base_size = 11)
#
#   ############################################################
#   ## Combine
#   ############################################################
#
#   combined <- (p1 | p2) +
#     patchwork::plot_layout(guides = "collect") &
#     ggplot2::theme(legend.position = "right")
#
#   return(combined)
# }
# method="lfmm_F"
plot.ld_scan <- function(x, method) {

  if (!inherits(x, "ld_scan"))
    stop("Object must be of class 'ld_scan'.")

  if (is.null(x$result[[1]]$qq_data))
    stop("QQ data not stored (run ld_scan with full = TRUE).")
  #res <- x$result[[1]]

  dt <- data.table::copy(x$result[[method]]$qq_data)

  # ----------------------------------------------------------
  # Thin lower quantiles to reduce overplotting
  # Keep:
  #   - all points above thinning threshold
  #   - every kth point below threshold
  # ----------------------------------------------------------

  thin_prop <- 0.7   # lower 40% thinned
  thin_step <- 100    # keep every 10th point

  n <- nrow(dt)
  cutoff <- floor(n * thin_prop)

  idx_keep <- c(
    seq(1, cutoff, by = thin_step),
    seq(cutoff + 1, n)
  )

  dt <- dt[idx_keep]
  required <- c("F_prime_obs", "null_perm", "F_obs")
  if (!all(required %in% colnames(dt)))
    stop("Required QQ components not found.")

  ############################################################
  ## Panel 1 — Quantile decomposition
  ############################################################

  dt_ld <- data.table::data.table(
    q = dt$null_true,
    ymin = pmin(dt$null_true, dt$null_perm),
    ymax = pmax(dt$null_true, dt$null_perm),
    type = "LD-induced distortion"
  )

  dt_signal <- data.table::data.table(
    q = dt$null_true,
    ymin = pmin(dt$null_perm, dt$F_prime_obs),
    ymax = pmax(dt$null_perm, dt$F_prime_obs),
    type = "Selection-consistent excess"
  )

  dt_excess2 <- data.table::data.table(
    q    = dt$null_true,
    ymin = dt$null_true,
    ymax = dt$F_prime,
    type = "Selection-consistent excess"
  )

  # Keep only positive excess
  dt_excess2 <- dt_excess2[ymax > ymin]
  dt_ribbon <- data.table::rbindlist(list(dt_ld, dt_signal))
  dt_ribbon[,type:=factor(type,levels=c("Selection-consistent excess","LD-induced distortion"))]

  dt_lines <- data.table::rbindlist(list(
    data.table::data.table(
      q = dt$null_true,
      value = dt$null_true,
      curve = "Theoretical null"
    ),
    data.table::data.table(
      q = dt$null_true,
      value = dt$null_perm,
      curve = "Permutation null"
    ),
    data.table::data.table(
      q = dt$null_true,
      value = dt$F_prime_obs,
      curve = "Raw LD-weighted statistic"
    )
  ))

  p1 <- ggplot2::ggplot() +

    ggplot2::geom_ribbon(
      data = dt_ribbon,
      ggplot2::aes(q, ymin = ymin, ymax = ymax, fill = type),
      alpha = 0.4
    ) +

    ggplot2::geom_line(
      data = dt_lines,
      ggplot2::aes(q, value, colour = curve),
      linewidth = 1
    ) +

    ggplot2::scale_fill_manual(
      values = c(
        "LD-induced distortion" = "salmon",
        "Selection-consistent excess" = "#009E73"
      ),
      name = NULL
    ) +

    ggplot2::scale_color_manual(
      values = c(
        "Theoretical null" = "grey30",
        "Permutation null" = "grey50",
        "Raw LD-weighted statistic" = "#009E73"
      ),
      guide=NULL,
      name = NULL
    ) +

    ggplot2::labs(
      x = "Expected quantiles under null",
      y = "Observed quantiles"
    ) +

    ggplot2::coord_equal() +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      aspect.ratio = 1,
      legend.position.inside = c(0.1,0.8),
      legend.position = "inside",
      panel.grid = ggplot2::element_blank()
    )

  ############################################################
  ## Panel 2 — Statistic comparison
  ############################################################

  dt2 <- data.table::rbindlist(list(
    data.table::data.table(
      q = dt$null_true,
      value = dt$F_obs,
      curve = "Original F-value"
    ),
    # data.table::data.table(
    #   q = dt$null_true,
    #   value = dt$F_prime_obs,
    #   curve = "Raw LD-weighted statistic"
    # ),
    data.table::data.table(
      q = dt$null_true,
      value = dt$F_prime,
      curve = "Raw LD-weighted statistic/\nLD-scaled statistic (F\u2032)"
    ),
    data.table::data.table(
      q = dt$null_true,
      value = dt$null_true,
      curve = "Theoretical null"
    ),
    data.table::data.table(
      q = dt$null_true,
      value = dt$null_true,
      curve = "Permutation null"
    )
  ))
  p2 <- ggplot2::ggplot() +

    # Selection-consistent excess ribbon
    ggplot2::geom_ribbon(
      data = dt_excess2,
      ggplot2::aes(q, ymin = ymin, ymax = ymax, fill = type),
      alpha = 0.4
    ) +

    # Curves
    ggplot2::geom_line(
      data = dt2,
      ggplot2::aes(q, value, colour = curve),
      linewidth = 1
    ) +

    ggplot2::scale_fill_manual(
      values = c(
        "Selection-consistent excess" = "#009E73"
      ),
      guide=NULL,
      name = NULL
    ) +

    ggplot2::scale_color_manual(
      values = c(
        "Original F-value" = "steelblue",
        "Raw LD-weighted statistic/\nLD-scaled statistic (F\u2032)" = "#009E73",
        "Theoretical null" = "grey30",
        "Permutation null" = "grey50"
      ),
      name = NULL
    ) +

    ggplot2::labs(
      x = "Expected quantiles under null",
      y = "Observed quantiles"
    ) +

    ggplot2::coord_equal() +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.grid = ggplot2::element_blank()
    )


  ############################################################
  ## Combine with shared legend
  ############################################################

  combined <- (p1+ggtitle("Raw LD-weighted statistic vs. permuted null") | p2+ggtitle("LD-scaled statistic (F\u2032) vs. original F-value")) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "right")
  #method <- "lfmm"
  combined <- combined + patchwork::plot_annotation(title = method)

  return(combined)
}


#' @export
print.ld_scan <- function(x, ...) {

  cat("\nLD-scaled scan result\n")
  cat("----------------------\n")

  if (!is.null(x$call))
    print(x$call)

  n <- nrow(x$F_vals)
  cat("\nNumber of SNPs:", n, "\n")

  if (is.null(x$n_inds)) {
    warning("n_inds not stored in object; cannot compute genomic inflation.")
    invisible(x)
    return(invisible(x))
  }

  df2 <- x$n_inds - 2
  f_median_null <- stats::qf(0.5, df1 = 1, df2 = df2)

  F_mat <- as.matrix(x$F_vals)
  method_names <- colnames(F_mat)
  if (is.null(method_names))
    method_names <- paste0("Method_", seq_len(ncol(F_mat)))

  for (i in seq_along(method_names)) {

    name_i <- method_names[i]
    F_obs  <- F_mat[, i]

    # Original lambda
    lambda_raw <- stats::median(F_obs, na.rm = TRUE) / f_median_null

    # LD-scaled lambda
    F_prime <- x$result[[i]]$F_prime
    lambda_scaled <- stats::median(F_prime, na.rm = TRUE) / f_median_null

    min_q <- min(x$result[[i]]$q_prime, na.rm = TRUE)

    cat("\n", name_i, "\n", sep = "")
    cat("  Min q-value:         ", signif(min_q, 4), "\n")
    cat("  λ (raw F):           ", round(lambda_raw, 3), "\n")
    cat("  λ (LD-scaled F′):    ", round(lambda_scaled, 3), "\n")
  }

  cat("\n")

  invisible(x)
}

# print.ld_scan <- function(x, ...) {
#   cat("\nLD-scaled scan result\n")
#   cat("----------------------\n")
#   if (!is.null(x$call)) print(x$call)
#
#   cat("\nNumber of SNPs:", nrow(x$F_vals), "\n")
#   for(i in seq_along(x$result)){
#     cat("Min q-value","for",names(x$F_vals)[i],":", signif(min(x$result[[i]]$q_prime, na.rm=TRUE), 4), "\n")
#   }
#
#   invisible(x)
# }
