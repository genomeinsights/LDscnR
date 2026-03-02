#' Compute LD-based summary statistics from histogram representation
#'
#' @param ld_structure Object of class "ld_structure"
#' @param method "ld_int" (kernel integration) or "ld_w" (hard truncation)
#' @param eps Tail tolerance for LD_int
#' @param d_window Physical cutoff for ld_w (if NULL and method="ld_w",
#'        it is derived from eps)
#' @param rho distance from rho; d_window = rho / (a * (1 - rho))
#' @param cores Number of cores
#'
#' @export
compute_ld_summary <- function(
    ld_structure,
    method = c("ld_int", "ld_w"),
    eps = 0.005,
    d_window = NULL,
    rho = 0.95
) {

  method <- match.arg(method)

  results <- parallel_apply(ld_structure$by_chr, function(chr_obj) {

    hist_list <- chr_obj$hist_obj
    a <- chr_obj$decay_sum$a

    if (method == "ld_int") {

      d_star <- d_from_rho(a, 1-rho)

      return(sapply(hist_list, function(h)
        integrate_ld_int_from_hist(h, a, d_star)
      ))
    }

    if (method == "ld_w") {

      if (is.null(d_window)) {
        d_window <- d_from_rho(a, rho)
      }

      return(sapply(hist_list, function(h)
        compute_ldw_interp_from_hist(h, d_window)
      ))
    }

  }, cores = cores)

  unlist(results)
}
compute_ld_summary <- function(
    ld_structure,
    method = c("ld_int", "ld_w"),
    rho,
    cores = 1
) {

  method <- match.arg(method)

  results <- parallel_apply(ld_structure$by_chr, function(chr_obj) {

    hist_list <- chr_obj$hist_obj

    sapply(hist_list, function(hist_mat) {

      if (is.null(hist_mat)) return(NA_real_)

      rho_vals <- as.numeric(rownames(hist_mat))
      keep <- rho_vals <= rho

      if (!any(keep)) return(NA_real_)

      h_sub <- hist_mat[keep, , drop = FALSE]

      # Collapse r² distributions
      if (method == "ld_w") {

        # Hard truncation
        r2_counts <- colSums(h_sub)
        return(weighted_median_from_counts(r2_counts))

      } else {

        # Uniform kernel in rho-space
        rho_sub <- rho_vals[keep]

        if (length(rho_sub) > 1) {
          delta_rho <- c(diff(rho_sub), tail(diff(rho_sub), 1))
        } else {
          delta_rho <- 1
        }

        w_shell <- delta_rho
        w_shell <- w_shell / sum(w_shell)

        # Weight counts per shell before collapsing
        weighted_counts <- sweep(h_sub, 1, w_shell, FUN = "*")

        r2_counts <- colSums(weighted_counts)

        weighted_median_from_counts(r2_counts)
      }
    })

  }, cores = cores)

  unlist(results)
}

weighted_median_from_counts <- function(counts) {

  total <- sum(counts)
  if (total == 0) return(NA_real_)

  r2_vals <- as.numeric(names(counts))
  cs <- cumsum(counts)

  r2_vals[which(cs >= total / 2)[1]]
}

compute_ld_int_from_shells <- function(
    shell_list,
    a
) {

  weighted_median <- function(x, w) {
    if (length(x) == 0) return(NA_real_)
    ord <- order(x)
    x <- x[ord]
    w <- w[ord]
    w <- w / sum(w)
    cs <- cumsum(w)
    x[which(cs >= 0.5)[1]]
  }

  sapply(shell_list, function(shell_dt) {

    if (is.null(shell_dt) || nrow(shell_dt) == 0)
      return(NA_real_)

    # Kernel weights
    w <- a / (1 + a * shell_dt$d_mid)^2
    w <- w * shell_dt$delta_d

    if (sum(w) == 0) return(NA_real_)

    weighted_median(shell_dt$median_r2, w)
  })
}

###### rests #####
integrate_ld_int_from_hist <- function(hist_mat, a, d_star) {

  if (is.null(hist_mat)) return(NA_real_)

  dist_vals <- as.numeric(rownames(hist_mat))
  keep <- dist_vals <= d_star

  if (!any(keep)) return(NA_real_)

  h_sub <- hist_mat[keep, , drop = FALSE]
  dist_vals_sub <- dist_vals[keep]

  # Kernel weights (derivative of decay curve)
  w_d <- a / (1 + a * dist_vals_sub)^2

  # Distance increments
  if (length(dist_vals_sub) > 1) {
    d_mid <- sqrt(breaks[-1] * breaks[-length(breaks)])
    delta_d <- diff(breaks)
  } else {
    delta_d <- 1
  }

  d_mid <- sqrt(breaks[-1] * breaks[-length(breaks)])
  delta_d <- diff(breaks)
  if (sum(w_use) == 0) return(NA_real_)

  w_use <- w_use / sum(w_use)

  # Median per distance shell
  shell_medians <- apply(h_sub, 1, median_from_counts)

  weighted_median(shell_medians, w_use)
}

compute_ldw_interp_from_hist <- function(hist_mat, d_window) {

  if (is.null(hist_mat)) return(NA_real_)

  dist_vals <- as.numeric(rownames(hist_mat))
  if (length(dist_vals) < 1) return(NA_real_)

  # Identify bin
  k <- findInterval(d_window, dist_vals)

  if (k == 0) return(NA_real_)

  # Estimate bin width
  if (length(dist_vals) > 1) {
    dist_unit <- dist_vals[2] - dist_vals[1]
  } else {
    dist_unit <- d_window
  }

  # Fully included bins
  if (k > 1) {
    h_full <- hist_mat[seq_len(k - 1), , drop = FALSE]
    full_counts <- colSums(h_full)
  } else {
    full_counts <- rep(0, ncol(hist_mat))
  }

  # Partial bin
  h_partial <- hist_mat[k, , drop = FALSE]
  d_bin_start <- dist_vals[k]

  frac <- (d_window - d_bin_start) / dist_unit
  frac <- max(min(frac, 1), 0)

  partial_counts <- frac * as.numeric(h_partial)

  r2_counts <- full_counts + partial_counts

  median_from_counts(r2_counts)
}

weighted_median_from_counts <- function(counts) {

  total <- sum(counts)
  if (total == 0) return(NA_real_)

  r2_vals <- as.numeric(names(counts))
  cs <- cumsum(counts)

  r2_vals[which(cs >= total / 2)[1]]
}

median_from_counts <- function(counts_row) {

  total <- sum(counts_row)
  if (total == 0) return(NA_real_)

  r2_vals <- as.numeric(names(counts_row))
  cs <- cumsum(counts_row)

  r2_vals[which(cs >= total / 2)[1]]
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
weighted_median <- function(x, w) {

  ord <- order(x)
  x <- x[ord]
  w <- w[ord]

  w <- w / sum(w)
  cw <- cumsum(w)

  x[which(cw >= 0.5)[1]]
}

derive_ld_radius <- function(a, eps) {
  (1 / a) * (1 / eps - 1)
}

#### rests
# summarize_hist_shell <- function(hist_row) {
#
#   type <- match.arg(type)
#
#   total <- sum(hist_row)
#   if (total == 0) return(NA_real_)
#
#   r2_vals <- as.numeric(names(hist_row))
#
#   cs <- cumsum(hist_row)
#   idx <- which(cs >= total / 2)[1]
#   return(r2_vals[idx])
#
# }
#
#
# integrate_ld_kernel <- function(hist_mat,
#                                 a,
#                                 b,
#                                 d_star) {
#
#   if (is.null(hist_mat)) return(NA_real_)
#
#   dist_vals <- as.numeric(rownames(hist_mat))
#   keep <- dist_vals <= d_star
#
#   if (!any(keep)) return(NA_real_)
#
#   h_sub <- hist_mat[keep, , drop = FALSE]
#   dist_vals_sub <- as.numeric(rownames(h_sub))
#
#   # Kernel weights
#   w_d <- a / (1 + a * dist_vals_sub)^2
#
#   if (length(dist_vals_sub) > 1) {
#     delta_d <- c(diff(dist_vals_sub),
#                  tail(diff(dist_vals_sub), 1))
#   } else {
#     delta_d <- 1
#   }
#
#   w_use <- w_d * delta_d
#   if (sum(w_use) == 0) return(NA_real_)
#   w_use <- w_use / sum(w_use)
#
#   # Shell summaries
#   shell_vals <- apply(h_sub, 1, summarize_hist_shell,
#                       b = b)
#
#   sum(shell_vals * w_use)
# }
#
#
# integrate_ld_kernel_median <- function(hist_mat,
#                                        a,
#                                        b,
#                                        d_star) {
#
#   if (is.null(hist_mat)) return(NA_real_)
#
#   dist_vals <- as.numeric(rownames(hist_mat))
#   keep <- dist_vals <= d_star
#
#   if (!any(keep)) return(NA_real_)
#
#   h_sub <- hist_mat[keep, , drop = FALSE]
#   dist_vals_sub <- as.numeric(rownames(h_sub))
#
#   # decay weights
#   w_d <- a / (1 + a * dist_vals_sub)^2
#
#   if (length(dist_vals_sub) > 1) {
#     delta_d <- c(diff(dist_vals_sub),
#                  tail(diff(dist_vals_sub), 1))
#   } else {
#     delta_d <- 1
#   }
#
#   w_use <- w_d * delta_d
#   if (sum(w_use) == 0) return(NA_real_)
#
#   # shell summaries
#   shell_vals <- apply(h_sub, 1, summarize_hist_shell,
#                       b = b)
#
#   weighted_median(shell_vals, w_use)
# }
#
# #d_window
# compute_ld_rhow_from_hist <- function(hist_mat,
#                                       d_window) {
#
#   if (is.null(hist_mat)) return(NA_real_)
#
#   dist_vals <- as.numeric(rownames(hist_mat))
#   keep <- dist_vals <= d_window
#
#   if (!any(keep)) return(NA_real_)
#
#   h_sub <- hist_mat[keep, , drop = FALSE]
#
#   # Collapse all rows into one distribution
#   combined <- colSums(h_sub)
#
#   summarize_hist_shell(combined)
# }
#
#
# uniform_thin_to_distance <- function(pos,
#                                      d_star,
#                                      k_max,
#                                      tol = 1e-8) {
#
#   pos <- sort(pos)
#
#   # required spacing
#   target_spacing <- d_star / k_max
#
#   # current median spacing
#   med_spacing <- median(diff(pos))
#
#   # if already sufficient, no thinning needed
#   if (k_max * med_spacing >= d_star) {
#     return(seq_along(pos))
#   }
#
#   keep <- integer()
#   last_kept_pos <- -Inf
#
#   for (i in seq_along(pos)) {
#     if ((pos[i] - last_kept_pos) >= (target_spacing - tol)) {
#       keep <- c(keep, i)
#       last_kept_pos <- pos[i]
#     }
#   }
#
#   keep
# }
#
#
# compute_shell_medians <- function(dt,
#                                   d_star,
#                                   n_dist_target) {
#
#   if (nrow(dt) == 0) return(NULL)
#
#   dt <- dt[d <= d_star]
#   if (nrow(dt) == 0) return(NULL)
#
#   dist_unit <- d_star / n_dist_target
#   dt[, dist_idx := floor(d / dist_unit)]
#
#   shell_dt <- dt[, .(
#     median_r2 = median(r2)
#   ), by = dist_idx]
#
#   shell_dt
# }
#
# integrate_ld_kernel_median_shell <- function(shell_dt,
#                                              a,
#                                              d_star) {
#
#   if (is.null(shell_dt) || nrow(shell_dt) == 0)
#     return(NA_real_)
#
#   dist_unit <- d_star / max(shell_dt$dist_idx)
#
#   dist_vals <- shell_dt$dist_idx * dist_unit
#
#   w_d <- a / (1 + a * dist_vals)^2
#
#   if (length(dist_vals) > 1) {
#     delta_d <- c(diff(dist_vals), tail(diff(dist_vals), 1))
#   } else {
#     delta_d <- 1
#   }
#
#   w_use <- w_d * delta_d
#   if (sum(w_use) == 0) return(NA_real_)
#   w_use <- w_use / sum(w_use)
#
#   weighted_median(shell_dt$median_r2, w_use)
# }
#
# compute_ld_int_kernel_smooth <- function(dt, a, d_star, n_bins = 20) {
#
#   dt <- dt[d <= d_star]
#   if (nrow(dt) == 0) return(NA_real_)
#
#   # sort by distance
#   dt <- dt[order(d)]
#
#   # equal-width bins
#   breaks <- seq(0, d_star, length.out = n_bins + 1)
#   dt[, bin := cut(d, breaks = breaks, include.lowest = TRUE, labels = FALSE)]
#
#   shell_dt <- dt[, .(
#     median_r2 = median(r2)
#   ), by = bin]
#
#   # representative distance per bin
#   dist_vals <- breaks[shell_dt$bin]
#
#   # kernel weights
#   w_d <- a / (1 + a * dist_vals)^2
#
#   delta_d <- diff(breaks)
#   w_use <- w_d * delta_d[shell_dt$bin]
#   w_use <- w_use / sum(w_use)
#
#   weighted_median(shell_dt$median_r2, w_use)
# }
#
# compute_ldw_interp <- function(hist_mat, d_window) {
#
#   if (is.null(hist_mat)) return(NA_real_)
#
#   dist_vals <- as.numeric(rownames(hist_mat))
#   dist_unit <- dist_vals[2] - dist_vals[1]
#
#   # Identify bin index
#   k <- findInterval(d_window, dist_vals)
#
#   if (k == 0) return(NA_real_)
#
#   # Full bins
#   h_full <- hist_mat[seq_len(k-1), , drop = FALSE]
#
#   # Partial bin
#   h_partial <- hist_mat[k, , drop = FALSE]
#
#   d_bin_start <- dist_vals[k]
#   frac <- (d_window - d_bin_start) / dist_unit
#   frac <- max(min(frac, 1), 0)
#
#   # Combine
#   r2_counts <- colSums(h_full) + frac * as.numeric(h_partial)
#
#   total <- sum(r2_counts)
#   if (total == 0) return(NA_real_)
#
#   r2_vals <- as.numeric(colnames(hist_mat))
#   cs <- cumsum(r2_counts)
#
#   r2_vals[which(cs >= total/2)[1]]
# }
