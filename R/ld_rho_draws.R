#' Perform multiple LD-window (rho_w) draws with OR detection
#'
#' Repeats LD-scaling across multiple LD-decay quantiles (rho_w),
#' and for each run performs multiple OR parameter draws.
#'
#' @param ld_struct Object of class "ld_structure".
#' @param decay_obj Object of class "ld_decay".
#' @param SNP_ids Vector of SNP IDs in correct order.
#' @param F_vals Matrix or data.frame of F statistics.
#' @param n_rho Number of LD-window draws.
#' @param rho_w_lim List with elements `min` and `max`.
#' @param n_or_draws Number of OR parameter draws per rho_w.
#' @param rho_d_lim List(min,max) for rho_d.
#' @param rho_ld_lim List(min,max) for rho_ld.
#' @param alpha_lim List(min,max) for alpha exponent.
#' @param lmin_lim List(min,max) for minimum OR size.
#' @param n_inds Number of individuals.
#' @param seed Optional random seed.
#'
#' @return Object of class "ld_rho_draws".
#' @export
ld_rho_draws <- function(ld_struct,
                         decay_obj,
                         SNP_ids,
                         F_vals,
                         n_rho = 10,
                         rho_w_lim,
                         n_or_draws = 25,
                         rho_d_lim,
                         rho_ld_lim,
                         alpha_lim,
                         lmin_lim,
                         n_inds,
                         seed = NULL) {

  if (!is.null(seed))
    set.seed(seed)

  if (!is.list(rho_w_lim) || !all(c("min","max") %in% names(rho_w_lim)))
    stop("rho_w_lim must be a list with min and max.")

  if (rho_w_lim$min >= rho_w_lim$max)
    stop("rho_w_lim$min must be < rho_w_lim$max.")

  draws <- vector("list", n_rho)
  rho_values <- numeric(n_rho)

  for (i in seq_len(n_rho)) {

    cat("rho_w draw", i, "..")

    rho_w <- runif(1, rho_w_lim$min, rho_w_lim$max)
    rho_values[i] <- rho_w

    scan <- ld_scan(
      ld_struct = ld_struct,
      decay_obj = decay_obj,
      SNP_ids   = SNP_ids,
      F_vals    = F_vals,
      rho_w     = rho_w,
      n_inds    = n_inds,
      full      = FALSE
    )

    draws[[i]] <- or_draws(
      scan_obj   = scan,
      ld_struct  = ld_struct,
      decay_obj  = decay_obj,
      n_draws    = n_or_draws,
      rho_d_lim  = rho_d_lim,
      rho_ld_lim = rho_ld_lim,
      alpha_lim  = alpha_lim,
      lmin_lim   = lmin_lim
    )
  }

  structure(
    list(
      rho_w_values = rho_values,
      rho_w_lim    = rho_w_lim,
      draws        = draws,
      n_rho        = n_rho,
      n_or_draws   = n_or_draws
    ),
    class = "ld_rho_draws"
  )
}

#' @export
print.ld_rho_draws <- function(x, ...) {

  cat("\nLD rho-window draws\n")
  cat("--------------------\n")
  cat("Number of rho_w draws:", x$n_rho, "\n")
  cat("OR draws per rho_w:", x$n_or_draws, "\n")
  cat("rho_w range:",
      paste0("[", x$rho_w_lim$min, ", ", x$rho_w_lim$max, "]"), "\n\n")

  invisible(x)
}
