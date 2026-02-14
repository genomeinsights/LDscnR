#' Perform multiple LD-window (rho_w) draws with OR detection
#'
#' Repeats LD-scaling across multiple LD-decay quantiles (rho_w),
#' and for each run performs multiple OR parameter draws.
#'
#' @param ld_struct Object of class "ld_structure".
#' @param decay_obj Object of class "ld_decay".
#' @param SNP_ids Vector of SNP IDs in correct order.
#' @param n_inds Number of individuals.
#' @param F_vals Matrix or data.frame of F statistics.
#' @param q_vals Matrix or data.frame of q values (optional).
#' @param n_rho Number of LD-window draws.
#' @param rho_w_lim List(min,max) for rho_w.
#' @param n_or_draws Number of OR parameter draws per rho_w.
#' @param rho_d_lim List(min,max) for rho_d.
#' @param rho_ld_lim List(min,max) for rho_ld.
#' @param alpha_lim List(min,max) for alpha exponent.
#' @param lmin_lim List(min,max) for minimum OR size.
#' @param seed Optional random seed.
#'
#' @return Object of class "ld_rho_draws".
#' @export
ld_rho_draws <- function(ld_struct,
                         decay_obj,
                         SNP_ids,
                         n_inds,
                         F_vals,
                         q_vals=NULL,
                         n_rho = 25,
                         rho_w_lim = list(min=0.8,max=0.99),
                         n_or_draws = 25,
                         rho_d_lim=list(min=0.5,max=0.999),
                         rho_ld_lim=list(min=0.9,max=0.999),
                         alpha_lim=list(min=1.31,max=4),
                         lmin_lim=list(min=1,max=10),
                         n_cores=8,
                         seed = NULL) {

  if (!is.null(seed))
    set.seed(seed)

  rho_values <- runif(n_rho, rho_w_lim$min, rho_w_lim$max)

  run_one <- function(i) {

    cat("rho_w draw", i, "..\n")

    scan <- ld_scan(
      ld_struct = ld_struct,
      decay_obj = decay_obj,
      SNP_ids   = SNP_ids,
      F_vals    = F_vals,
      rho_w     = rho_values[i],
      n_inds    = n_inds,
      full      = FALSE
    )

    q_primes <- do.call(cbind,lapply(scan$result,function(x) x$q_prime))
    colnames(q_primes) <- paste0(colnames(q_primes),"_prime")
    qvals <- cbind(q_vals, q_primes)

    or_draws(
      q_vals     = qvals,
      ld_struct  = ld_struct,
      decay_obj  = decay_obj,
      n_draws    = n_or_draws,
      rho_d_lim  = rho_d_lim,
      rho_ld_lim = rho_ld_lim,
      alpha_lim  = alpha_lim,
      lmin_lim   = lmin_lim
    )
  }

  if (.Platform$OS.type == "unix") {
    draws <- parallel::mclapply(
      seq_len(n_rho),
      run_one,
      mc.cores = n_cores
    )
  } else {
    warning("Parallel rho_w draws not supported on Windows; using single core.")
    draws <- lapply(seq_len(n_rho), run_one)
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
  methods <- names((x$draws[[1]]$draws[[1]]$ORs))
  cat("\nLD rho-window draws\n")
  cat("--------------------\n")
  cat("Number of rho_w draws:", x$n_rho, "\n")
  cat("OR draws per rho_w:", x$n_or_draws, "\n")
  cat("rho_w range:",
      paste0("[", x$rho_w_lim$min, ", ", x$rho_w_lim$max, "]"), "\n\n")

  cat("Included methods:", paste(methods,collapse="/"), "\n")
  invisible(x)
}
