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
ld_rho_draws <- function(gds,
                         ld_struct,
                         n_inds,
                         F_vals,
                         q_vals=NULL,
                         n_rho = 25,
                         n_or_draws = 25,
                         rho_w_lim = list(min=0.8,max=0.99),
                         rho_d_lim=list(min=0.5,max=0.999),
                         rho_ld_lim=list(min=0.9,max=0.999),
                         alpha_lim=list(min=1.31,max=4),
                         lmin_lim=list(min=1,max=10),
                         cores=1,
                         seed = NULL,
                         mode = c("joint","per_method"),
                         use = c("robust","median")
                         ){

  if (!is.null(seed))
    set.seed(seed)

  rho_values <- runif(n_rho, rho_w_lim$min, rho_w_lim$max)
  ids <- .read_gds_ids(gds)

  # i <- 1
  # rho_w <- 0.9
  run_one <- function(i) {

    cat("rho_w draw", i, "..\n")

    scan <- ld_scan(
      ld_struct = ld_struct,
      SNP_ids   = ids$snp_id,
      F_vals    = F_vals,
      rho_w     = rho_values[i],
      n_inds    = n_inds,
      full      = TRUE,
      use       = ld_struct$params$use
    )

    q_primes <- do.call(cbind,lapply(scan$result,function(x) x$q_prime))
    colnames(q_primes) <- paste0(colnames(q_primes),"_prime")
    qvals <- cbind(q_vals, q_primes)
    q_min <- apply(qvals,1,min)
    #plot(-log10(q_min))
    ## pre-estimate all pairwise LD values for outliers
    idx <- which(q_min<1/10^alpha_lim$min)
    el  <- get_el(gds, idx, slide_win_ld = -1,by_chr = TRUE)


    out <- or_draws(
      el         = el,
      q_vals     = qvals,
      SNP_ids    = ids$snp_id,
      SNP_chr    = ids$snp_chr,
      ld_struct  = ld_struct,
      n_draws    = n_or_draws,
      rho_d_lim  = rho_d_lim,
      rho_ld_lim = rho_ld_lim,
      alpha_lim  = alpha_lim,
      lmin_lim   = lmin_lim,
      mode       = mode[1]
    )
    #mode="joint"

    cbind(rho_w=rho_values[i],out)
  }

  if (.Platform$OS.type == "unix") {
    draws <- rbindlist(parallel::mclapply(
      seq_len(n_rho),
      run_one,
      mc.cores = cores
    ))
  } else {
    warning("Parallel rho_w draws not supported on Windows; using single core.")
    draws <- rbindlist(lapply(seq_len(n_rho), run_one))
  }

  structure(
    list(draws=draws,
         n_rho=n_rho,
         n_or_draws=n_or_draws),
    class = "ld_rho_draws"
  )
}



#' @export
print.ld_rho_draws <- function(x, ...) {
  methods <- unique(x$draws$method)
  cat("\nLD rho-window draws\n")
  cat("--------------------\n")
  cat("Number of rho_w draws:", x$n_rho, "\n")
  cat("OR draws per rho_w:", x$n_or_draws, "\n")

  cat("Included methods:", paste(methods,collapse="/"), "\n")
  invisible(x)
}
