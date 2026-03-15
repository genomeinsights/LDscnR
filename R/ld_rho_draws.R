#' Perform multiple LD-window (rho_w) draws with OR detection
#'
#' Repeats LD-scaling across multiple LD-decay quantiles (rho_w),
#' and for each run performs multiple OR parameter draws.
#'
#' @param SNP_ids Vector of SNP IDs in correct order.
#' @param n_inds Number of individuals.
#' @param F_vals Matrix or data.frame of F statistics.
#' @param q_vals Matrix or data.frame of q values (optional).
#' @param n_rho_w Number of LD-window draws.
#' @param rho_w_lim List(min,max) for rho_w.
#' @param n_draws Number of OR parameter draws per rho_w.
#' @param rho_d_lim List(min,max) for rho_d.
#' @param rho_ld_lim List(min,max) for rho_ld.
#' @param alpha_lim List(min,max) for alpha exponent.
#' @param lmin_lim List(min,max) for minimum OR size.
#' @param seed Optional random seed.
#'
#' @return Object of class "ld_rho_draws".
#' @export
ld_rho_draws <- function(gds,
                         ld_decay,
                         n_inds,
                         F_vals     = NULL,
                         q_vals     = NULL,
                         C_scores   = NULL,
                         stat_type  = c("q","C"),
                         mode       = c("joint","per_method"),
                         n_rho_w    = 25,
                         ld_w_int   = NULL,
                         n_draws    = 100,
                         rho_w_lim  = list(min=0.8,max=0.99),
                         rho_d_lim  = list(min=0.5,max=0.999),
                         rho_ld_lim = list(min=0.9,max=0.999),
                         alpha_lim  = list(min=1.31,max=4),
                         lmin_lim   = list(min=1,max=10),
                         C_lim      = list(min=0,max=0.5),
                         cores=1

){

  #mode = "per_method"
  n_inds <- .get_n_inds(gds)

  ids <- .read_gds_ids(gds)

  if(!is.null(ld_w_int)){
    scan_int <- ld_scan(
      SNP_ids   = ids$snp_id,
      F_vals    = F_vals,
      ld_w      = ld_w_int,
      n_inds    = n_inds,
      full      = TRUE
    )

    q_primes_int <- do.call(cbind,lapply(scan_int$result,function(x) x$q_prime))
    colnames(q_primes_int) <- paste0(colnames(q_primes_int),"_prime_int")

  }else{
    q_primes_int = NULL
  }


  if(stat_type[1]=="q"){
    draws <- rbindlist(lapply(seq_len(n_rho_w),function(dr){
      #message("draw ",dr)

      if(!is.null(rho_w_lim)){
        rho_w <- runif(1, rho_w_lim$min,rho_w_lim$max)

        ld_w <- compute_ld_w(ld_decay,
                             rho = rho_w,
                             cores = cores)

        scan <- ld_scan(
          SNP_ids   = ids$snp_id,
          F_vals    = F_vals,
          ld_w      = ld_w,
          n_inds    = n_inds,
          full      = TRUE
        )


        q_primes <- do.call(cbind,lapply(scan$result,function(x) x$q_prime))
        colnames(q_primes) <- paste0(colnames(q_primes),"_prime")
        qvals <- cbind(q_vals, q_primes,q_primes_int)

      }else{
        rho_w = NA
        qvals <- cbind(q_vals, q_primes_int)
      }

      q_min <- apply(qvals,1,min)

      ## pre-estimate all pairwise LD values for outliers
      idx <- which(q_min<1/10^alpha_lim$min)

      if(length(idx)>0)
        el  <- get_el(gds, idx, slide_win_ld = -1,by_chr = TRUE)
      else
        el  <- data.table(matrix(0,0,0))

      #n_or_draws=100
      draws <- or_draws(
        el         = el,
        vals       = qvals,
        SNP_ids    = ids$snp_id,
        SNP_chr    = ids$snp_chr,
        ld_decay  = ld_decay,
        n_draws    = n_draws,
        rho_d_lim  = rho_d_lim,
        rho_ld_lim = rho_ld_lim,
        alpha_lim  = alpha_lim,
        lmin_lim   = lmin_lim,
        C_lim      = NULL,
        mode       = mode[1],
        stat_type  = stat_type[1],
        cores      = cores
      )

      draws[,rho_w:=rho_w]
      return(draws)
    }))
  }
  #stat_type <- "C"
  if(stat_type[1]=="C"){


      idx <- which(C_scores>0)

      if(length(idx)>0)
        el  <- get_el(gds, idx, slide_win_ld = -1,by_chr = TRUE)
      else
        el  <- data.table(matrix(0,0,0))

      draws <- or_draws(
        el         = el,
        vals       = C_scores,
        SNP_ids    = ids$snp_id,
        SNP_chr    = ids$snp_chr,
        ld_decay   = ld_decay,
        n_draws    = n_draws,
        rho_d_lim  = rho_d_lim,
        rho_ld_lim = rho_ld_lim,
        alpha_lim  = NULL,
        lmin_lim   = NULL,
        C_lim      = C_lim,
        mode       = "per_method",
        stat_type  = "C",
        cores      = cores
      )

      draws[,rho_w:=NA]
    }

  list(draws=draws[,.(draw_id, method, rho_w, rho_d, rho_ld, alpha, l_min, OR, OR_size)],class = "ld_rho_draws")
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
