#' Draw OR-detection parameters across LD-support thresholds
#'
#' Generates repeated draws of outlier-region (OR) detection parameters,
#' optionally across multiple LD-support thresholds (\eqn{\rho_w}).
#'
#' For each value of \eqn{\rho_w}, local LD support (\code{ld_w}) is
#' recomputed from the fitted LD-decay object, genome scan statistics are
#' optionally transformed, and \code{or_draws()} is used to sample OR
#' detection parameters and identify candidate outlier regions.
#'
#' Two modes of operation are supported:
#' \enumerate{
#'   \item \code{stat_type = "q"}: repeated draws are performed for each
#'   value of \code{rho}, after recomputing \code{ld_w} and deriving
#'   \code{q_prime} statistics from \code{F_vals}.
#'   \item \code{stat_type = "C"}: repeated draws are performed directly on
#'   precomputed consistency scores, without looping over \code{rho}.
#' }
#'
#' @param gds GDS file handle containing genotype data.
#' @param ld_decay Object of class \code{"ld_decay"} returned by
#'   \code{compute_LD_decay()}.
#' @param F_vals Matrix or data.frame of test statistics on the F scale,
#'   typically one column per method. Required when \code{stat_type = "q"}.
#' @param q_vals Optional matrix or data.frame of q-values corresponding to
#'   \code{F_vals}. If provided, these are combined with derived
#'   \code{q_prime} values before OR detection.
#' @param C_scores Matrix or data.frame of consistency scores. Required when
#'   \code{stat_type = "C"}.
#' @param stat_type Character string specifying which input statistics are
#'   used. One of \code{"q"} or \code{"C"}.
#' @param mode Character string controlling how OR draws are performed.
#'   Passed to \code{or_draws()}. Typically one of \code{"joint"} or
#'   \code{"per_method"}.
#' @param n_draws Number of OR parameter draws per \code{rho} value.
#' @param rho Numeric vector of LD-support thresholds (\eqn{\rho_w}) used to
#'   recompute \code{ld_w}. Only used when \code{stat_type = "q"}.
#' @param ld_ws Pre-calculated matrix of \code{ld_w}-values from \code{precalculate_ld_w}.
#'   If Null, they are calculated based on the \code{rho}.
#' @param rho_d_lim Named list with elements \code{min} and \code{max}
#'   defining the sampling range for the distance threshold parameter
#'   \code{rho_d}.
#' @param rho_ld_lim Named list with elements \code{min} and \code{max}
#'   defining the sampling range for the LD-threshold parameter
#'   \code{rho_ld}.
#' @param alpha_lim Named list with elements \code{min} and \code{max}
#'   defining the sampling range for the significance exponent
#'   \code{alpha}. Only used when \code{stat_type = "q"}.
#' @param lmin_lim Named list with elements \code{min} and \code{max}
#'   defining the sampling range for the minimum outlier-region size
#'   \code{l_min}. Only used when \code{stat_type = "q"}.
#' @param C_lim Named list with elements \code{min} and \code{max}
#'   defining the sampling range for the minimum consistency threshold.
#'   Only used when \code{stat_type = "C"}.
#' @param cores Number of CPU cores used for parallel steps.
#'
#' @return A list containing a data table of OR-detection draws. Each row
#' represents one parameter draw and includes:
#' \describe{
#'   \item{method}{Method or method combination used for OR detection.}
#'   \item{rho_w}{LD-support threshold used to compute \code{ld_w}.}
#'   \item{rho_d}{Distance-based clustering threshold.}
#'   \item{rho_ld}{LD-based clustering threshold.}
#'   \item{alpha}{Significance exponent used for q-value thresholding.}
#'   \item{l_min}{Minimum outlier-region size.}
#'   \item{OR}{Detected outlier regions.}
#'   \item{OR_size}{Number or size of detected outlier regions.}
#' }
#'
#' @details
#' When \code{stat_type = "q"}, the function first computes \code{ld_w} for
#' each requested \code{rho} value using \code{compute_ld_w()}, then runs
#' \code{ld_scan()} to derive \code{q_prime} statistics from \code{F_vals}.
#' These are combined with the supplied \code{q_vals}, and SNPs passing a
#' liberal prefilter are used to construct an LD edge list for repeated OR
#' detection.
#'
#' When \code{stat_type = "C"}, repeated OR detection is applied directly to
#' positive consistency scores. In this mode, \code{rho} is ignored and
#' \code{rho_w} is returned as \code{NA}.
#'
#' This function is primarily intended for sensitivity analyses over LD
#' window definitions and OR-detection parameters.
#'
#' @export
ld_rho_draws <- function(gds,
                         ld_decay,
                         F_vals     = NULL,
                         q_vals     = NULL,
                         C_scores   = NULL,
                         stat_type  = c("q","C"),
                         mode       = c("joint","per_method"),
                         n_draws    = 100,
                         rho        = seq(0.75,0.95,by=0.025),
                         ld_ws      = NULL,
                         rho_d_lim  = list(min=0.5,max=0.999),
                         rho_ld_lim = list(min=0.9,max=0.999),
                         alpha_lim  = list(min=1.31,max=4),
                         lmin_lim   = list(min=1,max=10),
                         C_lim      = list(min=0,max=0.5),
                         cores=1

){

  n_inds <- .get_n_inds(gds)

  ids <- .read_gds_ids(gds)


  #rho = 0.9

  if(is.null(ld_ws)){
    ld_ws <- precalculate_ld_w(rho = rho, ld_decay=ld_decay)
  }else{
    rho = colnames(ld_ws)
  }



  message("Getting draws")
  if(stat_type[1]=="q"){
    if(length(rho)>1){
      pb <- txtProgressBar(min = 0, max = length(rho)-1, style = 3)
      setTxtProgressBar(pb, 0)
    }
    draws <- rbindlist(lapply(rho,function(rh){



      scan <- ld_scan(
        SNP_ids   = ids$snp_id,
        F_vals    = F_vals,
        ld_w      = ld_ws[,which(colnames(ld_ws) == rh)],
        n_inds    = n_inds,
        full      = TRUE
      )

      q_primes <- do.call(cbind,lapply(scan$result,function(x) x$q_prime))
      colnames(q_primes) <- paste0(colnames(q_primes),"_prime")
      qvals <- cbind(q_vals, q_primes)


      ## pre-estimate all pairwise LD values for outliers
      q_min <- apply(qvals,1,min)
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
        ld_decay   = ld_decay,
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

      draws[,rho_w:=rh]

      if(length(rho)>1)
        setTxtProgressBar(pb, which(rho==rh))

      return(draws)
    }))
    if(length(rho)>1)
      close(pb)
  }



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

  out <- list(
    draws = draws[,.(method, rho_w, rho_d, rho_ld, alpha, l_min, OR, OR_size)],
    n_rho = if (stat_type[1] == "q") length(rho) else 1L,
    n_or_draws = n_draws,
    stat_type = stat_type[1],
    mode = mode[1]
  )
  class(out) <- "ld_rho_draws"
  out
}


#' Pre-calculates \code{ld_w}'s for \code{"ld_rho_draws"}
#'
#' Pre-calculates \code{ld_w}'s for use with \code{"ld_rho_draws"}. Saves time, especially if edge lists are saved to file by \code{"compute_LD_decay"}
#'
#' @param rho A vector o rho values. Can be a single value.
#' @param ld_decay Object of class \code{"ld_decay"} returned by
#'   \code{compute_LD_decay()}.
#'
#' @return matrix of \code{ld_w}'s for SNPs with each column corresponding to a different rho value
#'
#' @export
precalculate_ld_w <- function(rho,ld_decay){
  message("Precalculating ld_w")
  if(length(rho)>1){
    pb <- txtProgressBar(min = 0, max = length(rho)-1, style = 3)
    setTxtProgressBar(pb, 0)
  }
  #chr_obj <- ld_decay$by_chr[[1]]
  ld_ws <- do.call(rbind,lapply(ld_decay$by_chr, function(chr_obj) {

    a <- ld_decay$decay_sum[Chr==chr_obj$decay_sum$Chr,a_pred]
    b <- ld_decay$decay_sum[Chr==chr_obj$decay_sum$Chr,b]

    d_window <- d_from_rho(a, rho)

    if(is.null(chr_obj$el)) stop("No edge list present")

    if(is.character(chr_obj$el)) chr_obj$el <- fread(chr_obj$el,showProgress = FALSE)
    #make symmetric
    chr_obj$el <- data.table::rbindlist(list(
      chr_obj$el[, .(SNP = SNP1, pos = pos1, pos_other = pos2, r2, d)],
      chr_obj$el[, .(SNP = SNP2, pos = pos2, pos_other = pos1, r2, d)]
    ))

    #win = 300000
    ld_w <-do.call(cbind,lapply(d_window,function(win){
      ld_w <- chr_obj$el[d<win,.(r2_median=median(r2)),by=SNP]
      ld_w[match(chr_obj$snp_ids,ld_w$SNP)]$r2_median
    }))


    if(length(rho)>1)
      setTxtProgressBar(pb, which(names(ld_decay$by_chr)==chr_obj$decay_sum$Chr))

    return(ld_w)

  }))

  if(length(rho)>1)
    close(pb)

  colnames(ld_ws) <-  rho
  return(ld_ws)
}



#' Print method for \code{ld_rho_draws} objects
#'
#' Displays a concise summary of repeated OR-detection draws across
#' LD-support thresholds.
#'
#' @param x Object of class \code{"ld_rho_draws"}.
#'
#' @return The input object, invisibly.
#'
#' @export
print.ld_rho_draws <- function(x) {
  methods <- unique(x$draws$method)
  cat("\nLD rho-window draws\n")
  cat("--------------------\n")
  cat("Number of rho_w draws:", x$n_rho, "\n")
  cat("OR draws per rho_w:", x$n_or_draws, "\n")

  cat("Included methods:", paste(methods,collapse="/"), "\n")
  invisible(x)
}
