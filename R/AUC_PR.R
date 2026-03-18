#' Compute precision-recall summaries from simulated outlier regions
#'
#' Evaluates detected outlier regions (ORs) against known causal loci in
#' simulated data and computes per-draw precision, recall, and related
#' summary statistics.
#'
#' For each set of detected outlier regions, the function assigns each region
#' to a focal QTN using the SNP within that region with the highest
#' \code{rho_ld} value, and compares the resulting set of focal QTNs to the
#' known true QTNs.
#'
#' @param draw_ORs A list of detected outlier-region sets, typically one
#'   element per draw. Each element should itself be a list of regions, where
#'   each region is a character vector of SNP identifiers.
#' @param map_filt SNP annotation table for the simulated data. Must contain
#'   at least the columns \code{marker}, \code{p_Va}, \code{focal_QTN},
#'   \code{rho_ld}, and \code{Chr_type}.
#' @param p_Va_th Threshold used to define true focal QTNs from \code{p_Va}.
#' @param cores Number of CPU cores.
#'
#' @return A \code{data.table} with one row per OR draw and columns:
#' \describe{
#'   \item{TP}{Number of true-positive outlier regions.}
#'   \item{FP}{Number of false-positive outlier regions.}
#'   \item{FN}{Number of true focal QTNs not recovered by any detected region.}
#'   \item{precision}{Precision, computed as \code{TP / (TP + FP)}.}
#'   \item{recall}{Recall, computed as \code{TP / (TP + FN)}.}
#'   \item{PR}{Product of precision and recall.}
#'   \item{FP_ntrl}{Number of detected regions overlapping neutral chromosomes.}
#'   \item{OR_focals}{List-column containing the focal QTNs assigned to the
#'   detected outlier regions.}
#' }
#'
#' @details
#' True positives are defined at the level of focal QTNs rather than individual
#' SNPs. A detected outlier region is counted as a true positive if its assigned
#' focal QTN belongs to the set of true QTNs defined by \code{p_Va > p_Va_th}.
#'
#' The \code{PR} summary returned here is the product of precision and recall,
#' not the area under a precision-recall curve.
#' @export
get_PR <- function(draw_ORs,map,map_filt,p_Va_th=0.05,cores=1) {


  true <- unique(map_filt[p_Va > p_Va_th, focal_QTN])

  rbindlist(parallel_apply(draw_ORs,function(ORs){

    if (length(ORs) > 0) {

      #map_filt[marker %in% unlist(ORs)]
      # Identify focal per OR
      OR_focals <- vapply(ORs, function(or) {

        sub <- map_filt[marker %in% or]

        if (nrow(sub) == 0L) return(NA_character_)

        sub[which.max(rho_ld), focal_QTN]

      }, character(1))

      OR_focals <- unique(na.omit(OR_focals))

      # Confusion components
      TP <- sum(OR_focals %in% true)
      FP <- length(ORs) - TP
      FN <- length(setdiff(true, OR_focals))

      precision <- if ((TP + FP) > 0) TP / (TP + FP) else 0
      recall    <- if ((TP + FN) > 0) TP / (TP + FN) else 0
      PR        <- precision * recall

      # Neutral FP check (vectorized)
      neutral_markers <- unique(map[Chr_type == "ntrl", marker])

      FP_ntrl <- sum(vapply(ORs, function(or) {
        any(or %in% neutral_markers)
      }, logical(1)))

      return(data.table(
        TP = TP,
        FP = FP,
        FN = FN,
        precision = precision,
        recall = recall,
        PR = PR,
        FP_ntrl = FP_ntrl,
        OR_focals = list(list(OR_focals))
      ))

    } else {

      return(data.table(
        TP = 0,
        FP = 0,
        FN = length(true),
        precision = 0,
        recall = 0,
        PR = 0,
        FP_ntrl = 0,
        OR_focals = list(list(character()))
      ))
    }
  },cores=cores))


}

#' Summarize OR performance using cumulative precision-recall trajectories
#'
#' Estimates the area under a cumulative precision-recall trajectory derived
#' from repeated bootstrap resampling of OR evaluation results.
#'
#' For each method, the function repeatedly resamples \code{PR} values,
#' computes the running maximum across draws, averages these trajectories
#' across bootstrap replicates, and calculates the area under the resulting
#' curve.
#'
#' @param PR_data A \code{data.table} containing at least the columns
#'   \code{method} and \code{PR}, typically returned by \code{get_PR()} and
#'   combined across methods.
#' @param n1 Number of bootstrap replicates used to estimate the average
#'   cumulative trajectory.
#' @param n2 Number of sampled draws per bootstrap replicate.
#' @param cores Number of CPU cores.
#'
#' @return A list with two components:
#' \describe{
#'   \item{data}{A \code{data.table} containing the mean cumulative PR
#'   trajectory for each method.}
#'   \item{AUC}{A \code{data.table} summarizing raw and normalized AUC values
#'   for each method.}
#' }
#'
#' @details
#' The normalized AUC rescales the observed area relative to the theoretical
#' maximum cumulative trajectory of 1 across draws.
#'
#' This summary is intended for simulation-based benchmarking of OR detection
#' performance across methods.
#' @export
get_AUC_OR <- function(PR_data, n1=500, n2=500,cores=1){

  ## get running maximum

  cummax_PR <- rbindlist(parallel_apply(1:n1, function(x){
    rbindlist(lapply(PR_data[,unique(method)],function(meth){
      cbind(indx=1:n2,method=meth,PR_data[method==meth][sample(1:.N,n2,replace = TRUE),.(value=cummax(PR),rep=x)])
    }))
  },cores = cores))

  AUC_norm <- list()
  AUC <- list()

  ## get AUC for each method
  for(meth in cummax_PR[,unique(method)]){
    PR_star <- cummax_PR[ method==meth,mean(value),by=.(method,indx)][,V1]
    k <- seq_along(PR_star)
    AUC_pracma <- pracma::trapz(k, PR_star)
    PR_star_max <- pracma::trapz(k, c(0,rep(1,(max(k)-1))))
    AUC_norm[[meth]] <- AUC_pracma / PR_star_max
    AUC[[meth]] <- AUC_pracma
  }

  data_auc <- cummax_PR[,.(mean_PR=mean(value)),by=.(method,indx)]

  max_PR <- cummax_PR[,max(value),by=method]
  max_PR <- max_PR[match(names(AUC),max_PR$method),V1]
  AUC <- data.table(method=names(AUC),AUC=unlist(AUC),AUC_norm=unlist(AUC_norm),max_PR)

  return(list(data=data_auc,AUC=AUC))
}

