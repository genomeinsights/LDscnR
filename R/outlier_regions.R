#' Detect outlier regions from SNP-level statistics
#'
#' Identifies outlier regions (ORs) by clustering outlier SNPs that are both
#' physically close and in sufficiently strong linkage disequilibrium (LD).
#'
#' SNPs are first classified as outliers according to a user-supplied threshold.
#' Outlier SNPs are then grouped into connected components within each chromosome
#' using pairwise LD edges that satisfy both a distance threshold and an
#' LD threshold. These thresholds are derived from chromosome-specific LD-decay
#' parameters.
#'
#' Two detection modes are supported:
#' \enumerate{
#'   \item \code{per_method}: detect outlier regions separately for each column
#'   of \code{vals},
#'   \item \code{joint}: take the union of outlier SNPs across methods before
#'   clustering, and return a single joint set of outlier regions.
#' }
#'
#' @param el LD edge list containing pairwise SNP LD values and genomic distances.
#' @param vals Matrix or data.frame of SNP-level statistics, with SNPs in rows and
#'   methods in columns.
#' @param ld_decay Object of class \code{"ld_decay"} containing chromosome-specific
#'   LD-decay summaries.
#' @param SNP_ids Vector of SNP identifiers in the same order as \code{vals}.
#' @param SNP_chr Vector of chromosome identifiers in the same order as \code{vals}.
#' @param sign_th Significance threshold used to classify SNPs as outliers.
#' @param sign_if Character string indicating the outlier rule:
#'   \code{"less"} for lower-tail statistics (e.g. q-values),
#'   \code{"greater"} for upper-tail statistics (e.g. consistency scores).
#' @param rho_d Relative distance threshold used to derive the maximum allowed
#'   physical distance between SNPs in the same outlier region.
#' @param rho_ld Relative LD threshold used to derive the minimum required
#'   pairwise LD for SNPs in the same outlier region.
#' @param l_min Minimum number of SNPs required for an outlier region.
#' @param mode One of \code{"per_method"} or \code{"joint"}.
#' @param ret_table Logical; if \code{TRUE}, return a long-format table instead
#'   of an \code{ld_or} object.
#'
#' @return
#' If \code{ret_table = FALSE}, an object of class \code{"ld_or"} containing:
#' \describe{
#'   \item{ORs}{Named list of detected outlier regions. Each region is a character
#'   vector of SNP identifiers.}
#'   \item{params}{List of parameters used for detection.}
#' }
#'
#' If \code{ret_table = TRUE}, a long \code{data.table} with one row per SNP-region
#' assignment.
#'
#' @details
#' For each chromosome, the function retrieves chromosome-specific LD-decay
#' parameters from \code{ld_decay$decay_sum}. The relative thresholds
#' \code{rho_d} and \code{rho_ld} are converted to an absolute distance threshold
#' and an absolute LD threshold using \code{d_from_rho()} and \code{ld_from_rho()}.
#'
#' Outlier regions are defined as connected components in an undirected graph,
#' where nodes are outlier SNPs and edges connect SNP pairs satisfying both
#' thresholds.
#'
#' The expected direction of significance depends on the statistic type:
#' q-like measures typically use \code{sign_if = "less"}, whereas
#' C-like measures typically use \code{sign_if = "greater"}.
#'
#' @export
detect_or <- function(el,
                      vals,
                      ld_decay,
                      SNP_ids,
                      SNP_chr,
                      sign_th = 0.05,
                      sign_if = c("less", "greater"),
                      rho_d = 0.99,
                      rho_ld = 0.99,
                      l_min = 1,
                      mode = c("per_method","joint"),
                      ret_table = FALSE) {


  if (is.null(colnames(vals)))
    stop("`vals` must have column names.")

  methods <- colnames(vals)

  mode <- match.arg(mode)

  #stat_type <- infer_stat_type(methods)

  outlier_fun <- switch(
    sign_if,
    less    = function(x) x < sign_th,
    greater = function(x) x > sign_th
  )

  chrs <- unique(SNP_chr)

  # Outlier rule
  outlier_fun <- switch(
    sign_if[1],
    less    = function(x) x < sign_th,
    greater = function(x) x > sign_th
  )


  # ------------------------------------------------------------
  # Helper: build ORs from SNP set
  # ------------------------------------------------------------
  build_or_from_snps <- function(outliers) {

    if (length(outliers) == 0)
      return(list())

    idx <- which(SNP_ids %in% outliers)

    or_list <- list()

    for (ch in chrs) {
    #cat(ch," -- ")

      out_chr <- SNP_ids[idx][SNP_chr[idx] == ch]

      if (length(out_chr) < l_min)
        next

      # thresholds
      a_chr  <- ld_decay$decay_sum[Chr == ch, a_pred]
      b_chr  <- ld_decay$decay_sum[Chr == ch, b]
      c_chr  <- ld_decay$decay_sum[Chr == ch, c_pred]


      d_th  <- d_from_rho(a_chr, rho = rho_d)
      ld_th <- ld_from_rho(b_chr, c_chr, rho = rho_ld)

      ed <- el[
        d < d_th &
          r2 > ld_th &
          SNP1 %in% out_chr &
          SNP2 %in% out_chr,
        .(SNP1, SNP2)
      ]

      if (nrow(ed) == 0)
        next

      g     <- igraph::graph_from_data_frame(ed, directed = FALSE)
      comps <- igraph::components(g)

      ors_chr <- split(names(comps$membership), comps$membership)
      ors_chr <- ors_chr[vapply(ors_chr, length, integer(1)) >= l_min]

      if (length(ors_chr) > 0)
        or_list <- c(or_list, ors_chr)
    }

    or_list
  }

  # ------------------------------------------------------------
  # PER-METHOD MODE (existing behavior)
  # ------------------------------------------------------------
  if (mode == "per_method") {
    #qs <- unlist(q_vals[,1])

    out <- apply(as.matrix(vals), 2, function(qs) {

      outliers <- SNP_ids[!is.na(qs) & outlier_fun(qs)]
      build_or_from_snps(outliers)
    })

  }

  # ------------------------------------------------------------
  # JOINT MODE (UNION BEFORE CLUSTERING)
  # ------------------------------------------------------------
  if (mode == "joint") {

    # union of outliers across methods
    qs <- apply(vals, 1, min)

    outliers <- SNP_ids[!is.na(qs) & outlier_fun(qs)]

    #union_outliers <- unique(unlist(outlier_matrix))

    joint_or <- build_or_from_snps(outliers)

    # return as single named element
    out <- list(Joint = joint_or)
  }

  out <- structure(
    list(
      ORs     = out,
      params= list (sign_th = sign_th,
                    sign_if = sign_if,
                    rho_d   = rho_d,
                    rho_ld  = rho_ld,
                    l_min   = l_min,
                    mode    = mode)

    ),
    class = "ld_or"
  )

  if (ret_table)
    return(or_to_table(out))

  out
}



#' Sample outlier regions across detection parameters
#'
#' Repeatedly draws outlier-region detection parameters from user-specified
#' ranges and applies \code{detect_or()} to generate a flat table of outlier
#' region (OR) results across parameter combinations.
#'
#' This function is intended for sensitivity analyses and robustness assessment
#' of OR detection. Each draw samples clustering thresholds and significance
#' parameters, detects outlier regions, and stores the resulting ORs together
#' with the sampled parameter values.
#'
#' @param el LD edge list containing pairwise SNP LD values and genomic distances.
#' @param vals Matrix or data.frame of SNP-level statistics, with SNPs in rows and
#'   methods in columns.
#' @param SNP_ids Vector of SNP identifiers in the same order as \code{vals}.
#' @param SNP_chr Vector of chromosome identifiers in the same order as \code{vals}.
#' @param ld_decay Object of class \code{"ld_decay"} containing chromosome-specific
#'   LD-decay summaries.
#' @param n_draws Number of parameter draws.
#' @param stat_type Type of statistic used for OR detection; one of \code{"q"}
#'   or \code{"C"}.
#' @param mode Detection mode passed to \code{detect_or()}, either
#'   \code{"per_method"} or \code{"joint"}. Multiple modes may be supplied.
#' @param rho_d_lim Named list with elements \code{min} and \code{max} defining
#'   the sampling range for the relative distance threshold \code{rho_d}.
#' @param rho_ld_lim Named list with elements \code{min} and \code{max} defining
#'   the sampling range for the relative LD threshold \code{rho_ld}.
#' @param alpha_lim Named list with elements \code{min} and \code{max} defining
#'   the sampling range for the significance exponent when
#'   \code{stat_type = "q"}.
#' @param C_lim Named list with elements \code{min} and \code{max} defining
#'   the sampling range for the consistency threshold when
#'   \code{stat_type = "C"}.
#' @param lmin_lim Named list with elements \code{min} and \code{max} defining
#'   the sampling range for the minimum region size \code{l_min}.
#' @param cores Number of CPU cores.
#'
#' @return A \code{data.table} with one row per draw and method, containing:
#' \describe{
#'   \item{draw_id}{Identifier of the parameter draw.}
#'   \item{method}{Method name or joint-analysis label.}
#'   \item{rho_d}{Sampled relative distance threshold.}
#'   \item{rho_ld}{Sampled relative LD threshold.}
#'   \item{alpha}{Sampled significance threshold parameter. For q-based input,
#'   this is the actual q-value threshold derived from the sampled exponent.}
#'   \item{l_min}{Sampled minimum outlier-region size.}
#'   \item{OR}{List-column containing detected outlier regions.}
#'   \item{OR_size}{Number of detected outlier regions.}
#' }
#'
#' @details
#' For q-based statistics, the significance threshold is sampled on a
#' log-scale by drawing an exponent and converting it to
#' \eqn{1 / 10^{\alpha}}. For C-based statistics, thresholds are sampled
#' directly from \code{C_lim}.
#'
#' If no LD edges are available, the function returns empty OR results for all
#' draws.
#'
#' This is an internal workflow helper but may also be useful for advanced users
#' exploring OR parameter sensitivity.
or_draws <- function(el,
                     vals,
                     SNP_ids,
                     SNP_chr,
                     ld_decay,
                     n_draws    = 25,
                     stat_type  = c("q","C"),
                     mode       = c("per_method","joint"),
                     rho_d_lim  = list(min=0.5,max=0.999),
                     rho_ld_lim = list(min=0.9,max=0.999),
                     alpha_lim  = list(min=1.31,max=4),
                     C_lim      = list(min=0,max=0.5),
                     lmin_lim   = list(min=1,max=10),
                     cores      = 1
                     ) {

  stat_type <- match.arg(stat_type)

  if (is.null(colnames(vals)))
    stop("Values must have column names.")
  #mod <- "joint"

  out <- rbindlist(lapply(mode,function(mod){

    methods <- colnames(vals)
    if(cores==1) cat("replicate:")
    out <- rbindlist(parallel_apply(seq_len(n_draws),function(i) {
      if(cores==1) cat(i, "- ")

      rho_d  <- runif(1, rho_d_lim$min, rho_d_lim$max)
      rho_ld <- runif(1, rho_ld_lim$min, rho_ld_lim$max)

      l_min    <- if (stat_type[1] == "C") 1 else sample(seq(lmin_lim$min, lmin_lim$max), 1)
      alpha    <- if (stat_type[1] == "C") runif(1, C_lim$min, C_lim$max) else 1 / 10^(runif(1, alpha_lim$min, alpha_lim$max))
      sign_if  <- if (stat_type[1] == "C") "greater" else "less"


      if(length(el)==0)
        return(data.table(
          draw_id = i,
          method  = methods,
          rho_d   = rho_d,
          rho_ld  = rho_ld,
          alpha   = alpha,
          l_min   = l_min,
          OR      = replicate(length(methods), list(list()), simplify = FALSE),
          OR_size = 0L
        ))
      #mod = "joint"
      or_obj <- detect_or(
        el        = el,
        vals      = vals,
        SNP_ids   = SNP_ids,
        SNP_chr   = SNP_chr,
        ld_decay  = ld_decay,
        sign_th   = alpha,
        sign_if   = sign_if,
        rho_d     = rho_d,
        rho_ld    = rho_ld,
        l_min     = l_min,
        ret_table = FALSE,
        mode      = mod
      )

      OR_methods <- names(or_obj$ORs)
      if (all(lengths(or_obj$ORs) == 0)) {
        data.table(
          draw_id = i,
          method  = OR_methods,
          rho_d   = rho_d,
          rho_ld  = rho_ld,
          alpha   = alpha,
          l_min   = l_min,
          OR      = replicate(length(OR_methods), list(list()), simplify = FALSE),
          OR_size = 0L
        )
        #meth  <- OR_methods[1]
      }else{
        rbindlist(lapply(OR_methods, function(meth) {

          or_list <- or_obj$ORs[[meth]]
          if (length(or_list) == 0) return(
            data.table(
              draw_id = i,
              method  = meth,
              rho_d   = rho_d,
              rho_ld  = rho_ld,
              alpha   = alpha,
              l_min   = l_min,
              OR      = list(list()),
              OR_size = 0L
            ))

          data.table(
            draw_id = i,
            method  = meth,
            rho_d   = rho_d,
            rho_ld  = rho_ld,
            alpha   = alpha,
            l_min   = l_min,
            OR      = list(or_list),
            OR_size = length(or_list)
          )

        }), fill = TRUE)
      }

    },cores=cores))
    if(cores==1) cat("\n")
    return(out)
  }))


}


#' Convert outlier regions to long format
#'
#' Converts an object of class \code{"ld_or"} into a long-format table with one
#' row per SNP assigned to an outlier region.
#'
#' @param or_obj Object of class \code{"ld_or"}.
#'
#' @return A \code{data.table} with columns:
#' \describe{
#'   \item{SNP}{SNP identifier.}
#'   \item{method}{Method or analysis mode under which the region was detected.}
#'   \item{OR_id}{Outlier-region identifier.}
#' }
#'
#' @details
#' Region identifiers are constructed as \code{<Chr>_OR<i>}, where the
#' chromosome is inferred from the first SNP in each region.
#'
#' If no outlier regions are present, an empty \code{data.table} is returned.
#'
#' @export
or_to_table <- function(or_obj) {

  if (!inherits(or_obj, "ld_or"))
    stop("Input must be class 'ld_or'.")

  res <- list()

  for (method in names(or_obj$ORs)) {

    or_list <- or_obj$ORs[[method]]
    if (length(or_list) == 0) next

    for (i in seq_along(or_list)) {

      snps <- or_list[[i]]

      # Extract chromosome from first SNP
      chr <- sub(":.*", "", snps[1])

      or_name <- paste0(chr, "_OR", i)

      res[[length(res) + 1]] <- data.table::data.table(
        SNP    = snps,
        method = method,
        OR_id  = or_name
      )
    }
  }

  if (length(res) == 0)
    return(data.table::data.table(SNP=character(), method=character(), OR_id=character()))

  data.table::rbindlist(res)
}


#' Infer statistic type from method names
#'
#' Internal helper that determines whether a set of method names corresponds
#' to q-based statistics or C-based statistics.
#'
#' @param method_names Character vector of method names.
#'
#' @return Character scalar: \code{"q"} or \code{"C"}.
#'
#' @details
#' Method names beginning with \code{"C"} or containing \code{"_C"} are treated
#' as C-based. All others are treated as q-based. Mixed inputs are rejected.
#'
#' @keywords internal
infer_stat_type <- function(method_names) {

  has_C  <- grepl("^C", method_names) | grepl("_C", method_names)
  has_q  <- !has_C

  if (any(has_C) && any(has_q)) {
    stop("Cannot mix C-based and q-value-based methods in the same call.")
  }

  if (all(has_C)) return("C")
  if (all(has_q)) return("q")

  stop("Unable to infer statistic type.")
}
