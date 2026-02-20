#' Run complete LD-aware genome scan pipeline
#'
#' @param geno Genotype matrix in 012-format. Expected genotypes are rounded to nearest integer for gds.
#' @param map  SNP metadata; must contain Chr, Pos, marker and columns for the relevant test statistics.
#' @param F_cols Vector of column names for F-statistics .
#' @param q_cols Optional vector of q-values.
#' @param rho_w LD-decay quantile for single LD-scaled scan.
#' @param n_rho Number of LD-decay draws for robustness.
#' @param rho_w_lim Range of LD-decay quantiles for robustness.
#' @param rho_d_lim Range of LD-decay quantiles for robustness.
#' @param rho_ld_lim Range of LD-decay quantiles for robustness.'
#' @param keep_draws keep draws.
#' @param keep_ld_str keep ld_structure.
#' @param keep_ld_decay keep decay.
#' @param cores N cores.
#' @param detect_OR Logical; whether to detect ORs.
#' @param sign_th Significance threshold for OR detection.
#' @param sign_if "less" (q-values) or "greater" (C-score).
#' @param ... Additional arguments passed to ld_rho_draws().
#'
#' @export
LDscn_pipeline <- function(geno = NULL,
                        map,
                        ## for LD-structure estimation
                        q = 0.95,
                        n_sub = 5000,
                        window_size = 1e6,
                        step_size = 5e5,
                        dist_unit = 5000,
                        K_target = 1000,
                        n_win = 1000,
                        prob_robust = 0.5,
                        use="robust",
                        ## input F and q-values
                        F_cols,
                        q_cols=NULL,
                        ## Draw parameters
                        n_rho = 40,
                        n_or_draws = 25,
                        rho_w_lim = list(min=0.5,max=0.95),
                        rho_d_lim=list(min=0.9,max=0.999),
                        rho_ld_lim=list(min=0.5,max=0.999),
                        alpha_lim=list(min=1.31,max=4),
                        lmin_lim=list(min=1,max=10),
                        cores=1,
                        max_rho=0.99,
                        ## output
                        verbose=TRUE) {

  call <- match.call()

  if(cores>1)
    message("You are using multiple cores, make sure you have enough RAM")
  # ------------------------------------------------------------
  # 1. GDS handling
  # ------------------------------------------------------------

  SNP_ids <- map$marker
  n_inds <- nrow(geno)

  message("Creating gds file")

  gds_path <- tempfile(fileext = ".gds")
  gds <- create_gds_from_geno(geno, map, gds_path)
  on.exit({ SNPRelate::snpgdsClose(gds); unlink(gds_path) }, add = TRUE)

  # ------------------------------------------------------------
  # 2. LD decay + structure
  # ------------------------------------------------------------

  ld_struct <-  compute_ld_structure(gds,
                                     q           = q,
                                     n_sub       = n_sub,
                                     window_size = window_size,
                                     step_size   = step_size,
                                     cores       = cores,
                                     dist_unit   = dist_unit,
                                     K_target    = K_target,
                                     n_win       = n_win,
                                     prob_robust = prob_robust,
                                     use         = use,
                                     max_rho     = max_rho)

  #plot(ld_struct,type = "decay")
  # ------------------------------------------------------------
  # 4. Multi-rho robustness
  # ------------------------------------------------------------

  if (verbose)
    message("Getting draws")

  draws <- ld_rho_draws(
    ld_struct  = ld_struct,
    F_vals     = map[,..F_cols],
    q_vals     = map[,..q_cols],
    n_inds     = n_inds,
    SNP_ids    = SNP_ids,
    n_rho      = n_rho,
    n_or_draws = n_or_draws,
    rho_w_lim  = rho_w_lim,
    rho_d_lim  = rho_d_lim,
    rho_ld_lim = rho_ld_lim,
    alpha_lim  = alpha_lim,
    lmin_lim   = lmin_lim,
    cores      = cores,
    mode       = c("per_method","joint")
  )


  consistency <- consistency_score(draws)

  out <- list(
    ld_struct   = ld_struct,
    draws       = ld_struct,
    consistency = consistency,
    call        = call
  )

  class(out) <- "LDscn"
  out
}
