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
                        F_cols,
                        q_cols=NULL,
                        slide_win_ld=1000,
                        n_rho = 40,
                        n_or_draws = 25,
                        rho_w_lim = list(min=0.8,max=0.99),
                        rho_d_lim=list(min=0.9,max=0.999),
                        rho_ld_lim=list(min=0.5,max=0.999),
                        alpha_lim=list(min=1.31,max=4),
                        lmin_lim=list(min=1,max=10),
                        cores=1,
                        verbose=TRUE,
                        keep_draws = FALSE,
                        keep_ld_str = TRUE,
                        keep_ld_decay = TRUE) {

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

  if (verbose)
    message("Estimating LD-structure")

  ld_str <- compute_ld_structure(gds,slide_win_ld = slide_win_ld,cores=1)


  if (verbose)
    print(ld_str)

  decay  <- ld_decay(gds,ld_struct = ld_str,cores=cores)

  if (verbose)
    print(decay)

  # ------------------------------------------------------------
  # 4. Multi-rho robustness
  # ------------------------------------------------------------

  if (verbose)
    message("Getting draws")

  draws <- ld_rho_draws(
    ld_struct  = ld_str,
    decay_obj  = decay,
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
    cores      = cores
  )

  consistency <- consistency_score(draws)

  out <- list(
    decay       = if (keep_ld_decay) decay else NULL,
    ld_struct   = if (keep_ld_str) ld_str else NULL,
    draws       = if (keep_draws) draws else NULL,
    consistency = consistency,
    call        = call
  )

  class(out) <- "LDscn"
  out
}
