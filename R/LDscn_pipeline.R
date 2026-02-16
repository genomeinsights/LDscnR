#' Run complete LD-aware genome scan pipeline
#'
#' @param geno Optional genotype matrix.
#' @param map  SNP metadata (must contain Chr, Pos, marker).
#' @param gds  Optional open GDS object.
#' @param F_vals Matrix/data.table of F-statistics.
#' @param n_inds Number of individuals.
#' @param rho_w LD-decay quantile for single LD-scaled scan.
#' @param n_rho Number of LD-decay draws for robustness.
#' @param rho_w_lim Range of LD-decay quantiles for robustness.
#' @param detect_OR Logical; whether to detect ORs.
#' @param sign_th Significance threshold for OR detection.
#' @param sign_if "less" (q-values) or "greater" (C-score).
#' @param ... Additional arguments passed to ld_rho_draws().
#'
#' @export
LDscn_pipeline <- function(geno = NULL,
                        map,
                        gds = NULL,
                        F_vals,
                        q_vals,
                        n_rho = 40,
                        n_or_draws = 25,
                        rho_w_lim = list(min=0.8,max=0.99),
                        rho_d_lim=list(min=0.5,max=0.999),
                        rho_ld_lim=list(min=0.9,max=0.999),
                        alpha_lim=list(min=1.31,max=4),
                        lmin_lim=list(min=1,max=10)
                        ) {

  call <- match.call()

  # ------------------------------------------------------------
  # 1. GDS handling
  # ------------------------------------------------------------

  if (is.null(gds)) {
    if (is.null(geno))
      stop("Provide either 'gds' or 'geno'.")

    gds_path <- tempfile(fileext = ".gds")
    gds <- create_gds_from_geno(geno, map, gds_path)
  }


  SNP_ids <- .read_gds_ids(gds)$snp_id
  n_inds <- .get_n_inds(gds)
  # ------------------------------------------------------------
  # 2. LD decay + structure
  # ------------------------------------------------------------

  if (verbose)
    message("Estimating background LD and LD-decay")

  decay  <- ld_decay(gds)

  if (verbose)
    print(decay)

  if (verbose)
    message("Estimating LD-structure")
  ld_str <- compute_ld_structure(gds)

  if (verbose)
    print(ld_str)
  # ------------------------------------------------------------
  # 4. Multi-rho robustness
  # ------------------------------------------------------------

  if (verbose)
    message("Getting draws")
  draws <- ld_rho_draws(
    ld_struct  = ld_str,
    decay_obj  = decay,
    F_vals     = F_vals,
    q_vals     = q_vals,
    n_inds     = n_inds,
    SNP_ids    = SNP_ids,
    n_rho      = 40,
    n_or_draws = 25,
    rho_w_lim  = rho_w_lim,
    rho_d_lim  = rho_d_lim,
    rho_ld_lim = rho_ld_lim,
    alpha_lim  = alpha_lim,
    lmin_lim   = lmin_lim
  )

  consistency <- consistency_score(draws, combine = TRUE)

  out <- list(
    decay       = decay,
    ld_struct   = ld_str,
    draws       = draws,
    consistency = consistency,
    call        = call
  )

  class(out) <- "ld_analysis"
  out
}
