#' @keywords internal
"_PACKAGE"

#' @import data.table
#' @importFrom graphics abline axis barplot legend lines matplot par plot.new points text title
#' @importFrom stats coef median na.omit nls nls.control predict setNames
#' @importFrom utils setTxtProgressBar txtProgressBar
NULL

## data.table columns referenced by non-standard evaluation across the package,
## declared here so R CMD check does not flag them as undefined globals.
utils::globalVariables(c(
  ".", ".SD", ".data", "..pos_cols", "..rec_cols",
  "CL_col", "Chr1", "Chr2", "SNP", "SNP1", "SNP2", "Var1", "Var2",
  "a", "a_pred", "c_pred", "cand", "chr_size", "cl_rank", "cluster", "contrast",
  "core_snp", "d", "emp_p", "emp_q", "end", "group_id", "grp", "is_core", "ldw0",
  "median_ld", "members", "n", "n_sig", "n_snp_chr", "n_snps", "n_w_used",
  "pos1", "pos2", "pv0", "q_snp", "r2_eMLG", "r2_median", "regime",
  "rho_slide_pred", "rho_slide_raw", "sig", "sig_stat", "significant",
  "slide_snp", "start", "strata",
  ## added for ld_outlier_test()/ld_outlier_perm()/ld_region_rotation()
  "from", "to", "run", "n_markers", "occupancy", "n_units", "unit_id", "chr", "w"))
