#' Simulated example dataset for LDscnR
#'
#' A simulated genome scan dataset including:
#' \itemize{
#'   \item \code{GTs}: genotype matrix (individuals × SNPs)
#'   \item \code{SNP_res}: SNP-level statistics and annotations
#' }
#'
#' @format A list with components:
#' \describe{
#'   \item{GTs}{Numeric matrix of genotypes}
#'   \item{SNP_res}{data.table containing SNP annotations and test statistics}
#' }
#' @source Simulated data generated for LDscnR manuscript
"sim_ex"
