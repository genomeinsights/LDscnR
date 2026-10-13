#' Simulated example dataset for LDscnR
#'
#' A simulated genome scan dataset including:
#' \itemize{
#'   \item \code{GTs}: genotype matrix (individuals × SNPs)
#'   \item \code{map}: marker annotations, QTN status and test statistics
#' }
#'
#' @format A list with components:
#' \describe{
#'   \item{GTs}{Numeric matrix of genotypes; columns are in the row order of \code{map}}
#'   \item{map}{data.table with one row per marker: \code{Chr}, \code{Pos}, \code{marker},
#'     \code{type} (\code{"QTN"} or \code{"ntrl"}), \code{p_Va} and marker-wise test results}
#' }
#' @source Simulated data generated for LDscnR manuscript
"sim_ex"
