#' Pre-compute the EMMAX mixed-model rotation for fast repeated scans
#'
#' Eigendecomposes the kinship once and rotates the genotype matrix into the
#' whitening basis, so a genome-wide EMMAX scan for a new phenotype
#' ([emmax_fast()]) reduces to a per-SNP whitened F test. Numerically identical
#' to [emmax()] (per-SNP p-values correlate ~1) but roughly 25x faster per
#' phenotype. Use it when many phenotypes are scanned against the same
#' genotypes + kinship, e.g. the surrogate phenotypes of a structure-aware null
#' ([structured_null()]).
#'
#' @param GTs Numeric genotype dosage matrix, individuals in rows and SNPs in
#'   columns (column names = markers).
#' @param K Kinship / GRM matrix (individuals x individuals), aligned to `GTs`.
#'
#' @return An `emmax_prep` list consumed by [emmax_fast()]: the scaled kinship,
#'   its eigenbasis, the rotated genotypes, and the residual degrees of freedom.
#' @seealso [emmax_fast()], [emmax()]
#' @export
emmax_setup <- function(GTs, K) {
  n <- nrow(GTs)
  one_n <- matrix(1 / n, n, n)
  Kn <- (n - 1) / sum((diag(n) - one_n) * K) * K        # scale kinship to match emmax()
  Xo <- matrix(1, n, 1)
  ev <- eigen(Kn, symmetric = TRUE)
  list(Kn = Kn, Xo = Xo, eigR = emma.eigen.R.wo.Z(Kn, Xo),
       V = ev$vectors, lam = ev$values,
       Xt = crossprod(ev$vectors, GTs), xot = as.numeric(crossprod(ev$vectors, Xo)),
       df2 = n - 2L)
}

#' Fast EMMAX association scan for one phenotype
#'
#' Genome-wide EMMAX p-values for phenotype `y`, using a rotation pre-computed by
#' [emmax_setup()]. Estimates the variance components by REML once, whitens
#' response and genotypes, and returns per-SNP F-test p-values (df 1 and n-2).
#' Identical to the p-values from [emmax()] at a fraction of the cost.
#'
#' @param prep An `emmax_prep` object from [emmax_setup()].
#' @param y Numeric phenotype vector, length = number of individuals (a 0/1
#'   ecotype/case vector is treated as a continuous response, as in EMMAX).
#'
#' @return Numeric vector of per-SNP p-values, one per column of the genotype
#'   matrix passed to [emmax_setup()].
#' @seealso [emmax_setup()]
#' @export
emmax_fast <- function(prep, y) {
  re  <- emma.REMLE(y, prep$Xo, prep$Kn, eig.R = prep$eigR)
  wv  <- 1 / sqrt(re$vg * prep$lam + re$ve)              # whitening weights
  yt  <- as.numeric(crossprod(prep$V, y)) * wv
  xo  <- prep$xot * wv
  Xtw <- prep$Xt * wv
  a2  <- sum(xo * xo)
  ra  <- yt - xo * (sum(xo * yt) / a2)                  # residual under the null (intercept only)
  RSS0 <- sum(ra * ra)
  Rb  <- Xtw - outer(xo, as.numeric(crossprod(xo, Xtw)) / a2)
  RSSf <- RSS0 - as.numeric(crossprod(ra, Rb))^2 / colSums(Rb * Rb)
  stats::pf((RSS0 / RSSf - 1) * prep$df2, 1, prep$df2, lower.tail = FALSE)
}
