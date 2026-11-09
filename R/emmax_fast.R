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
#' Covariates enter the model exactly: the fixed effects are \eqn{X = [1, C]},
#' REML is fitted with that \eqn{X}, and each SNP is tested by generalised least
#' squares conditional on \eqn{X} (Frisch-Waugh-Lovell in the \eqn{V^{-1}}
#' metric), with `n - ncol(X) - 1` residual degrees of freedom. This is not the
#' same as regressing the phenotype on `covar` first and scanning the residuals:
#' that leaves the genotypes unadjusted, which shrinks the effect of every SNP
#' that is correlated with a covariate and loses power for exactly those SNPs.
#'
#' @param GTs Numeric genotype dosage matrix, individuals in rows and SNPs in
#'   columns (column names = markers).
#' @param K Kinship / GRM matrix (individuals x individuals), aligned to `GTs`.
#' @param covar Optional covariates, aligned to the rows of `GTs`: a numeric
#'   vector, a factor, a numeric matrix, or a data.frame. It is turned into a
#'   design matrix once, with `model.matrix(~ ., data = as.data.frame(covar))`
#'   without its intercept column, so factor and character columns become
#'   indicator columns. Columns with zero variance are dropped. The intercept
#'   is always added, so do not include one. Together with the intercept the
#'   remaining columns must have full rank, and `covar` may not contain `NA`. A
#'   prep object holds one covariate set, reused for every phenotype passed to
#'   [emmax_fast()]; build another prep for a different set. `NULL` (default)
#'   fits an intercept only, exactly as before covariates were supported.
#'
#' @details The model is \eqn{y = X\beta + x\gamma + g + e}, with
#'   \eqn{g \sim N(0, \sigma_g^2 K_n)} and \eqn{e \sim N(0, \sigma_e^2 I)}, where
#'   \eqn{X = [1, C]} and \eqn{K_n} is `K` rescaled as in [emmax()]. The variance
#'   components are estimated by REML under the null (\eqn{\gamma = 0}) with
#'   that \eqn{X}, once per phenotype.
#'
#' @return An `emmax_prep` list consumed by [emmax_fast()]: the scaled kinship
#'   `Kn`, its eigenbasis (`V`, `lam`), the rotated genotypes `Xt`, the design
#'   `X` (`[1, covar]`) and its rotation `Xt_X`, the REML eigensystem `eigR`
#'   for that design, the residual degrees of freedom `df2`
#'   (`n - ncol(X) - 1`) and `covar_names`. `Xo` and `xot` (the intercept and
#'   its rotation) are kept for code written against earlier versions.
#' @seealso [emmax_fast()], [emmax()]
#' @export
emmax_setup <- function(GTs, K, covar = NULL) {
  n <- nrow(GTs)
  one_n <- matrix(1 / n, n, n)
  Kn <- (n - 1) / sum((diag(n) - one_n) * K) * K        # scale kinship to match emmax()
  Xo <- matrix(1, n, 1)
  X <- .emmax_design(covar, n)
  ev <- eigen(Kn, symmetric = TRUE)
  list(Kn = Kn, Xo = Xo, eigR = emma.eigen.R.wo.Z(Kn, X),
       V = ev$vectors, lam = ev$values,
       Xt = crossprod(ev$vectors, GTs), xot = as.numeric(crossprod(ev$vectors, Xo)),
       df2 = n - ncol(X) - 1L,
       X = X, Xt_X = crossprod(ev$vectors, X), covar_names = colnames(X)[-1L])
}

## The fixed-effect design [1, covar] for emmax_setup() and emmax(). Factors
## and character columns become indicator columns and zero-variance columns are
## dropped; the result must have full column rank, because REML and the GLS
## projection both assume it.
.emmax_design <- function(covar, n) {
  X <- matrix(1, n, 1L, dimnames = list(NULL, "(Intercept)"))
  if (is.null(covar)) return(X)
  df <- as.data.frame(covar)
  if (is.null(dim(covar)) && ncol(df) == 1L) names(df) <- "covar"
  if (nrow(df) != n)
    stop("`covar` has ", nrow(df), " rows but there are ", n, " individuals.")
  if (anyNA(df)) stop("`covar` must not contain NA.")
  C <- stats::model.matrix(~ ., data = df)[, -1L, drop = FALSE]
  C <- C[, apply(C, 2L, function(v) any(v != v[1L])), drop = FALSE]
  X <- cbind(X, C)
  qx <- qr(X)
  if (qx$rank < ncol(X)) {
    bad <- colnames(X)[qx$pivot[seq.int(qx$rank + 1L, ncol(X))]]
    stop("`covar` is rank-deficient together with the intercept: ",
         paste(bad, collapse = ", "), " is a linear combination of the other ",
         "columns. Drop it.")
  }
  if (n - ncol(X) - 1L < 1L)
    stop("Too many covariates for ", n, " individuals: no residual degrees of freedom left.")
  X
}

#' Fast EMMAX association scan for one phenotype
#'
#' Genome-wide EMMAX p-values for phenotype `Y`, using a rotation pre-computed by
#' [emmax_setup()]. Estimates the variance components by REML once, whitens
#' response and genotypes, and returns per-SNP F-test p-values (df 1 and
#' `prep$df2`, which is `n - 2` without covariates). Identical to the p-values
#' from [emmax()] at a fraction of the cost.
#'
#' With covariates (`emmax_setup(covar = )`), each SNP is tested conditional on
#' them by exact GLS: the whitened design \eqn{X^*} is projected out of the
#' whitened phenotype once, and per SNP only \eqn{\|Q'G^*_j\|^2} is added to
#' the intercept-only cost, where \eqn{Q} is an orthonormal basis of
#' \eqn{X^*}.
#'
#' @param prep An `emmax_prep` object from [emmax_setup()].
#' @param Y Numeric phenotype vector, length = number of individuals (a 0/1
#'   ecotype/case vector is treated as a continuous response, as in EMMAX). A
#'   matrix with one phenotype per column (e.g. permuted phenotypes) is scanned
#'   column by column, each with its own REML fit.
#'
#' @return For a vector `Y`, a numeric vector of per-SNP p-values, one per
#'   column of the genotype matrix passed to [emmax_setup()]. For a matrix `Y`,
#'   a SNPs x `ncol(Y)` matrix. Either carries the attributes `df2` (residual
#'   degrees of freedom) and `covar_names` (the covariate columns of the model,
#'   `character(0)` without covariates).
#' @seealso [emmax_setup()]
#' @param y Deprecated. Use `Y`.
#' @export
emmax_fast <- function(prep, Y,
                       y = NULL) {
  if (!missing(y)) Y <- .renamed_arg(y, "y", "Y")
  covar_names <- if (is.null(prep$covar_names)) character(0) else prep$covar_names
  p <- if (is.matrix(Y)) {
    out <- vapply(seq_len(ncol(Y)), function(b) .emmax_fast1(prep, Y[, b]),
                  numeric(ncol(prep$Xt)))
    dimnames(out) <- list(colnames(prep$Xt), colnames(Y))
    out
  } else .emmax_fast1(prep, Y)
  attr(p, "df2") <- prep$df2
  attr(p, "covar_names") <- covar_names
  p
}

## One phenotype. An intercept-only prep keeps the original arithmetic, so its
## p-values are bit-for-bit what they were before covariates were supported. A
## SNP with no variation left to test (monomorphic, or explained completely by
## the covariates) gets NA.
.emmax_fast1 <- function(prep, Y) {
  X   <- if (is.null(prep$X)) prep$Xo else prep$X     # preps from before 0.9.3 have no X
  re  <- emma.REMLE(Y, X, prep$Kn, eig.R = prep$eigR)
  wv  <- 1 / sqrt(re$vg * prep$lam + re$ve)              # whitening weights
  yt  <- as.numeric(crossprod(prep$V, Y)) * wv
  Xtw <- prep$Xt * wv
  if (ncol(X) == 1L) {
    xo  <- prep$xot * wv
    a2  <- sum(xo * xo)
    ra  <- yt - xo * (sum(xo * yt) / a2)                # residual under the null (intercept only)
    RSS0 <- sum(ra * ra)
    Rb  <- Xtw - outer(xo, as.numeric(crossprod(xo, Xtw)) / a2)
    RSSf <- RSS0 - as.numeric(crossprod(ra, Rb))^2 / colSums(Rb * Rb)
    p <- stats::pf((RSS0 / RSSf - 1) * prep$df2, 1, prep$df2, lower.tail = FALSE)
    p[is.nan(p)] <- NA_real_
    return(p)
  }
  ## Exact GLS conditional on X = [1, covar] (Frisch-Waugh-Lovell): with Q an
  ## orthonormal basis of the whitened X and M = I - QQ', the SNP effect is
  ## G'My / G'MG. M is idempotent, so the numerator needs only My, and
  ## G'MG = ||G||^2 - ||Q'G||^2.
  Q   <- qr.Q(qr(prep$Xt_X * wv))
  ra  <- yt - as.numeric(Q %*% crossprod(Q, yt))
  RSS0 <- sum(ra * ra)
  num <- as.numeric(crossprod(ra, Xtw))
  den <- colSums(Xtw * Xtw) - colSums(crossprod(Q, Xtw)^2)
  den[den <= 1e-10 * colSums(Xtw * Xtw)] <- NA_real_
  RSSf <- RSS0 - num^2 / den
  p <- stats::pf((RSS0 / RSSf - 1) * prep$df2, 1, prep$df2, lower.tail = FALSE)
  p[is.nan(p)] <- NA_real_
  p
}
