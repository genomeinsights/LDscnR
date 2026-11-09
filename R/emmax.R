## EMMAX association test + EMMA REML machinery.
##
## Original implementation by Zitong Li (lizitong1985@gmail.com), from:
##   Li Z, Kemppainen P, Rastas P, Merila J. "Linkage disequilibrium
##   clustering-based approach for association mapping with tightly linked
##   genome-wide data." Molecular Ecology Resources (2018).
## Included in LDscnR with permission as the downstream association step for the
## LD-pruned marker set; only namespace-qualification of stats/MASS/parallel
## calls has been added -- the numerical algorithm is unchanged, except for
## emmax(covar_method = "gls"), the exact covariate model added in 0.9.3 (the
## original covariate handling is covar_method = "residualise"). The emma.*
## helpers are the EMMA REML routines (Kang et al. 2008) and are internal.

## Session state: the `Covar` deprecation warning is given once per session.
.emmax_state <- new.env(parent = emptyenv())

#' EMMAX single-locus association test accounting for relatedness
#'
#' Efficient mixed-model association (EMMAX) for a quantitative or binary
#' phenotype, correcting for genome-wide relatedness through a kinship matrix
#' `K`. This is a natural downstream use of the LD-pruned representative marker
#' set produced by [ld_prune_and_eMLG()]: build `K` from the pruned markers
#' (e.g. `SNPRelate::snpgdsGRM`) and then test individual markers, cluster
#' representatives, or eMLG consensus dosages for association.
#'
#' @param Y Numeric phenotype vector (length `n` individuals).
#' @param X Numeric genotype matrix, individuals (rows) by markers/units
#'   (columns), coded as dosages.
#' @param K Pre-computed `n` by `n` kinship / genomic relationship matrix.
#' @param B Optional number of permutations for a permutation-adjusted p-value.
#'   `NULL` (default) skips permutation; if set, the `MASS` package is required.
#' @param binary Logical; if `TRUE` and `B` is set, permutation phenotypes are
#'   dichotomised to match the case/control ratio of `Y`.
#' @param Covar Deprecated. Use `covar`. Supplying `Covar` selects
#'   `covar_method = "residualise"` (unless `covar_method` is given), so scripts
#'   written for earlier versions give the same numbers as before.
#' @param cores Number of cores for the permutation loop.
#' @param covar Optional covariates: a vector, factor, matrix or data.frame
#'   with one row per individual (see [emmax_setup()]'s `covar` for how it is
#'   turned into a design matrix).
#' @param covar_method How `covar` enters the test. `"gls"` (default) fits it
#'   as fixed effects alongside the intercept: REML uses \eqn{X = [1, C]} and
#'   every marker is tested conditional on \eqn{X} by generalised least
#'   squares, on `n - ncol(X) - 1` residual degrees of freedom; this is the
#'   same model as [emmax_setup()] + [emmax_fast()]. `"residualise"` is the
#'   method of earlier versions: `Y` is replaced by the residuals of
#'   `lm(Y ~ covar)` (a covariate that does not vary is ignored), the
#'   genotypes are left unadjusted and the test uses `n - 2` degrees of
#'   freedom. It loses power for markers correlated with a covariate and is
#'   slightly anti-conservative; use it only to reproduce old results.
#'
#' @return A list with `F` (per-column F statistics), `pval` (F-test p-values),
#'   `Rsq` (per-column \eqn{R^2}), and -- when `B` is supplied -- `pval.corr`
#'   (permutation-adjusted p-values).
#'
#' @references
#' Li Z, Kemppainen P, Rastas P, Merila J (2018). Linkage disequilibrium
#' clustering-based approach for association mapping with tightly linked
#' genome-wide data. \emph{Molecular Ecology Resources}. Original implementation
#' by Zitong Li.
#'
#' @export
emmax <- function(Y, X, K, B = NULL, binary = FALSE, Covar = NULL, cores = 1,
                  covar = NULL, covar_method = c("gls", "residualise")) {
  if (!is.null(Covar)) {
    if (!is.null(covar)) stop("Supply `covar` only; `Covar` is its deprecated alias.")
    if (!isTRUE(.emmax_state$covar_warned)) {
      warning("`Covar` is deprecated: it residualises `Y` on the covariates ",
              "(covar_method = \"residualise\"). Use `covar`, which fits them ",
              "exactly (covar_method = \"gls\").", call. = FALSE)
      .emmax_state$covar_warned <- TRUE
    }
    covar <- Covar
    if (missing(covar_method)) covar_method <- "residualise"
  }
  covar_method <- match.arg(covar_method)

  n <- length(Y)
  m <- ncol(X)
  nbchunks <- 2
  stopifnot(ncol(K) == n)
  stopifnot(nrow(K) == n)
  stopifnot(nrow(X) == n)
  stopifnot(nbchunks >= 2)

  ## FIXED EFFECTS: intercept plus covariates, tested exactly ("gls"), or the
  ## intercept only with Y residualised on the covariates ("residualise")
  if (covar_method == "residualise") {
    if (!is.null(covar) & length(unique(covar)) > 1) {
      Y <- stats::resid(stats::lm(Y ~ covar))
    }
    Xo <- matrix(1, n, 1)
  } else {
    Xo <- .emmax_design(covar, n)
  }

  ## K MATRIX NORMALISATION
  K_norm <- (n - 1) / sum((diag(n) - matrix(1, n, n) / n) * K) * K

  ## NULL MODEL + EMMAX association test
  null <- emma.REMLE(Y, Xo, K_norm)

  M <- solve(chol(null$vg * K_norm + null$ve * diag(n)))
  Y_t <- crossprod(M, Y)
  Xo_t <- crossprod(M, Xo)

  RSS <- list()
  for (j in 1:(nbchunks - 1)) {
    X_t <- crossprod(M, X[, ((j - 1) * round(m / nbchunks) + 1):(j * round(m / nbchunks))])
    RSS[[j]] <- apply(X_t, 2, function(x) sum(stats::lsfit(cbind(Xo_t, x), Y_t, intercept = FALSE)$residuals^2))
    rm(X_t)
  }
  X_t <- crossprod(M, X[, ((j) * round(m / nbchunks) + 1):(m)])
  RSS[[nbchunks]] <- apply(X_t, 2, function(x) sum(stats::lsfit(cbind(Xo_t, x), Y_t, intercept = FALSE)$residuals^2))
  rm(X_t, j)

  RSSf <- unlist(RSS)
  RSS_H0 <- rep(sum(stats::lsfit(Xo_t, Y_t, intercept = FALSE)$residuals^2), m)
  df1 <- 1
  df2 <- n - ncol(Xo) - df1
  R2 <- 1 - 1 / (RSS_H0 / RSSf)
  F <- (RSS_H0 / RSSf - 1) * df2 / df1
  pval <- stats::pf(F, df1, df2, lower.tail = FALSE)

  ## optional permutation test
  if (!is.null(B)) {
    if (!requireNamespace("MASS", quietly = TRUE))
      stop("Permutation testing (B) requires the 'MASS' package.")
    Fo <- sort(abs(F))
    Ord <- order(F)
    sigma <- K_norm * null$vg + diag(rep(1, n)) * null$ve
    Ynew <- MASS::mvrnorm(B, rep(0, ncol(sigma)), sigma)

    if (binary) {
      Min <- min(table(Y))
      Ynew <- apply(Ynew, 1, function(x) {
        Ord <- order(x)
        y <- rep(NA, length(Ord))
        y[Ord <= Min] <- 1
        y[Ord > Min] <- 0
        return(y)
      })
      Ynew <- t(Ynew)
    }

    FR <- do.call(rbind, parallel::mclapply(1:B, function(k) {
      y <- Ynew[k, ]
      null <- emma.REMLE(y, Xo, K_norm)
      M <- solve(chol(null$vg * K_norm + null$ve * diag(n)))
      Y_t <- crossprod(M, y)
      Xo_t <- crossprod(M, Xo)
      RSS <- list()
      for (j in 1:(nbchunks - 1)) {
        X_t <- crossprod(M, X[, ((j - 1) * round(m / nbchunks) + 1):(j * round(m / nbchunks))])
        RSS[[j]] <- apply(X_t, 2, function(x) sum(stats::lsfit(cbind(Xo_t, x), Y_t, intercept = FALSE)$residuals^2))
        rm(X_t)
      }
      X_t <- crossprod(M, X[, ((j) * round(m / nbchunks) + 1):(m)])
      RSS[[nbchunks]] <- apply(X_t, 2, function(x) sum(stats::lsfit(cbind(Xo_t, x), Y_t, intercept = FALSE)$residuals^2))
      rm(X_t, j)
      RSSf <- unlist(RSS)
      RSS_H0 <- rep(sum(stats::lsfit(Xo_t, Y_t, intercept = FALSE)$residuals^2), m)
      df1 <- 1
      df2 <- n - ncol(Xo) - df1
      (RSS_H0 / RSSf - 1) * df2 / df1
    }, mc.cores = cores))

    Qmat <- t(apply(FR, 1, cummax))
    Padj <- apply(t(matrix(rep(Fo, B), m)) < Qmat, 2, mean)
    o <- order(Ord)
    out <- Padj[o]
    names(out) <- names(pval)
    list("F" = F, "pval" = pval, "pval.corr" = out, "Rsq" = R2)
  } else {
    list("F" = F, "pval" = pval, "Rsq" = R2)
  }
}

## ---- EMMA REML machinery (Kang et al. 2008); internal ---------------------
## Only the no-Z (Z = NULL) path is retained -- the sole path emmax() exercises.
## The original had an unreachable Z-branch that referenced helpers it never
## defined; dropping it keeps only verified, exercised code.

emma.REMLE <- function(y, X, K, ngrids = 100, llim = -10, ulim = 10,
                       esp = 1e-10, eig.R = NULL) {
  n <- length(y)
  t <- nrow(K)
  q <- ncol(X)

  stopifnot(ncol(K) == t)
  stopifnot(nrow(X) == n)

  if (det(crossprod(X, X)) == 0) {
    warning("X is singular")
    return(list(REML = 0, delta = 0, ve = 0, vg = 0))
  }

  if (is.null(eig.R)) eig.R <- emma.eigen.R.wo.Z(K, X)
  etas <- crossprod(eig.R$vectors, y)

  logdelta <- (0:ngrids) / ngrids * (ulim - llim) + llim
  m <- length(logdelta)
  delta <- exp(logdelta)
  Lambdas <- matrix(eig.R$values, n - q, m) + matrix(delta, n - q, m, byrow = TRUE)
  Etasq <- matrix(etas * etas, n - q, m)
  dLL <- 0.5 * delta * ((n - q) * colSums(Etasq / (Lambdas * Lambdas)) / colSums(Etasq / Lambdas) - colSums(1 / Lambdas))

  optlogdelta <- vector(length = 0)
  optLL <- vector(length = 0)
  if (dLL[1] < esp) {
    optlogdelta <- append(optlogdelta, llim)
    optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim, eig.R$values, etas))
  }
  if (dLL[m - 1] > 0 - esp) {
    optlogdelta <- append(optlogdelta, ulim)
    optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim, eig.R$values, etas))
  }

  for (i in 1:(m - 1)) {
    if ((dLL[i] * dLL[i + 1] < 0 - esp * esp) && (dLL[i] > 0) && (dLL[i + 1] < 0)) {
      r <- stats::uniroot(emma.delta.REML.dLL.wo.Z, lower = logdelta[i], upper = logdelta[i + 1],
                          lambda = eig.R$values, etas = etas)
      optlogdelta <- append(optlogdelta, r$root)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root, eig.R$values, etas))
    }
  }

  maxdelta <- exp(optlogdelta[which.max(optLL)])
  maxLL <- max(optLL)
  maxva <- sum(etas * etas / (eig.R$values + maxdelta)) / (n - q)
  maxve <- maxva * maxdelta

  return(list(REML = maxLL, delta = maxdelta, ve = maxve, vg = maxva))
}

emma.eigen.R.wo.Z <- function(K, X) {
  n <- nrow(X)
  q <- ncol(X)
  S <- diag(n) - X %*% solve(crossprod(X, X)) %*% t(X)
  eig <- eigen(S %*% (K + diag(1, n)) %*% S, symmetric = TRUE)
  stopifnot(!is.complex(eig$values))
  return(list(values = eig$values[1:(n - q)] - 1, vectors = eig$vectors[, 1:(n - q)]))
}

emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <- exp(logdelta)
  etasq <- etas * etas
  ldelta <- lambda + delta
  return(0.5 * (nq * sum(etasq / (ldelta * ldelta)) / sum(etasq / ldelta) - sum(1 / ldelta)))
}

emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <- exp(logdelta)
  return(0.5 * (nq * (log(nq / (2 * pi)) - 1 - log(sum(etas * etas / (lambda + delta)))) - sum(log(lambda + delta))))
}
