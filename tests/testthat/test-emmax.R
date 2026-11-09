# Tests for emmax(): mixed-model association F statistics.

make_case <- function(seed = 42, n = 60, m = 40, causal_beta = 0.9) {
  set.seed(seed)
  X <- matrix(stats::rbinom(n * m, 2, 0.3), n, m)
  # marker 1 is causal; rest are noise
  Y <- X[, 1] * causal_beta + stats::rnorm(n)
  K <- X %*% t(X)
  K <- K / mean(diag(K))
  list(X = X, Y = Y, K = K)
}

test_that("emmax returns well-formed output and detects a causal marker", {
  d <- make_case()
  out <- emmax(d$Y, d$X, d$K)

  expect_named(out, c("F", "pval", "Rsq"), ignore.order = TRUE)
  expect_length(out$F, ncol(d$X))
  expect_length(out$pval, ncol(d$X))

  # F non-negative, p-values in [0, 1]
  expect_true(all(out$F >= 0))
  expect_true(all(out$pval >= 0 & out$pval <= 1))

  # the causal marker (column 1) is the strongest, by a clear margin
  expect_equal(which.max(out$F), 1L)
  expect_gt(out$F[1], max(out$F[-1]))
})

test_that("emmax F and p-value are internally consistent (F-test)", {
  d <- make_case()
  out <- emmax(d$Y, d$X, d$K)
  n <- length(d$Y)
  expect_equal(out$pval, stats::pf(out$F, 1, n - 2, lower.tail = FALSE))
})

test_that("a null phenotype yields no strong association", {
  d <- make_case()
  set.seed(7)
  y_null <- stats::rnorm(length(d$Y))          # unrelated to genotypes
  out <- emmax(y_null, d$X, d$K)
  # no marker should look strongly associated
  expect_lt(max(out$F), 25)
  expect_gt(min(out$pval), 0)
})

# ---- emmax_setup(covar =) / emmax_fast(): exact GLS covariates ---------------

## A kinship from background markers that are not scanned, a small set of
## covariates, and genotypes, some of which correlate with the first covariate.
make_cov_case <- function(seed = 1, n = 120, m = 60) {
  set.seed(seed)
  bg <- matrix(stats::rbinom(n * 400, 2, 0.3), n)
  bg <- scale(bg)
  K <- tcrossprod(bg) / ncol(bg) + diag(1e-3, n)
  C <- cbind(c1 = stats::rnorm(n), c2 = stats::rbinom(n, 1, 0.4))
  G <- matrix(stats::rbinom(n * m, 2, 0.3), n, m,
              dimnames = list(NULL, paste0("s", seq_len(m))))
  G[, 1:10] <- G[, 1:10] + round(C[, 1])           # correlated with c1
  g <- as.numeric(t(chol(0.5 * K)) %*% stats::rnorm(n))
  Y <- 0.4 * G[, 11] + 0.8 * C[, 1] - 0.5 * C[, 2] + g + stats::rnorm(n)
  list(G = G, K = K, C = C, Y = Y, n = n)
}

## emmax_fast() as it was before covariates were supported (0.9.2, 936ea2e).
emmax_fast_old <- function(prep, Y) {
  re  <- emma.REMLE(Y, prep$Xo, prep$Kn, eig.R = prep$eigR)
  wv  <- 1 / sqrt(re$vg * prep$lam + re$ve)
  yt  <- as.numeric(crossprod(prep$V, Y)) * wv
  xo  <- prep$xot * wv
  Xtw <- prep$Xt * wv
  a2  <- sum(xo * xo)
  ra  <- yt - xo * (sum(xo * yt) / a2)
  RSS0 <- sum(ra * ra)
  Rb  <- Xtw - outer(xo, as.numeric(crossprod(xo, Xtw)) / a2)
  RSSf <- RSS0 - as.numeric(crossprod(ra, Rb))^2 / colSums(Rb * Rb)
  stats::pf((RSS0 / RSSf - 1) * prep$df2, 1, prep$df2, lower.tail = FALSE)
}

## The reference implementation from the analyses (EpiSel tc_ld_3sp), which
## refits REML with X = [1, C] and projects the whitened X out of y and G.
emmax_fast_cov_ref <- function(prep, y, cov) {
  X  <- cbind(1, cov); re <- emma.REMLE(y, X, prep$Kn)
  wv <- 1 / sqrt(re$vg * prep$lam + re$ve)
  yt <- as.numeric(crossprod(prep$V, y)) * wv; Q <- qr.Q(qr(crossprod(prep$V, X) * wv))
  ra <- yt - Q %*% crossprod(Q, yt); RSS0 <- sum(ra^2)
  Rb <- prep$Xt * wv; Rb <- Rb - Q %*% crossprod(Q, Rb)
  df2 <- nrow(X) - ncol(X) - 1
  RSSf <- RSS0 - as.numeric(crossprod(ra, Rb))^2 / colSums(Rb * Rb)
  stats::pf((RSS0 / RSSf - 1) * df2, 1, df2, lower.tail = FALSE)
}

## |t| recovered from an F(1, df2) p-value.
abs_t <- function(p, df2) sqrt(stats::qf(p, 1, df2, lower.tail = FALSE))

test_that("covar = NULL reproduces the intercept-only scan exactly", {
  d <- make_cov_case()
  prep <- emmax_setup(d$G, d$K)
  p <- emmax_fast(prep, d$Y)
  expect_identical(as.numeric(p), as.numeric(emmax_fast_old(prep, d$Y)))
  expect_identical(names(p), colnames(d$G))
  expect_equal(prep$df2, d$n - 2L)
  expect_equal(attr(p, "df2"), d$n - 2L)
  expect_identical(attr(p, "covar_names"), character(0))
})

test_that("covariate scan matches a direct GLS fit at the REML variance ratio", {
  d <- make_cov_case()
  prep <- emmax_setup(d$G, d$K, covar = d$C)
  p <- emmax_fast(prep, d$Y)
  expect_equal(attr(p, "df2"), d$n - 4L)
  expect_identical(attr(p, "covar_names"), c("c1", "c2"))

  ## oracle 1: the reference implementation used for the published numbers
  expect_equal(as.numeric(p), as.numeric(emmax_fast_cov_ref(prep, d$Y, d$C)),
               tolerance = 1e-8)

  ## oracle 2: lm() on V^(-1/2)-transformed data, V = vg K + ve I from REML
  X <- cbind(1, d$C)
  re <- emma.REMLE(d$Y, X, prep$Kn)
  ev <- eigen(re$vg * prep$Kn + re$ve * diag(d$n), symmetric = TRUE)
  W <- ev$vectors %*% (t(ev$vectors) / sqrt(ev$values))
  ys <- W %*% d$Y; Xs <- W %*% X
  for (j in c(1, 5, 11, 30)) {
    xs <- W %*% d$G[, j]
    fit <- summary(stats::lm(ys ~ Xs + xs - 1))$coefficients["xs", ]
    expect_equal(p[[j]], fit[["Pr(>|t|)"]], tolerance = 1e-8)
    expect_equal(abs_t(p[[j]], d$n - 4L), abs(fit[["t value"]]), tolerance = 1e-8)
  }
})

test_that("REML with X = [1, C] agrees with an independently coded REML", {
  d <- make_cov_case()
  prep <- emmax_setup(d$G, d$K, covar = d$C)
  X <- cbind(1, d$C); n <- d$n; q <- ncol(X)
  reml_ll <- function(logdelta) {
    H <- prep$Kn + exp(logdelta) * diag(n)
    Hi <- solve(H); XHX <- crossprod(X, Hi %*% X)
    P <- Hi - Hi %*% X %*% solve(XHX, crossprod(X, Hi))
    -0.5 * (as.numeric(determinant(H)$modulus) + as.numeric(determinant(XHX)$modulus) +
              (n - q) * log(as.numeric(crossprod(d$Y, P %*% d$Y))))
  }
  opt <- stats::optimize(reml_ll, c(-10, 10), maximum = TRUE, tol = 1e-10)
  re <- emma.REMLE(d$Y, X, prep$Kn, eig.R = prep$eigR)
  expect_equal(log(re$delta), opt$maximum, tolerance = 1e-4)
})

test_that("results are invariant to affine recoding of the covariates", {
  d <- make_cov_case()
  p1 <- emmax_fast(emmax_setup(d$G, d$K, d$C), d$Y)
  C2 <- cbind(a = 3 * d$C[, 1] + 5, b = -2 * d$C[, 2] + d$C[, 1] + 1)
  p2 <- emmax_fast(emmax_setup(d$G, d$K, C2), d$Y)
  expect_equal(as.numeric(p1), as.numeric(p2), tolerance = 1e-8)

  ## a factor gives the same model as its indicator column
  df <- data.frame(c1 = d$C[, 1], grp = factor(ifelse(d$C[, 2] == 1, "B", "A")))
  p3 <- emmax_fast(emmax_setup(d$G, d$K, df), d$Y)
  expect_equal(as.numeric(p1), as.numeric(p3), tolerance = 1e-8)
  expect_identical(attr(p3, "covar_names"), c("c1", "grpB"))
})

test_that("rank-deficient or malformed covariates give an informative error", {
  d <- make_cov_case()
  expect_error(emmax_setup(d$G, d$K, cbind(d$C, dup = 2 * d$C[, 1])),
               "rank-deficient.*dup")
  expect_error(emmax_setup(d$G, d$K, cbind(d$C, k = 7)), "rank-deficient.*k")
  expect_error(emmax_setup(d$G, d$K, d$C[-1, ]), "rows but there are")
  C_na <- d$C; C_na[3, 1] <- NA
  expect_error(emmax_setup(d$G, d$K, C_na), "NA")
})

test_that("conditioning on a correlated covariate gains power over residualising y", {
  set.seed(11)
  n <- 200; m <- 60; r2 <- 0.6
  bg <- scale(matrix(stats::rbinom(n * 400, 2, 0.3), n))
  K <- tcrossprod(bg) / ncol(bg) + diag(1e-3, n)
  C <- stats::rnorm(n)
  Z <- matrix(stats::rnorm(n * m), n)
  G <- sqrt(r2) * C + sqrt(1 - r2) * Z                    # R^2(x ~ C) = 0.6
  colnames(G) <- paste0("s", seq_len(m))
  p_cov <- emmax_setup(G, K, covar = C); p_int <- emmax_setup(G, K)
  L <- t(chol(0.3 * K))
  t_exact <- t_yonly <- numeric(m)
  for (j in seq_len(m)) {
    y <- 0.3 * G[, j] + C + as.numeric(L %*% stats::rnorm(n)) + stats::rnorm(n)
    t_exact[j] <- abs_t(emmax_fast(p_cov, y)[[j]], p_cov$df2)
    t_yonly[j] <- abs_t(emmax_fast(p_int, stats::resid(stats::lm(y ~ C)))[[j]], p_int$df2)
  }
  ratio <- mean(t_exact) / mean(t_yonly)
  expect_gt(ratio, 1.3)                                   # ~ 1 / sqrt(1 - 0.6) = 1.58
  expect_lt(ratio, 1.9)
})

test_that("covariate scan is calibrated under the null", {
  set.seed(5)
  n <- 150; m <- 2000
  bg <- scale(matrix(stats::rbinom(n * 400, 2, 0.3), n))
  K <- tcrossprod(bg) / ncol(bg) + diag(1e-3, n)
  C <- cbind(stats::rnorm(n), stats::rbinom(n, 1, 0.5))
  G <- matrix(stats::rbinom(n * m, 2, 0.3), n, m) + round(0.7 * C[, 1])
  y <- as.numeric(C %*% c(1, -0.7)) + as.numeric(t(chol(0.5 * K)) %*% stats::rnorm(n)) +
    stats::rnorm(n)
  p <- emmax_fast(emmax_setup(G, K, covar = C), y)
  expect_gt(stats::ks.test(as.numeric(p), "punif")$p.value, 0.01)
})

test_that("a matrix of phenotypes equals one call per column", {
  d <- make_cov_case()
  prep <- emmax_setup(d$G, d$K, covar = d$C)
  set.seed(3)
  Y <- cbind(d$Y, replicate(3, sample(d$Y)))
  P <- emmax_fast(prep, Y)
  expect_equal(dim(P), c(ncol(d$G), 4L))
  for (b in seq_len(ncol(Y)))
    expect_identical(unname(P[, b]), unname(as.numeric(emmax_fast(prep, Y[, b]))))
  expect_equal(attr(P, "df2"), d$n - 4L)
})

test_that("emmax(Covar =) fits the covariates exactly, like emmax_fast()", {
  d <- make_cov_case()
  out <- emmax(d$Y, d$G, d$K, Covar = d$C)
  expect_equal(out$pval, stats::pf(out$F, 1, d$n - 4L, lower.tail = FALSE))
  p <- emmax_fast(emmax_setup(d$G, d$K, d$C), d$Y)
  expect_equal(unname(out$pval), as.numeric(p), tolerance = 1e-6)
  ## a constant covariate is ignored, as before
  expect_equal(emmax(d$Y, d$G, d$K, Covar = rep(1, d$n))$pval,
               emmax(d$Y, d$G, d$K)$pval)
})
