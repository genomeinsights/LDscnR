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
