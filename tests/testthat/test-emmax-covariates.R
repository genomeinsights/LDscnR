# Tests for exact fixed-effect covariates in emmax_setup()/emmax_fast()/emmax():
# GLS conditional on X = [1, covar], REML with that X, df2 = n - ncol(X) - 1.

## A kinship from background markers that are not scanned, two covariates, and
## genotypes, the first ten of which correlate with the first covariate.
make_cov_case <- function(seed = 1, n = 120, m = 60) {
  set.seed(seed)
  bg <- scale(matrix(stats::rbinom(n * 400, 2, 0.3), n))
  K <- tcrossprod(bg) / ncol(bg) + diag(1e-3, n)
  C <- cbind(c1 = stats::rnorm(n), c2 = stats::rbinom(n, 1, 0.4))
  G <- matrix(stats::rbinom(n * m, 2, 0.3), n, m,
              dimnames = list(NULL, paste0("s", seq_len(m))))
  G[, 1:10] <- G[, 1:10] + round(C[, 1])
  g <- as.numeric(t(chol(0.5 * K)) %*% stats::rnorm(n))
  Y <- 0.4 * G[, 11] + 0.8 * C[, 1] - 0.5 * C[, 2] + g + stats::rnorm(n)
  list(G = G, K = K, C = C, Y = Y, n = n)
}

## The reference implementation from the analyses (EpiSel tc_ld_3sp/03_tc_tests.R),
## which refits REML with X = [1, C] and projects the whitened X out of y and G.
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

## V^(-1/2) at the REML variance components for design X.
whitener <- function(prep, Y, X) {
  re <- emma.REMLE(Y, X, prep$Kn)
  ev <- eigen(re$vg * prep$Kn + re$ve * diag(nrow(X)), symmetric = TRUE)
  ev$vectors %*% (t(ev$vectors) / sqrt(ev$values))
}

## |t| recovered from an F(1, df2) p-value.
abs_t <- function(p, df2) sqrt(stats::qf(p, 1, df2, lower.tail = FALSE))

test_that("1. no covariates reproduces the stored 0.9.2 p-values bit for bit", {
  fx <- readRDS(test_path("fixtures", "emmax_0.9.2.rds"))
  prep <- emmax_setup(fx$G, fx$K)
  p <- emmax_fast(prep, fx$Y)
  expect_identical(as.numeric(p), as.numeric(fx$p_fast))
  expect_identical(names(p), names(fx$p_fast))
  expect_equal(prep$df2, nrow(fx$G) - 2L)
  expect_equal(attr(p, "df2"), nrow(fx$G) - 2L)
  expect_identical(attr(p, "covar_names"), character(0))
})

test_that("2. covariate scan matches an explicit GLS fit at the REML variance components", {
  d <- make_cov_case()
  prep <- emmax_setup(d$G, d$K, covar = d$C)
  p <- emmax_fast(prep, d$Y)
  expect_equal(attr(p, "df2"), d$n - 4L)
  expect_identical(attr(p, "covar_names"), c("c1", "c2"))

  ## the reference implementation behind the published conditioning numbers
  expect_equal(as.numeric(p), as.numeric(emmax_fast_cov_ref(prep, d$Y, d$C)),
               tolerance = 1e-8)

  ## lm(y* ~ 0 + X* + x*) on whitened data
  X <- cbind(1, d$C)
  W <- whitener(prep, d$Y, X)
  ys <- W %*% d$Y; Xs <- W %*% X
  for (j in c(1, 5, 11, 30)) {
    xs <- W %*% d$G[, j]
    fit <- summary(stats::lm(ys ~ 0 + Xs + xs))$coefficients["xs", ]
    expect_equal(p[[j]], fit[["Pr(>|t|)"]], tolerance = 1e-8)
    expect_equal(abs_t(p[[j]], d$n - 4L), abs(fit[["t value"]]), tolerance = 1e-8)
  }
})

test_that("3. FWL: equals the test of residualised y* on residualised x*", {
  d <- make_cov_case()
  prep <- emmax_setup(d$G, d$K, covar = d$C)
  p <- emmax_fast(prep, d$Y)
  X <- cbind(1, d$C); df2 <- d$n - ncol(X) - 1
  W <- whitener(prep, d$Y, X)
  Xs <- W %*% X
  ry <- stats::resid(stats::lm(W %*% d$Y ~ 0 + Xs))
  for (j in c(2, 11, 40)) {
    rx <- stats::resid(stats::lm(W %*% d$G[, j] ~ 0 + Xs))
    b <- sum(rx * ry) / sum(rx * rx)
    rss <- sum((ry - b * rx)^2)
    Fj <- b^2 * sum(rx * rx) / (rss / df2)
    expect_equal(p[[j]], stats::pf(Fj, 1, df2, lower.tail = FALSE), tolerance = 1e-8)
  }
})

test_that("4. a 4-level factor covariate costs 3 degrees of freedom", {
  d <- make_cov_case()
  set.seed(9)
  f <- factor(sample(c("a", "b", "c", "d"), d$n, replace = TRUE))
  prep <- emmax_setup(d$G, d$K, covar = f)
  expect_equal(prep$df2, d$n - 5L)
  expect_identical(prep$covar_names, c("covarb", "covarc", "covard"))
  p <- emmax_fast(prep, d$Y)
  expect_equal(attr(p, "df2"), d$n - 5L)
  expect_equal(as.numeric(p), as.numeric(emmax_fast_cov_ref(prep, d$Y, model.matrix(~ f)[, -1])),
               tolerance = 1e-8)
})

test_that("5. invariance to intercept-like and affine covariates; rank deficiency errors", {
  d <- make_cov_case()
  p1 <- emmax_fast(emmax_setup(d$G, d$K, d$C), d$Y)

  ## a covariate that is a linear function of the intercept is dropped
  p_const <- emmax_fast(emmax_setup(d$G, d$K, cbind(d$C, k = 7)), d$Y)
  expect_identical(as.numeric(p_const), as.numeric(p1))
  expect_identical(attr(p_const, "covar_names"), c("c1", "c2"))
  p_none <- emmax_fast(emmax_setup(d$G, d$K, rep(3, d$n)), d$Y)
  expect_identical(as.numeric(p_none), as.numeric(emmax_fast(emmax_setup(d$G, d$K), d$Y)))

  ## affine recoding of C spans the same column space
  C2 <- cbind(a = 3 * d$C[, 1] + 5, b = -2 * d$C[, 2] + d$C[, 1] + 1)
  expect_equal(as.numeric(emmax_fast(emmax_setup(d$G, d$K, C2), d$Y)), as.numeric(p1),
               tolerance = 1e-8)

  ## a data.frame with a factor gives the same model as its indicator column
  df <- data.frame(c1 = d$C[, 1], grp = ifelse(d$C[, 2] == 1, "B", "A"))
  p3 <- emmax_fast(emmax_setup(d$G, d$K, df), d$Y)
  expect_equal(as.numeric(p3), as.numeric(p1), tolerance = 1e-8)
  expect_identical(attr(p3, "covar_names"), c("c1", "grpB"))

  expect_error(emmax_setup(d$G, d$K, cbind(d$C, dup = 2 * d$C[, 1] - 1)),
               "rank-deficient.*dup")
  expect_error(emmax_setup(d$G, d$K, d$C[-1, ]), "rows but there are")
  C_na <- d$C; C_na[3, 1] <- NA
  expect_error(emmax_setup(d$G, d$K, C_na), "NA")
  expect_error(emmax_setup(d$G[1:4, ], d$K[1:4, 1:4], cbind(1:4, c(1, 3, 2, 5))),
               "degrees of freedom")
})

test_that("5b. a SNP the covariates explain completely, or a monomorphic one, gives NA", {
  d <- make_cov_case()
  G <- cbind(d$G, mono = 1, same = d$C[, 2])
  p <- emmax_fast(emmax_setup(G, d$K, d$C), d$Y)
  expect_true(is.na(p[["mono"]]) && !is.nan(p[["mono"]]))
  expect_true(is.na(p[["same"]]))
  expect_false(anyNA(p[colnames(d$G)]))
  p0 <- emmax_fast(emmax_setup(cbind(d$G, mono = 0), d$K), d$Y)
  expect_false(is.nan(p0[["mono"]]))
})

test_that("6. covar_method = 'residualise' and the Covar alias reproduce old emmax() exactly", {
  fx <- readRDS(test_path("fixtures", "emmax_0.9.2.rds"))
  old <- fx$emmax_Covar
  new <- emmax(fx$Y, fx$G, fx$K, covar = fx$covar, covar_method = "residualise")
  expect_identical(new$F, old$F)
  expect_identical(new$pval, old$pval)
  expect_identical(new$Rsq, old$Rsq)

  .emmax_state$covar_warned <- NULL
  expect_warning(alias <- emmax(fx$Y, fx$G, fx$K, Covar = fx$covar), "deprecated")
  expect_identical(alias$pval, old$pval)
  expect_silent(emmax(fx$Y, fx$G, fx$K, Covar = fx$covar))   # warned once per session
  expect_error(emmax(fx$Y, fx$G, fx$K, Covar = fx$covar, covar = fx$covar), "only")
})

test_that("6b. emmax(covar =) defaults to exact GLS, the same model as emmax_fast()", {
  d <- make_cov_case()
  out <- emmax(d$Y, d$G, d$K, covar = d$C)
  expect_equal(out$pval, stats::pf(out$F, 1, d$n - 4L, lower.tail = FALSE))
  p <- emmax_fast(emmax_setup(d$G, d$K, d$C), d$Y)
  expect_equal(unname(out$pval), as.numeric(p), tolerance = 1e-6)
  ## without covariates the two methods are the same test
  expect_identical(emmax(d$Y, d$G, d$K)$pval,
                   emmax(d$Y, d$G, d$K, covar_method = "residualise")$pval)
})

test_that("7. permuting y with C fixed is calibrated for SNPs unrelated to y", {
  skip_on_cran()
  set.seed(21)
  n <- 100; m <- 80
  bg <- scale(matrix(stats::rbinom(n * 400, 2, 0.3), n))
  K <- tcrossprod(bg) / ncol(bg) + diag(1e-3, n)
  C <- stats::rnorm(n)
  G <- matrix(stats::rbinom(n * m, 2, 0.3), n, m)
  G[, 1:20] <- G[, 1:20] + round(1.5 * C)                  # block strongly aligned with C
  y <- C + as.numeric(t(chol(0.5 * K)) %*% stats::rnorm(n)) + stats::rnorm(n)
  prep <- emmax_setup(G, K, covar = C)
  P <- emmax_fast(prep, replicate(200, sample(y)))
  expect_equal(dim(P), c(m, 200L))
  rate <- mean(P[21:m, ] < 0.05)
  expect_gt(rate, 0.035)
  expect_lt(rate, 0.065)
  expect_lt(abs(mean(P[1:20, ] < 0.05) - 0.05), 0.02)
})

test_that("conditioning on a correlated covariate gains power over residualising y", {
  set.seed(11)
  n <- 200; m <- 60; r2 <- 0.6
  bg <- scale(matrix(stats::rbinom(n * 400, 2, 0.3), n))
  K <- tcrossprod(bg) / ncol(bg) + diag(1e-3, n)
  C <- stats::rnorm(n)
  G <- sqrt(r2) * C + sqrt(1 - r2) * matrix(stats::rnorm(n * m), n)   # R^2(x ~ C) = 0.6
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

test_that("a single null phenotype gives uniform p-values over ~2,000 SNPs", {
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

test_that("a prep built before covariates were supported still scans", {
  d <- make_cov_case()
  prep <- emmax_setup(d$G, d$K)
  old <- prep[c("Kn", "Xo", "eigR", "V", "lam", "Xt", "xot", "df2")]
  expect_identical(as.numeric(emmax_fast(old, d$Y)), as.numeric(emmax_fast(prep, d$Y)))
})

test_that("structured_null() surrogates are orthogonal to [1, covar, y]", {
  d <- make_cov_case(m = 20)
  prep <- emmax_setup(d$G, d$K, covar = d$C)
  Xy <- cbind(prep$X, d$Y)
  set.seed(4)
  n <- d$n; eK <- eigen(d$K, symmetric = TRUE); Lv <- pmax(eK$values, 0)
  s <- as.numeric(eK$vectors %*% (sqrt(Lv) * stats::rnorm(n)))
  r <- as.numeric(qr.resid(qr(Xy), s))
  expect_lt(max(abs(crossprod(Xy, r))), 1e-8)
  ## and the function itself runs with covariates
  ld_ws <- matrix(0.5, ncol(d$G), 2, dimnames = list(colnames(d$G), c("0.5", "0.9")))
  nl <- structured_null(d$Y, d$G, d$K, ld_ws, B = 3, prep = prep)
  expect_s3_class(nl, "ld_null")
})

test_that("8. 9sp lineage: exact GLS reproduces 52/20 consensus and 24/14 Simes", {
  skip_on_cran()
  bundle <- Sys.getenv("LDSCNR_9SP_BUNDLE")
  scan <- Sys.getenv("LDSCNR_9SP_SCAN")
  skip_if(!nzchar(bundle) || !file.exists(bundle) || !file.exists(scan),
          "set LDSCNR_9SP_BUNDLE and LDSCNR_9SP_SCAN to module_9sp/out/{02_bundle/bundle,03_EMMAX/scan}.rds")
  b <- readRDS(bundle); sc <- readRDS(scan)
  lin <- factor(b$pheno$lineage)
  P <- sc$consensus$test$params
  tst <- function(p, stat) ld_outlier_test(b$stage1, b$map, p, statistic = stat,
    size_floor = P$size_floor, alpha = P$alpha, assembly = P$assembly, GTs = b$GTs,
    LD_decay = b$LD_decay, score_threshold = P$score_threshold,
    distance_threshold = P$distance_threshold)
  um <- ld_unit_matrix(b$GTs, b$stage1, b$map, size_floor = P$size_floor,
                       repr = "consensus_dosage")
  pu <- emmax_fast(emmax_setup(um, b$GRM, covar = lin), b$eco)
  tu <- tst(pu, "unit")
  expect_equal(sum(tu$units$significant), 52)
  expect_equal(nrow(tu$regions), 20)
  rm(um); gc()
  pm <- emmax_fast(emmax_setup(b$GTs, b$GRM, covar = lin), b$eco)
  tm <- tst(pm, "simes")
  expect_equal(sum(tm$units$significant), 24)
  expect_equal(nrow(tm$regions), 14)
})
