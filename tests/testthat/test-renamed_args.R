## 0.9.2 renamed arguments so that one parameter has one name across the
## package. The old names keep working with a deprecation warning.

test_that("old argument names still work and warn", {
  ds <- data.table::data.table(Chr = "1", a = 1e-4, b = 0.02, c = 0.9)

  new <- score_thresholds(ds, r2_rho = 0.6, dist_rho = 0.9)
  expect_warning(old <- score_thresholds(ds, rho_r2 = 0.6), "`rho_r2` is deprecated")
  expect_equal(old$r2min, new$r2min)
  expect_warning(old <- score_thresholds(ds, rho_d = 0.9), "`rho_d` is deprecated")
  expect_equal(old$dmax, new$dmax)

  set.seed(1)
  p <- runif(200)
  w <- runif(200)
  expect_warning(old <- rmsc_threshold(p, w, fdr = 0.1), "`fdr` is deprecated")
  expect_equal(old, rmsc_threshold(p, w, alpha = 0.1))
})

test_that("the new names give the same result as the defaults they replaced", {
  ds <- data.table::data.table(Chr = "1", a = 1e-4, b = 0.02, c = 0.9)
  expect_silent(score_thresholds(ds))
  expect_equal(score_thresholds(ds),
               suppressWarnings(score_thresholds(ds, rho_r2 = 0.75, rho_d = 0.95)))
})

test_that("get_el() and compute_LD_decay() share the \"corr\" default", {
  expect_identical(formals(get_el)$ld_method, "corr")
  expect_identical(formals(compute_LD_decay)$ld_method, "corr")
})
