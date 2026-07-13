# Regression tests for the LD-decay core workflow.
# These exercise the retained public surface end-to-end on the bundled
# `sim_ex` dataset: create_gds_from_geno() -> compute_LD_decay() -> helpers.

# Build an ld_decay object with small, fast settings. The temporary GDS file
# is closed and removed before returning, so callers get a self-contained object.
build_decay <- function(keep_el = FALSE) {
  skip_if_not_installed("SNPRelate")
  data(sim_ex, package = "LDscnR")

  gds_path <- tempfile(fileext = ".gds")
  gds <- create_gds_from_geno(sim_ex$GTs, sim_ex$map, gds_path)
  on.exit({
    SNPRelate::snpgdsClose(gds)
    unlink(gds_path)
  })

  compute_LD_decay(
    gds,
    n_win_decay    = 5,
    max_SNPs_decay = 2000,
    slide          = 1000,
    keep_el        = keep_el,
    cores          = 1
  )
}

test_that("compute_LD_decay returns a well-formed ld_decay object", {
  ld <- build_decay()

  expect_s3_class(ld, "ld_decay")
  expect_named(
    ld,
    c("by_chr", "decay_sum", "decay_model", "recommendation", "params"),
    ignore.order = TRUE
  )

  # one row of decay parameters per chromosome
  data(sim_ex, package = "LDscnR")
  n_chr <- length(unique(sim_ex$map$Chr))
  expect_equal(nrow(ld$decay_sum), n_chr)

  # core decay parameters are present and finite
  expect_true(all(c("a", "b", "c") %in% names(ld$decay_sum)))
  expect_true(all(is.finite(ld$decay_sum$a)))
})

test_that("print.ld_decay runs without error", {
  ld <- build_decay()
  expect_output(print(ld), "ld_decay")
})

test_that("keep_el stores edge lists per chromosome", {
  ld <- build_decay(keep_el = TRUE)
  expect_false(is.null(ld$by_chr[[1]]$el))
})

test_that("d_from_rho and ld_from_rho behave monotonically", {
  a <- 1e-5

  # larger rho => larger physical window
  expect_gt(d_from_rho(a, 0.99), d_from_rho(a, 0.90))
  expect_equal(d_from_rho(a, 0.90), 0.90 / (a * (1 - 0.90)))

  # higher rho = farther / more decayed, so expected r^2 falls toward
  # background b (rho -> 1) and rises toward short-range c (rho -> 0)
  b <- 0.1
  cc <- 1
  expect_lt(ld_from_rho(b, cc, 0.99), ld_from_rho(b, cc, 0.90))
  expect_lt(ld_from_rho(b, cc, 0.50), cc)
  expect_gt(ld_from_rho(b, cc, 0.50), b)
})
