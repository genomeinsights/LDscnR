# Regression tests for the LD-decay core workflow.
# These exercise the retained public surface end-to-end on the bundled
# `sim_ex` dataset: create_gds_from_geno() -> compute_LD_decay() -> helpers.

# Build an ld_decay object with small, fast settings. The temporary GDS file
# is closed and removed before returning, so callers get a self-contained object.
build_decay <- function(keep_el = FALSE, seed = NULL, n_sub_bg = 5000) {
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
    cores          = 1,
    seed           = seed,
    n_sub_bg       = n_sub_bg
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

test_that("in-place ld_w is the same whether edges are reused or rebuilt", {
  skip_if_not_installed("SNPRelate")
  data(sim_ex, package = "LDscnR")

  gds_path <- tempfile(fileext = ".gds")
  gds <- create_gds_from_geno(sim_ex$GTs, sim_ex$map, gds_path)
  on.exit({ SNPRelate::snpgdsClose(gds); unlink(gds_path) })

  RHO  <- c(0.5, 0.95)
  args <- list(gds, n_win_decay = 5, max_SNPs_decay = Inf, slide = 1000,
               cores = 1, ld_w_rho = RHO)

  ## same seed: the decay fit subsamples pairs, so the curve (and hence a_pred,
  ## which ld_w is defined against) is only comparable within a seed
  set.seed(7); kept  <- do.call(compute_LD_decay, c(args, keep_el = TRUE))
  set.seed(7); rebui <- do.call(compute_LD_decay, c(args, keep_el = FALSE))

  ## keep_el = FALSE holds no edges: ld_w rebuilt them per chromosome from gds
  expect_null(rebui$by_chr[[1]]$el)
  expect_false(is.null(kept$by_chr[[1]]$el))

  expect_equal(rebui$ld_ws, kept$ld_ws)
  expect_equal(rebui$decay_sum, kept$decay_sum)
  expect_equal(dim(rebui$ld_ws), c(nrow(sim_ex$map), length(RHO)))
})

test_that("a chromosome whose decay fit fails is dropped cleanly, not left as a NULL hole", {
  skip_if_not_installed("SNPRelate")
  data(sim_ex, package = "LDscnR")

  gds_path <- tempfile(fileext = ".gds")
  gds <- create_gds_from_geno(sim_ex$GTs, sim_ex$map, gds_path)
  on.exit({ SNPRelate::snpgdsClose(gds); unlink(gds_path) })

  ## n_win_decay = 1 yields far fewer than the 5 windows summarize_decay() needs,
  ## so every chromosome fails the gate. That used to leave NULL elements in the
  ## pre-allocated by_chr list, which consumers dereferenced into an opaque
  ## data.table error ("RHS of == is length 0") far from the real cause.
  expect_error(
    suppressWarnings(compute_LD_decay(gds, n_win_decay = 1, max_SNPs_decay = Inf,
                                      slide = 1000, cores = 1, keep_el = FALSE)),
    "failed for every chromosome"
  )

  ## and a normal run leaves no NULL entries behind
  ld <- compute_LD_decay(gds, n_win_decay = 5, max_SNPs_decay = Inf,
                         slide = 1000, cores = 1, keep_el = FALSE)
  expect_false(any(vapply(ld$by_chr, is.null, logical(1))))
  expect_setequal(names(ld$by_chr), ld$decay_sum$Chr)
})


test_that("seed makes the fit reproducible, and restores the caller's RNG stream", {
  ## Why this matters beyond tidiness: b is estimated from a subsample and feeds
  ## every decay-relative threshold via ld_from_rho(), so an unseeded pair of
  ## runs differs by enough to move stage-1 cluster counts and GRM sizes by
  ## ~0.7% -- larger than some real settings effects, which makes unseeded A/B
  ## comparisons in that range unsound.
  a <- build_decay(seed = 42)
  b <- build_decay(seed = 42)
  expect_equal(a$decay_sum$b, b$decay_sum$b)
  expect_equal(a$decay_sum$a, b$decay_sum$a)
  expect_equal(a$decay_sum$c, b$decay_sum$c)

  ## A different seed must be a genuinely different draw, or the test above
  ## would pass against a `seed` argument that is silently ignored. That needs
  ## n_sub_bg BELOW the marker count: sim_ex has 900 markers against the default
  ## subsample of 5000, so by default the background step takes everything and
  ## is deterministic whatever the seed. Which is worth knowing in itself --
  ## the seed only bites on data large enough to actually be subsampled.
  d1 <- build_decay(seed = 42, n_sub_bg = 200)
  d2 <- build_decay(seed = 7,  n_sub_bg = 200)
  expect_false(isTRUE(all.equal(d1$decay_sum$b, d2$decay_sum$b)))
  expect_equal(d1$decay_sum$b, build_decay(seed = 42, n_sub_bg = 200)$decay_sum$b)

  ## the caller's stream must be untouched: seeding the decay fit should not
  ## silently re-align later draws in a pipeline that also samples
  set.seed(99); before <- runif(3)
  set.seed(99); invisible(build_decay(seed = 42)); after <- runif(3)
  expect_equal(before, after)
})

test_that("co-located markers do not break the decay stratification", {
  ## Regression: two markers at the same bp give distance 0, log(0) = -Inf, and
  ## the strata breaks become non-finite -- which surfaced far from the cause as
  ## "'from' must be a finite number" from seq(), and made compute_LD_decay()
  ## fail outright on a dense simulated map (95k loci in 0.79 Mb, ~600 such
  ## pairs per chromosome). LDscnR::: because the helper is internal.
  set.seed(1)
  sub <- data.table::data.table(d = c(0, 0, 10^runif(200, 1, 5)), r2 = runif(202))
  expect_silent(res <- suppressMessages(
    LDscnR:::subsample_pairs_for_decay(sub, max_pairs = 50, n_strata = 5)))
  expect_true(nrow(res) > 0)
  expect_true(all(res$d > 0))          # the zero-distance pairs are gone

  ## degenerate case: everything at one distance, so there is nothing to stratify
  flat <- data.table::data.table(d = rep(500, 30), r2 = runif(30))
  expect_silent(r2 <- suppressMessages(
    LDscnR:::subsample_pairs_for_decay(flat, max_pairs = 10, n_strata = 5)))
  expect_equal(nrow(r2), 10)

  ## and nothing usable at all returns empty rather than erroring
  zero <- data.table::data.table(d = c(0, 0, NA_real_), r2 = runif(3))
  expect_equal(nrow(suppressMessages(
    LDscnR:::subsample_pairs_for_decay(zero, max_pairs = 5, n_strata = 5))), 0L)
})
