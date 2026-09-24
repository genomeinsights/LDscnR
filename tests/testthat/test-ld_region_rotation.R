## Fixture: three regions on two chromosomes of very different length. Region
## 2 sits right at the end of chr1 (from=100, to=150, chr1 len=200) and region
## 3's span (55) exactly equals chr2's length (55), the tightest possible fit.
make_regions <- function() {
  data.table::data.table(chr = c("1", "1", "2"), from = c(10, 100, 5), to = c(20, 150, 60))
}
make_len <- function() data.table::data.table(chr = c("1", "2"), len = c(200, 55))
make_annotation <- function() data.table::data.table(chr = "1", start = 15, end = 16)

test_that("within: relocated regions never exceed their chromosome's bounds", {
  R <- make_regions()
  len_of <- c("1" = 200, "2" = 55)
  relocate <- LDscnR:::.build_relocator(R, len_of, "within")
  set.seed(1)
  draws <- data.table::rbindlist(lapply(1:2000, function(i) relocate()))
  expect_true(all(draws$from >= 0))
  expect_true(all(draws[chr == "1", to] <= 200))
  expect_true(all(draws[chr == "2", to] <= 55))
  ## span-preserving: every draw's width equals the corresponding region's
  ## observed width (recycled every 3 rows, one per region, in row order)
  expect_equal(draws$to - draws$from, rep(R$to - R$from, 2000))
})

test_that("genome: chromosome is drawn proportional to placement space, not uniformly", {
  ## 400 equal-span regions, one long and one short chromosome; the
  ## proportion landing on the long chromosome should track
  ## (len_long - w) / ((len_long - w) + (len_short - w)), not 0.5.
  len_of <- c("1" = 1000, "2" = 100)
  w <- 5
  R <- data.table::data.table(chr = rep("1", 400), from = 1, to = 1 + w)
  relocate <- LDscnR:::.build_relocator(R, len_of, "genome")
  set.seed(7)
  draws <- data.table::rbindlist(lapply(1:300, function(i) relocate()))
  observed_frac <- mean(draws$chr == "1")
  expected_frac <- (1000 - w) / ((1000 - w) + (100 - w))
  expect_equal(observed_frac, expected_frac, tolerance = 0.02)
  expect_false(isTRUE(all.equal(observed_frac, 0.5, tolerance = 0.02)))
})

test_that("genome: relocated regions never exceed their assigned chromosome's bounds", {
  R <- make_regions()
  len_of <- c("1" = 200, "2" = 55)
  relocate <- LDscnR:::.build_relocator(R, len_of, "genome")
  set.seed(3)
  draws <- data.table::rbindlist(lapply(1:2000, function(i) relocate()))
  draws[, len := len_of[chr]]
  expect_true(all(draws$from >= 0))
  expect_true(all(draws$to <= draws$len))
  expect_equal(draws$to - draws$from, rep(R$to - R$from, 2000))
})

test_that("a region longer than its own chromosome errors under 'within'", {
  R_bad <- data.table::data.table(chr = "2", from = 0, to = 100)  # chr2 len = 55
  expect_error(
    ld_region_rotation(R_bad, make_annotation(), make_len(), scheme = "within", n_rotations = 10),
    "only 55 bp long"
  )
})

test_that("a region longer than every chromosome errors under 'genome'", {
  R_bad <- data.table::data.table(chr = "1", from = 0, to = 250)  # longer than both chr1=200, chr2=55
  expect_error(
    ld_region_rotation(R_bad, make_annotation(), make_len(), scheme = "genome", n_rotations = 10),
    "does not fit on any chromosome"
  )
})

test_that("a region exactly as long as a chromosome is placeable, not an error", {
  R_tight <- data.table::data.table(chr = "2", from = 0, to = 55)  # chr2 len = 55, exact fit
  expect_silent(
    out <- ld_region_rotation(R_tight, make_annotation(), make_len(), scheme = "within", n_rotations = 50)
  )
  expect_s3_class(out, "ld_region_rotation")
})

test_that("zero regions is handled without error and gives a well-defined result", {
  R0 <- make_regions()[0]
  for (sc in c("within", "genome")) {
    out <- ld_region_rotation(R0, make_annotation(), make_len(), scheme = sc, n_rotations = 200, seed = 1)
    expect_identical(out$observed, 0L)
    expect_identical(out$null_mean, 0)
    expect_identical(out$fold, 0)
    expect_identical(out$p, 1)
    expect_equal(nrow(out$per_region), 0L)
  }
})

test_that("results are reproducible for a fixed seed", {
  R <- make_regions()
  o1 <- ld_region_rotation(R, make_annotation(), make_len(), scheme = "genome", n_rotations = 500, seed = 42)
  o2 <- ld_region_rotation(R, make_annotation(), make_len(), scheme = "genome", n_rotations = 500, seed = 42)
  expect_identical(o1$null_mean, o2$null_mean)
  expect_identical(o1$p, o2$p)
  o3 <- ld_region_rotation(R, make_annotation(), make_len(), scheme = "genome", n_rotations = 500, seed = 43)
  expect_false(identical(o1$null_mean, o3$null_mean))
})

test_that("a region overlapping several annotation intervals is counted once", {
  R <- data.table::data.table(chr = "1", from = 0, to = 100)
  A <- data.table::data.table(chr = "1", start = c(10, 20, 30), end = c(15, 25, 35))
  L <- data.table::data.table(chr = "1", len = 1000)
  out <- ld_region_rotation(R, A, L, scheme = "within", n_rotations = 10, seed = 1)
  expect_identical(out$observed, 1L)
  expect_identical(sum(out$per_region$on_peak), 1L)
})

test_that("one and several discoveries both produce sensible fold/p", {
  L <- data.table::data.table(chr = "1", len = 1000)
  A <- data.table::data.table(chr = "1", start = c(50, 500), end = c(60, 510))

  R_one <- data.table::data.table(chr = "1", from = 45, to = 65)
  out_one <- ld_region_rotation(R_one, A, L, scheme = "within", n_rotations = 2000, seed = 1)
  expect_identical(out_one$observed, 1L)
  expect_true(out_one$p >= 0 && out_one$p <= 1)

  R_several <- data.table::data.table(chr = "1", from = c(45, 495), to = c(65, 515))
  out_several <- ld_region_rotation(R_several, A, L, scheme = "within", n_rotations = 2000, seed = 1)
  expect_identical(out_several$observed, 2L)
  expect_true(out_several$p >= 0 && out_several$p <= 1)
})

test_that("print method reports observed, fold and p without erroring", {
  out <- ld_region_rotation(make_regions(), make_annotation(), make_len(), n_rotations = 50, seed = 1)
  expect_output(print(out), "observed 1/3")
})
