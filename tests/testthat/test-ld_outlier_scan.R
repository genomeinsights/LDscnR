## Same real, seeded stage1 fixture as test-ld_outlier_test.R.
build_real_stage1 <- function(seed = 123) {
  skip_if_not_installed("SNPRelate")
  data(sim_ex, package = "LDscnR")
  map <- data.table::as.data.table(sim_ex$map)[, .(Chr, Pos, marker)]
  GTs <- sim_ex$GTs
  colnames(GTs) <- map$marker
  set.seed(seed)
  gds_path <- tempfile(fileext = ".gds")
  gds <- create_gds_from_geno(sim_ex$GTs, sim_ex$map, gds_path)
  on.exit({ SNPRelate::snpgdsClose(gds); unlink(gds_path) })
  ld <- compute_LD_decay(gds, n_win_decay = 5, max_SNPs_decay = 2000, slide = 1000,
                         keep_el = TRUE, cores = 1)
  stage1 <- ld_complexity_reduction(map, ld, rho = 0.5, cores = 1)
  list(GTs = GTs, map = map, LD_decay = ld, stage1 = stage1)
}

p_obs_for_cluster <- function(d, cl_id, lo = 1e-8, hi = 0.8) {
  mk <- d$stage1$clusters[CL_id == cl_id, members][[1]]
  p <- stats::setNames(rep(hi, nrow(d$map)), d$map$marker)
  p[mk] <- lo
  p
}

test_that("ld_outlier_scan() with only p_obs matches calling ld_outlier_test() directly", {
  d <- build_real_stage1()
  p <- p_obs_for_cluster(d, cl_id = 8)

  scan <- ld_outlier_scan(d$stage1, d$map, p, statistic = "simes", size_floor = 8L,
                          assembly = "physical", verbose = FALSE)
  direct <- ld_outlier_test(d$stage1, d$map, p, statistic = "simes", size_floor = 8L,
                            assembly = "physical")

  expect_s3_class(scan, "ld_outlier_scan")
  expect_identical(scan$test$units, direct$units)
  expect_identical(scan$test$regions, direct$regions)
  expect_null(scan$null)
  expect_null(scan$rotation)
})

test_that("supplying p_perm adds $null, matching a direct ld_outlier_perm() call", {
  d <- build_real_stage1()
  p <- p_obs_for_cluster(d, cl_id = 8)
  set.seed(1)
  p_perm <- matrix(stats::runif(nrow(d$map) * 10, 0.5, 1), nrow(d$map), 10,
                   dimnames = list(d$map$marker, NULL))

  scan <- ld_outlier_scan(d$stage1, d$map, p, p_perm = p_perm, statistic = "simes",
                          size_floor = 8L, assembly = "physical", perm_level = "units",
                          verbose = FALSE)
  expect_s3_class(scan$null, "ld_outlier_perm")

  direct_test <- ld_outlier_test(d$stage1, d$map, p, statistic = "simes",
                                 size_floor = 8L, assembly = "physical")
  direct_null <- ld_outlier_perm(direct_test, d$stage1, d$map, p_perm,
                                 level = "units", verbose = FALSE)
  expect_identical(scan$null$surrogates, direct_null$surrogates)
})

test_that("supplying annotation/chrom_lengths adds $rotation, matching a direct ld_region_rotation() call", {
  d <- build_real_stage1()
  p <- p_obs_for_cluster(d, cl_id = 8)
  chrom_lengths <- data.table::data.table(chr = c("Chr1", "Chr2", "Chr3"), len = c(2e6, 2e6, 2e6))
  annotation <- data.table::data.table(chr = "Chr1", start = 234000, end = 236000)

  scan <- ld_outlier_scan(d$stage1, d$map, p, statistic = "simes", size_floor = 8L,
                          assembly = "physical", annotation = annotation,
                          chrom_lengths = chrom_lengths, n_rotations = 200, verbose = FALSE)
  expect_s3_class(scan$rotation, "ld_region_rotation")

  direct_test <- ld_outlier_test(d$stage1, d$map, p, statistic = "simes",
                                 size_floor = 8L, assembly = "physical")
  direct_rot <- ld_region_rotation(direct_test$regions, annotation, chrom_lengths,
                                   n_rotations = 200, seed = 1L)
  expect_identical(scan$rotation$null_mean, direct_rot$null_mean)
  expect_identical(scan$rotation$p, direct_rot$p)
})

test_that("annotation without chrom_lengths errors before any computation", {
  d <- build_real_stage1()
  p <- p_obs_for_cluster(d, cl_id = 8)
  annotation <- data.table::data.table(chr = "Chr1", start = 1, end = 2)
  expect_error(
    ld_outlier_scan(d$stage1, d$map, p, statistic = "simes", size_floor = 8L,
                    assembly = "physical", annotation = annotation, verbose = FALSE),
    "chrom_lengths"
  )
})
