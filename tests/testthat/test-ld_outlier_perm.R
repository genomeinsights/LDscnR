## Same real, seeded stage1 fixture as test-ld_outlier_test.R (kept local to
## this file, matching this repo's convention of self-contained test files).
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

## B all-null (uniformly high) surrogates -- under any reasonable null
## construction the true positive discovery in `obs` should stand out.
build_obs_and_perm <- function(d, B = 20) {
  p <- p_obs_for_cluster(d, cl_id = 8)
  obs <- ld_outlier_test(d$stage1, d$map, p, statistic = "simes", size_floor = 8L,
                         assembly = "physical")
  set.seed(1)
  p_perm <- matrix(stats::runif(nrow(d$map) * B, 0.5, 1), nrow(d$map), B,
                   dimnames = list(d$map$marker, NULL))
  list(d = d, obs = obs, p_perm = p_perm, B = B)
}

test_that("level = 'units' and level = 'regions' both run and count sensibly, and need not agree", {
  x <- build_obs_and_perm(build_real_stage1())
  null_u <- ld_outlier_perm(x$obs, x$d$stage1, x$d$map, x$p_perm, B = x$B,
                            level = "units", verbose = FALSE)
  null_r <- ld_outlier_perm(x$obs, x$d$stage1, x$d$map, x$p_perm,
                            GTs = x$d$GTs, LD_decay = x$d$LD_decay, B = x$B,
                            level = "regions", verbose = FALSE)
  expect_identical(null_u$observed, sum(x$obs$units$significant))
  expect_identical(null_r$observed, nrow(x$obs$regions))
  expect_length(null_u$surrogates, x$B)
  expect_length(null_r$surrogates, x$B)
  expect_true(all(null_u$surrogates >= 0))
  expect_true(all(null_r$surrogates >= 0))
  ## all-null surrogates should essentially never beat one designed true discovery
  expect_true(mean(null_u$surrogates) < null_u$observed)
})

test_that("matrix, list and function forms of p_perm give identical surrogate counts for the same draws", {
  x <- build_obs_and_perm(build_real_stage1())

  null_mat <- ld_outlier_perm(x$obs, x$d$stage1, x$d$map, x$p_perm, B = x$B,
                              level = "units", verbose = FALSE)

  p_list <- lapply(seq_len(x$B), function(b) x$p_perm[, b])
  null_list <- ld_outlier_perm(x$obs, x$d$stage1, x$d$map, p_list, B = x$B,
                               level = "units", verbose = FALSE)

  p_fun <- function(b) x$p_perm[, b]
  null_fun <- ld_outlier_perm(x$obs, x$d$stage1, x$d$map, p_fun, B = x$B,
                              level = "units", verbose = FALSE)

  expect_identical(null_mat$surrogates, null_list$surrogates)
  expect_identical(null_mat$surrogates, null_fun$surrogates)
})

test_that("p_perm as a list infers B from its length without an explicit B", {
  x <- build_obs_and_perm(build_real_stage1())
  p_list <- lapply(seq_len(x$B), function(b) x$p_perm[, b])
  null_list <- ld_outlier_perm(x$obs, x$d$stage1, x$d$map, p_list, level = "units", verbose = FALSE)
  expect_length(null_list$surrogates, x$B)
})

test_that("a function form of p_perm requires an explicit B", {
  x <- build_obs_and_perm(build_real_stage1())
  p_fun <- function(b) x$p_perm[, b]
  expect_error(
    ld_outlier_perm(x$obs, x$d$stage1, x$d$map, p_fun, level = "units", verbose = FALSE),
    "`B`"
  )
})

test_that("an unrecognised p_perm shape errors clearly", {
  x <- build_obs_and_perm(build_real_stage1())
  expect_error(
    ld_outlier_perm(x$obs, x$d$stage1, x$d$map, "not a valid shape", level = "units", verbose = FALSE),
    "matrix.*list.*function"
  )
})

test_that("a named surrogate column in the wrong order is caught the same way p_obs is", {
  x <- build_obs_and_perm(build_real_stage1())
  shuffled <- sample(nrow(x$d$map))
  p_bad <- stats::setNames(x$p_perm[shuffled, 1], x$d$map$marker[shuffled])  # names shuffled too, still self-consistent but != map$marker order
  expect_error(
    ld_outlier_perm(x$obs, x$d$stage1, x$d$map, list(p_bad), B = 1L, level = "units", verbose = FALSE),
    "do not match `map\\$marker`"
  )
})
