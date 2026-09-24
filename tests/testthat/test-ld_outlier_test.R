## Real (not hand-built) stage1, from the bundled sim_ex dataset, the same
## way test-ld_complexity_reduction.R builds it -- deterministic because
## compute_LD_decay()'s internal subsampling reads the global RNG and has no
## seed of its own; set.seed() here pins it. `map` deliberately carries only
## marker/Chr/Pos (no ld_w_* column), the normal ld_outlier_test() input and
## the exact case the hidden ld_w_095 dependency broke.
build_real_stage1 <- function(seed = 123) {
  skip_if_not_installed("SNPRelate")
  data(sim_ex, package = "LDscnR")
  map <- data.table::as.data.table(sim_ex$map)[, .(Chr, Pos, marker)]
  ## sim_ex$GTs's own colnames ("1:5570") are a different naming convention
  ## from sim_ex$map$marker ("Chr1:5570") for the SAME markers in the SAME
  ## column/row order -- fine for create_gds_from_geno(), which only relies
  ## on position, but ld_unit_matrix()/ld_outlier_test() index GTs BY the
  ## marker names stage1's clusters carry (map$marker's convention), so GTs
  ## must be renamed to match before use here.
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

## All-high marker p-values except one cluster's members, set very low --
## Simes on that cluster survives BH; every other cluster does not.
p_obs_for_cluster <- function(d, cl_id, lo = 1e-8, hi = 0.8) {
  mk <- d$stage1$clusters[CL_id == cl_id, members][[1]]
  p <- stats::setNames(rep(hi, nrow(d$map)), d$map$marker)
  p[mk] <- lo
  p
}

test_that("stage1 from a map with only marker/Chr/Pos works under assembly = 'stage2_discovered' (the ld_w_095 regression)", {
  d <- build_real_stage1()
  expect_false(any(grepl("^ld_w", names(d$stage1$map_snp))))  # sanity: fixture has no ld_w column

  p <- p_obs_for_cluster(d, cl_id = 8)
  ## stage2_discovered's ld_prune_and_eMLG() pass prints its own progress --
  ## expect_no_error(), not expect_silent(), is the right bar here
  out <- expect_no_error(
    ld_outlier_test(d$stage1, d$map, p, statistic = "simes", size_floor = 8L,
                    assembly = "stage2_discovered", GTs = d$GTs, LD_decay = d$LD_decay)
  )
  expect_s3_class(out, "ld_outlier_test")
  expect_true(sum(out$units$significant) >= 1)
  expect_true(nrow(out$regions) >= 1)
})

test_that("no discoveries: uniformly high p-values give zero significant units and an empty, well-formed regions table", {
  d <- build_real_stage1()
  p <- stats::setNames(rep(0.9, nrow(d$map)), d$map$marker)
  out <- ld_outlier_test(d$stage1, d$map, p, statistic = "simes", size_floor = 8L,
                         assembly = "stage2_discovered", GTs = d$GTs, LD_decay = d$LD_decay)
  expect_identical(sum(out$units$significant), 0L)
  expect_identical(nrow(out$regions), 0L)
  expect_named(out$regions, c("Chr", "from", "to", "n_units", "n_markers", "occupancy"))

  out_phys <- ld_outlier_test(d$stage1, d$map, p, statistic = "simes", size_floor = 8L,
                              assembly = "physical")
  expect_identical(nrow(out_phys$regions), 0L)
})

test_that("one discovery produces exactly one region under both assembly rules", {
  d <- build_real_stage1()
  p <- p_obs_for_cluster(d, cl_id = 8)

  out_s2 <- ld_outlier_test(d$stage1, d$map, p, statistic = "simes", size_floor = 8L,
                            assembly = "stage2_discovered", GTs = d$GTs, LD_decay = d$LD_decay)
  expect_identical(sum(out_s2$units$significant), 1L)
  expect_identical(nrow(out_s2$regions), 1L)

  out_phys <- ld_outlier_test(d$stage1, d$map, p, statistic = "simes", size_floor = 8L,
                              assembly = "physical")
  expect_identical(nrow(out_phys$regions), 1L)
})

test_that("several discoveries on different chromosomes give at least that many regions under both assemblies", {
  d <- build_real_stage1()
  mk8   <- d$stage1$clusters[CL_id == 8, members][[1]]     # Chr1
  mk124 <- d$stage1$clusters[CL_id == 124, members][[1]]   # Chr2
  mk257 <- d$stage1$clusters[CL_id == 257, members][[1]]   # Chr3
  p <- stats::setNames(rep(0.8, nrow(d$map)), d$map$marker)
  p[c(mk8, mk124, mk257)] <- 1e-8

  out_s2 <- ld_outlier_test(d$stage1, d$map, p, statistic = "simes", size_floor = 8L,
                            assembly = "stage2_discovered", GTs = d$GTs, LD_decay = d$LD_decay)
  expect_identical(sum(out_s2$units$significant), 3L)
  expect_identical(nrow(out_s2$regions), 3L)  # three widely separated chromosomes -> never merge
  expect_setequal(out_s2$regions$Chr, c("Chr1", "Chr2", "Chr3"))

  out_phys <- ld_outlier_test(d$stage1, d$map, p, statistic = "simes", size_floor = 8L,
                              assembly = "physical")
  expect_identical(nrow(out_phys$regions), 3L)
})

test_that("statistic = 'unit' takes pre-computed per-unit p-values and matches Simes' significance call at a controlled extreme", {
  d <- build_real_stage1()
  units <- LDscnR:::.ld_outlier_units(d$stage1, d$map, 8L)
  p_unit <- stats::setNames(rep(0.8, nrow(units)), units$unit_id)
  target <- units[Chr == "Chr1"][1, unit_id]
  p_unit[as.character(target)] <- 1e-8

  out <- ld_outlier_test(d$stage1, d$map, p_unit, statistic = "unit", size_floor = 8L,
                         assembly = "physical")
  expect_identical(sum(out$units$significant), 1L)
  expect_identical(out$units[significant == TRUE, unit_id], target)
})

test_that("a named p_obs (statistic = 'simes') in the wrong order errors instead of silently mis-scoring", {
  d <- build_real_stage1()
  p <- p_obs_for_cluster(d, cl_id = 8)
  p_shuffled <- p[sample(names(p))]  # same names, same values, different order -> names no longer == map$marker
  expect_error(
    ld_outlier_test(d$stage1, d$map, p_shuffled, statistic = "simes", size_floor = 8L,
                    assembly = "physical"),
    "do not match `map\\$marker`"
  )
  ## an UNNAMED vector in map$marker order is the documented, still-legal contract
  expect_silent(
    ld_outlier_test(d$stage1, d$map, unname(p), statistic = "simes", size_floor = 8L,
                    assembly = "physical")
  )
})

test_that("a named p_obs (statistic = 'unit') in the wrong order errors instead of silently mis-scoring", {
  d <- build_real_stage1()
  units <- LDscnR:::.ld_outlier_units(d$stage1, d$map, 8L)
  p_unit <- stats::setNames(rep(0.8, nrow(units)), units$unit_id)
  p_shuffled <- p_unit[sample(names(p_unit))]
  expect_error(
    ld_outlier_test(d$stage1, d$map, p_shuffled, statistic = "unit", size_floor = 8L,
                    assembly = "physical"),
    "do not match the tested units"
  )
  expect_silent(
    ld_outlier_test(d$stage1, d$map, unname(p_unit), statistic = "unit", size_floor = 8L,
                    assembly = "physical")
  )
})

test_that("a right-length but wrong-order named p_obs is caught even when it happens to keep some names in place", {
  ## a weaker shuffle than full sample() -- swap just two entries -- still
  ## must be caught, confirming the check is an exact order comparison, not
  ## a coarser "are the same names present" set check
  d <- build_real_stage1()
  p <- p_obs_for_cluster(d, cl_id = 8)
  nm <- names(p)
  nm[c(1, 2)] <- nm[c(2, 1)]
  names(p) <- nm
  expect_error(
    ld_outlier_test(d$stage1, d$map, p, statistic = "simes", size_floor = 8L, assembly = "physical"),
    "do not match `map\\$marker`"
  )
})

test_that("duplicate map markers are rejected before any expensive computation", {
  d <- build_real_stage1()
  bad_map <- data.table::copy(d$map)
  bad_map$marker[2] <- bad_map$marker[1]
  p <- stats::setNames(rep(0.8, nrow(d$map)), d$map$marker)
  expect_error(
    ld_outlier_test(d$stage1, bad_map, unname(p), statistic = "simes", size_floor = 8L,
                    assembly = "physical"),
    "duplicate"
  )
})

test_that("markers referenced by stage1 but missing from map are rejected with an informative error", {
  d <- build_real_stage1()
  ## must drop a marker belonging to a cluster that clears size_floor -- a
  ## dropped singleton/small-cluster marker would never be looked at by
  ## .ld_outlier_units() in the first place (filtered out before the check).
  victim <- d$stage1$clusters[CL_id == 8, members][[1]][1]
  drop_i <- which(d$map$marker == victim)
  bad_map <- d$map[-drop_i]
  p <- stats::setNames(rep(0.8, nrow(d$map)), d$map$marker)
  expect_error(
    ld_outlier_test(d$stage1, bad_map, unname(p)[-drop_i], statistic = "simes", size_floor = 8L,
                    assembly = "physical"),
    "missing from"
  )
})

test_that("markers in the significant clusters missing from GTs are rejected before ld_prune_and_eMLG()", {
  d <- build_real_stage1()
  p <- p_obs_for_cluster(d, cl_id = 8)
  bad_GTs <- d$GTs[, colnames(d$GTs) != d$stage1$clusters[CL_id == 8, members][[1]][1]]
  expect_error(
    ld_outlier_test(d$stage1, d$map, p, statistic = "simes", size_floor = 8L,
                    assembly = "stage2_discovered", GTs = bad_GTs, LD_decay = d$LD_decay),
    "missing from"
  )
})

test_that("p_obs of the wrong length errors for both statistics", {
  d <- build_real_stage1()
  expect_error(
    ld_outlier_test(d$stage1, d$map, runif(3), statistic = "simes", size_floor = 8L, assembly = "physical"),
    "statistic = \"simes\""
  )
  units <- LDscnR:::.ld_outlier_units(d$stage1, d$map, 8L)
  expect_error(
    ld_outlier_test(d$stage1, d$map, runif(nrow(units) + 1), statistic = "unit", size_floor = 8L,
                    assembly = "physical"),
    "statistic = \"unit\""
  )
})
