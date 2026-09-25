## Same fixture as test-ld_unit_matrix.R: CL1+CL2 merge into one "F..." group
## covering two Stage-1 clusters; CL3 is correlated but too far away -> its
## own "F..." group; CL4 is unflagged -> passes through as its own "U..."
## group, a 1:1 map onto a single Stage-1 cluster.
build_stage1 <- function() {
  set.seed(1)
  n_ind <- 60
  latent1 <- rbinom(n_ind, 2, 0.35)
  mk <- function(latent, k, noise = 0.05) {
    sapply(seq_len(k), function(i) {
      pmin(pmax(latent + rbinom(n_ind, 1, noise) - rbinom(n_ind, 1, noise), 0), 2)
    })
  }
  CL1 <- mk(latent1, 3); colnames(CL1) <- paste0("cl1_", 1:3)
  CL2 <- mk(latent1, 3); colnames(CL2) <- paste0("cl2_", 1:3)
  CL3 <- mk(latent1, 3); colnames(CL3) <- paste0("cl3_", 1:3)
  CL4 <- matrix(rbinom(n_ind * 2, 2, 0.3), ncol = 2); colnames(CL4) <- paste0("cl4_", 1:2)
  GTs <- cbind(CL1, CL2, CL3, CL4)
  rownames(GTs) <- paste0("ind", seq_len(n_ind))

  map_snp <- data.table::data.table(
    marker = colnames(GTs), Chr = "Chr1",
    Pos = c(1000:1002, 1010:1012, 5000:5002, 2000:2001),
    CL_id = rep(1:4, c(3, 3, 3, 2)),
    n_loci = rep(c(3, 3, 3, 2), c(3, 3, 3, 2)),
    is_core = unlist(list(c(TRUE, FALSE, FALSE), c(TRUE, FALSE, FALSE),
                          c(TRUE, FALSE, FALSE), c(TRUE, FALSE))),
    median_ld = NA_real_,
    ld_w_095 = c(0.9, 0.05, 0.05,  0.9, 0.9, 0.9,  0.9, 0.9, 0.9,  0.05, 0.05)
  )
  clusters <- data.table::data.table(
    Chr = "Chr1", CL_id = 1:4,
    core_snp = c("cl1_1", "cl2_1", "cl3_1", "cl4_1"),
    median_ld = c(0.9, 0.9, 0.9, NA),
    n_snps = c(3, 3, 3, 2),
    members = list(colnames(CL1), colnames(CL2), colnames(CL3), colnames(CL4))
  )
  map <- map_snp[, .(marker, Chr, Pos)]
  list(GTs = GTs, map = map, stage1 = list(map_snp = map_snp, clusters = clusters))
}

build_prune_result <- function(d) {
  ld_prune_and_eMLG(
    GTs = d$GTs, stage1 = d$stage1, ld_w_col = "ld_w_095", ld_w_threshold = 0.5,
    score_threshold = 0.80, min_r2 = 0.2, distance_threshold = 100, cores = 1
  )
}

test_that("ld_group_matrix() errors if fill = FALSE is passed through best_snp_args", {
  d <- build_stage1()
  pr <- build_prune_result(d)
  expect_error(
    ld_group_matrix(d$GTs, pr, d$map, size_floor = 2L, best_snp_args = list(fill = FALSE)),
    "fill = TRUE"
  )
})

test_that("ld_group_matrix() reports prune_result's own groups, not stage1's units, with exact column alignment", {
  d <- build_stage1()
  pr <- build_prune_result(d)
  ## sanity on the fixture: CL1+CL2 merged into one 6-locus group (spans two
  ## Stage-1 clusters), CL4 passed through unflagged (1:1 with Stage-1 CL4)
  merged_id <- pr$groups[startsWith(group_id, "F") & n_loci == 6, group_id]
  expect_length(merged_id, 1)
  unflagged_id <- pr$groups[startsWith(group_id, "U") & n_loci == 2, group_id]
  expect_length(unflagged_id, 1)

  m <- ld_group_matrix(d$GTs, pr, d$map, size_floor = 2L)
  u <- attr(m, "units")

  ## column order matches attr(,"units") row order exactly
  expect_identical(colnames(m), u$unit_id)
  expect_identical(nrow(m), nrow(d$GTs))

  ## unit_id is prune_result's group_id (a character id, e.g. "F1"/"U3"),
  ## not stage1's own integer unit_id -- and its count/identity differs from
  ## what .ld_outlier_units() would build for the same stage1/size_floor,
  ## which is exactly why this cannot be handed to
  ## ld_outlier_test(statistic = "unit") as p_obs (see ld_group_matrix()'s docs).
  expect_type(u$unit_id, "character")
  expect_true(merged_id %in% u$unit_id)
  expect_true(unflagged_id %in% u$unit_id)
  stage1_units <- LDscnR:::.ld_outlier_units(d$stage1, d$map, 2L)
  expect_false(nrow(stage1_units) == nrow(u) &&
               identical(sort(as.character(stage1_units$unit_id)), sort(u$unit_id)))

  ## the merged group's reported span covers BOTH constituent Stage-1
  ## clusters (Pos 1000-1012), not just one of them
  merged_row <- u[unit_id == merged_id]
  expect_equal(merged_row$from, 1000)
  expect_equal(merged_row$to, 1012)
  expect_equal(merged_row$n_markers, 6L)

  ## every reported unit's best_marker is an actual column of GTs, and that
  ## column's genotype matches the returned matrix column exactly
  expect_true(all(u$best_marker %in% colnames(d$GTs)))
  for (id in u$unit_id) {
    bm <- u[unit_id == id, best_marker]
    expect_equal(unname(m[, id]), unname(d$GTs[, bm]), tolerance = 1e-8)
  }
})

test_that("ld_group_matrix()'s size_floor filters prune_result's groups by n_loci, independent of stage1's own filtering", {
  d <- build_stage1()
  pr <- build_prune_result(d)

  ## size_floor = 5 exceeds stage1's largest RAW cluster (3), which would
  ## make .ld_outlier_units() error -- but the merged 6-locus "F" group
  ## still clears it, so ld_group_matrix() must not fail here.
  m <- ld_group_matrix(d$GTs, pr, d$map, size_floor = 5L)
  u <- attr(m, "units")
  expect_true(all(u$n_markers >= 5L))
  expect_true(all(startsWith(u$unit_id, "F")))  # the 2-locus "U" group is filtered out
})

test_that("ld_group_matrix() errors informatively when no prune_result group clears size_floor", {
  d <- build_stage1()
  pr <- build_prune_result(d)
  expect_error(
    ld_group_matrix(d$GTs, pr, d$map, size_floor = 1000L),
    "size_floor"
  )
})

test_that("ld_group_matrix()'s column count need not equal ld_unit_matrix()'s unit count at the same size_floor", {
  ## The concrete failure mode this function's docs warn against: naively
  ## treating its columns as "the units" and feeding them to
  ## ld_outlier_test(statistic = "unit") would use the wrong test count.
  d <- build_stage1()
  pr <- build_prune_result(d)
  m_group <- ld_group_matrix(d$GTs, pr, d$map, size_floor = 2L)
  m_unit <- ld_unit_matrix(d$GTs, d$stage1, d$map, size_floor = 2L, repr = "consensus_dosage")
  expect_false(ncol(m_group) == ncol(m_unit) &&
               identical(sort(colnames(m_group)), sort(as.character(colnames(m_unit)))))
})
