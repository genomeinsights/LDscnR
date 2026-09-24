## Fixture shared in spirit with test-ld_prune_and_eMLG.R's build_stage1():
## CL1+CL2 are close together and correlated -> merge into one "F..." group
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

test_that("best_snp errors without prune_result", {
  d <- build_stage1()
  expect_error(
    ld_unit_matrix(d$GTs, d$stage1, d$map, size_floor = 2L, repr = "best_snp"),
    "prune_result"
  )
})

test_that("best_snp errors if fill = FALSE is passed through best_snp_args", {
  d <- build_stage1()
  pr <- build_prune_result(d)
  expect_error(
    ld_unit_matrix(d$GTs, d$stage1, d$map, size_floor = 2L, repr = "best_snp",
                   prune_result = pr, best_snp_args = list(fill = FALSE)),
    "fill = TRUE"
  )
})

test_that("best_snp reports prune_result's own groups, not stage1's units, with exact column alignment", {
  d <- build_stage1()
  pr <- build_prune_result(d)
  ## sanity on the fixture: CL1+CL2 merged into one 6-locus group (spans two
  ## Stage-1 clusters), CL4 passed through unflagged (1:1 with Stage-1 CL4)
  merged_id <- pr$groups[startsWith(group_id, "F") & n_loci == 6, group_id]
  expect_length(merged_id, 1)
  unflagged_id <- pr$groups[startsWith(group_id, "U") & n_loci == 2, group_id]
  expect_length(unflagged_id, 1)

  m <- ld_unit_matrix(d$GTs, d$stage1, d$map, size_floor = 2L, repr = "best_snp", prune_result = pr)
  u <- attr(m, "units")

  ## column order matches attr(,"units") row order exactly
  expect_identical(colnames(m), u$unit_id)
  expect_identical(nrow(m), nrow(d$GTs))

  ## unit_id is prune_result's group_id (a character id, e.g. "F1"/"U3"),
  ## not stage1's own integer unit_id
  expect_type(u$unit_id, "character")
  expect_true(merged_id %in% u$unit_id)
  expect_true(unflagged_id %in% u$unit_id)

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

test_that("best_snp's size_floor filters prune_result's groups by n_loci, independent of stage1's own filtering", {
  d <- build_stage1()
  pr <- build_prune_result(d)

  ## size_floor = 5 exceeds stage1's largest RAW cluster (3), which would
  ## make .ld_outlier_units() error -- but the merged 6-locus "F" group
  ## still clears it, so best_snp must not fail here.
  m <- ld_unit_matrix(d$GTs, d$stage1, d$map, size_floor = 5L, repr = "best_snp", prune_result = pr)
  u <- attr(m, "units")
  expect_true(all(u$n_markers >= 5L))
  expect_true(all(startsWith(u$unit_id, "F")))  # the 2-locus "U" group is filtered out
})

test_that("best_snp errors informatively when no prune_result group clears size_floor", {
  d <- build_stage1()
  pr <- build_prune_result(d)
  expect_error(
    ld_unit_matrix(d$GTs, d$stage1, d$map, size_floor = 1000L, repr = "best_snp", prune_result = pr),
    "size_floor"
  )
})

test_that("the other three repr options are unaffected by the best_snp changes", {
  d <- build_stage1()
  m_cons <- ld_unit_matrix(d$GTs, d$stage1, d$map, size_floor = 2L, repr = "consensus_dosage")
  expect_equal(nrow(m_cons), nrow(d$GTs))
  expect_true(all(c("unit_id", "Chr", "from", "to", "n_markers") %in% names(attr(m_cons, "units"))))

  v_rep <- ld_unit_matrix(d$GTs, d$stage1, d$map, size_floor = 2L, repr = "representative")
  expect_true(is.character(v_rep))
  expect_true(all(v_rep %in% colnames(d$GTs)))
})

test_that("consensus_dosage/eMLG error informatively when a tested unit's marker is missing from GTs", {
  d <- build_stage1()
  bad_GTs <- d$GTs[, colnames(d$GTs) != "cl1_1"]  # cl1_1 is CL1's core member, CL1 clears size_floor=2
  expect_error(
    ld_unit_matrix(bad_GTs, d$stage1, d$map, size_floor = 2L, repr = "consensus_dosage"),
    "missing from `colnames\\(GTs\\)`"
  )
  expect_error(
    ld_unit_matrix(bad_GTs, d$stage1, d$map, size_floor = 2L, repr = "eMLG"),
    "missing from `colnames\\(GTs\\)`"
  )
  ## "representative" never touches GTs, so it is unaffected
  expect_no_error(
    ld_unit_matrix(bad_GTs, d$stage1, d$map, size_floor = 2L, repr = "representative")
  )
})

test_that("duplicate map markers are rejected before any expensive computation", {
  d <- build_stage1()
  bad_map <- data.table::copy(d$map)
  bad_map$marker[2] <- bad_map$marker[1]
  expect_error(
    ld_unit_matrix(d$GTs, d$stage1, bad_map, size_floor = 2L, repr = "consensus_dosage"),
    "duplicate"
  )
})
