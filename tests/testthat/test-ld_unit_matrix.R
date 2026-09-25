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

test_that("repr is restricted to the three Stage-1-unit representations", {
  d <- build_stage1()
  expect_error(
    ld_unit_matrix(d$GTs, d$stage1, d$map, size_floor = 2L, repr = "best_snp"),
    "'arg' should be one of"
  )
})

test_that("the three repr options are unaffected by moving best_snp out to ld_group_matrix()", {
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
