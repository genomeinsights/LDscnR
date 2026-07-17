# Regression tests for ld_prune_and_eMLG(), the combined LD-pruning + eMLG
# generation function built on top of dynamic_cut_eMLG().

## A hand-built ld_complexity_reduction()-shaped stage1 object with 5
## Stage-1 clusters on one chromosome:
##  - CL1, CL2: correlated with each other (shared latent), close together
##    (Pos 1000-1012, gap 8bp) -- should merge under a small
##    distance_threshold.
##  - CL3: correlated with CL1/CL2's latent too, but far away (Pos 5000+,
##    gap ~3988bp) -- should NOT merge with them despite the correlation.
##  - CL1 has only ONE member (cl1_1) with high ld_w; the other two are
##    low -- tests that flagging is cluster-level (any member exceeds
##    threshold), not marker-level.
##  - CL4: independent, low ld_w, 2 loci -- stays unflagged.
##  - CL5: independent, low ld_w, but 6 loci -- used to test
##    min_n_loci_flag pulling a large-but-low-ld_w cluster into merging.
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
  CL5 <- matrix(rbinom(n_ind * 6, 2, 0.3), ncol = 6); colnames(CL5) <- paste0("cl5_", 1:6)

  GTs <- cbind(CL1, CL2, CL3, CL4, CL5)

  map_snp <- data.table::data.table(
    marker = colnames(GTs), Chr = "Chr1",
    Pos = c(1000:1002, 1010:1012, 5000:5002, 2000:2001, 3000:3005),
    CL_id = rep(1:5, c(3, 3, 3, 2, 6)),
    n_loci = rep(c(3, 3, 3, 2, 6), c(3, 3, 3, 2, 6)),
    is_core = unlist(list(c(TRUE, FALSE, FALSE), c(TRUE, FALSE, FALSE), c(TRUE, FALSE, FALSE),
                           c(TRUE, FALSE), c(TRUE, rep(FALSE, 5)))),
    median_ld = NA_real_,
    ld_w_095 = c(0.9, 0.05, 0.05,  0.9, 0.9, 0.9,  0.9, 0.9, 0.9,
                 0.05, 0.05,  0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
  )

  clusters <- data.table::data.table(
    Chr = "Chr1", CL_id = 1:5,
    core_snp = c("cl1_1", "cl2_1", "cl3_1", "cl4_1", "cl5_1"),
    median_ld = c(0.9, 0.9, 0.9, NA, NA),
    n_snps = c(3, 3, 3, 2, 6),
    members = list(colnames(CL1), colnames(CL2), colnames(CL3), colnames(CL4), colnames(CL5))
  )

  list(GTs = GTs, stage1 = list(map_snp = map_snp, clusters = clusters))
}

test_that("flagging is cluster-level, not marker-level, and distance restriction prevents merging a correlated-but-distant cluster", {
  d <- build_stage1()

  res <- ld_prune_and_eMLG(
    GTs = d$GTs, stage1 = d$stage1, ld_w_col = "ld_w_095", ld_w_threshold = 0.5,
    score_threshold = 0.80, min_r2 = 0.2, distance_threshold = 100, cores = 1
  )

  fg <- res$groups[startsWith(group_id, "F")]
  expect_equal(nrow(fg), 2)

  ## CL1+CL2 (close, correlated) merged into one 6-locus group
  merged <- fg[n_loci == 6]
  expect_equal(nrow(merged), 1)
  expect_setequal(merged$members[[1]], c(colnames(d$GTs)[1:3], colnames(d$GTs)[4:6]))

  ## CL3 (correlated with the same latent, but too far away) stayed separate
  distant <- fg[n_loci == 3]
  expect_equal(nrow(distant), 1)
  expect_setequal(distant$members[[1]], colnames(d$GTs)[7:9])

  ## CL1 was flagged even though only ONE of its 3 members exceeded
  ## ld_w_threshold -- confirms cluster-level, not marker-level, flagging
  expect_true(all(colnames(d$GTs)[1:3] %in% merged$members[[1]]))
})

test_that("unflagged clusters pass straight through unchanged", {
  d <- build_stage1()

  res <- ld_prune_and_eMLG(
    GTs = d$GTs, stage1 = d$stage1, ld_w_col = "ld_w_095", ld_w_threshold = 0.5,
    score_threshold = 0.80, min_r2 = 0.2, distance_threshold = 100, cores = 1
  )

  ug <- res$groups[startsWith(group_id, "U")]
  expect_equal(sort(ug$n_loci), c(2, 6))
  expect_equal(ug[n_loci == 2]$representative, "cl4_1")
  expect_equal(ug[n_loci == 2]$members[[1]], colnames(d$GTs)[10:11])

  ## every marker accounted for exactly once across all groups
  expect_setequal(unlist(res$groups$members), colnames(d$GTs))
  expect_equal(length(unlist(res$groups$members)), ncol(d$GTs))
})

test_that("min_n_loci_flag pulls a large low-ld_w cluster into the flagged/merge pathway", {
  d <- build_stage1()

  res <- ld_prune_and_eMLG(
    GTs = d$GTs, stage1 = d$stage1, ld_w_col = "ld_w_095", ld_w_threshold = 0.5,
    score_threshold = 0.80, min_r2 = 0.2, distance_threshold = 100,
    min_n_loci_flag = 5, cores = 1
  )

  ## CL5 (6 loci, low ld_w) is no longer in the unflagged bucket
  ug <- res$groups[startsWith(group_id, "U")]
  expect_false(any(vapply(ug$members, function(m) all(colnames(d$GTs)[12:17] %in% m), logical(1))))

  fg <- res$groups[startsWith(group_id, "F")]
  expect_true(any(vapply(fg$members, function(m) setequal(m, colnames(d$GTs)[12:17]), logical(1))))

  ## total markers covered is unchanged regardless of min_n_loci_flag
  expect_setequal(unlist(res$groups$members), colnames(d$GTs))
})

test_that("compute_unflagged_eMLG=FALSE and min_n_loci_eMLG trim the eMLG matrix but never change the pruned marker set", {
  d <- build_stage1()

  full <- ld_prune_and_eMLG(
    GTs = d$GTs, stage1 = d$stage1, ld_w_col = "ld_w_095", ld_w_threshold = 0.5,
    score_threshold = 0.80, min_r2 = 0.2, distance_threshold = 100, cores = 1
  )
  trimmed <- ld_prune_and_eMLG(
    GTs = d$GTs, stage1 = d$stage1, ld_w_col = "ld_w_095", ld_w_threshold = 0.5,
    score_threshold = 0.80, min_r2 = 0.2, distance_threshold = 100,
    compute_unflagged_eMLG = FALSE, min_n_loci_eMLG = 5, cores = 1
  )

  expect_equal(sort(full$pruned), sort(trimmed$pruned))
  expect_equal(full$groups$group_id, trimmed$groups$group_id)
  expect_equal(full$groups$n_loci, trimmed$groups$n_loci)

  ## trimmed eMLG matrix only keeps groups with n_loci >= 5
  expect_true(all(trimmed$groups[has_eMLG == TRUE]$n_loci >= 5))
  expect_lt(ncol(trimmed$eMLG), ncol(full$eMLG))
})

test_that("cores=1 and cores>1 produce identical output", {
  d <- build_stage1()

  res_c1 <- ld_prune_and_eMLG(
    GTs = d$GTs, stage1 = d$stage1, ld_w_col = "ld_w_095", ld_w_threshold = 0.5,
    score_threshold = 0.80, min_r2 = 0.2, distance_threshold = 100, cores = 1
  )
  res_c2 <- ld_prune_and_eMLG(
    GTs = d$GTs, stage1 = d$stage1, ld_w_col = "ld_w_095", ld_w_threshold = 0.5,
    score_threshold = 0.80, min_r2 = 0.2, distance_threshold = 100, cores = 2
  )

  expect_identical(res_c1$eMLG, res_c2$eMLG)
  expect_identical(res_c1$groups, res_c2$groups)
  expect_identical(res_c1$pruned, res_c2$pruned)
})
