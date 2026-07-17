# Regression tests for dynamic_cut_eMLG(), the quality-gated average-linkage
# tree walk that consolidates cluster-level eMLGs. See the file header in
# R/dynamic_cut_eMLG.R for the linkage comparison and the two historical
# bugs these tests target directly.

test_that("correlated clusters merge into one group, an independent one stays separate", {
  set.seed(1)
  n <- 80
  latent <- rbinom(n, 2, 0.3)
  A <- pmin(pmax(latent + rbinom(n, 1, 0.05) - rbinom(n, 1, 0.05), 0), 2)
  B <- pmin(pmax(latent + rbinom(n, 1, 0.05) - rbinom(n, 1, 0.05), 0), 2)
  C <- rbinom(n, 2, 0.3)  # independent

  eMLG <- cbind(A = A, B = B, C = C)
  n_loci <- c(A = 1, B = 1, C = 1)

  res <- dynamic_cut_eMLG(eMLG, n_loci, threshold = 0.80, min_r2 = 0.2)

  sizes <- sort(lengths(lapply(res, `[[`, "members")))
  expect_equal(sizes, c(1, 2))

  merged <- res[[which.max(lengths(lapply(res, `[[`, "members")))]]
  expect_setequal(merged$members, c("A", "B"))
  expect_gte(merged$score, 0.80)
})

test_that("two near-monomorphic, uncorrelated clusters do NOT merge just because their weighted average happens to round cleanly (bug #2 regression)", {
  ## deterministic construction: a single deviation from 0 (A) and a single
  ## deviation from 2 (B) at different positions gives pairwise r2 ~ 2.5e-05
  ## (matches the historically observed 1.8e-05) while the resulting
  ## weighted-average consensus still scores a PERFECT 1.0 on score_eMLG --
  ## verified directly before writing this test. score_eMLG alone would
  ## wrongly accept this merge; the min_r2 gate must reject it.
  n <- 200
  A <- rep(0, n); A[1] <- 1
  B <- rep(2, n); B[2] <- 1

  pair_r2 <- suppressWarnings(cor(A, B))^2
  expect_lt(pair_r2, 0.001)

  candidate_score <- score_eMLG(weighted_row_mean(polarize_genotypes(cbind(A, B)), c(1, 1)))
  expect_equal(candidate_score, 1)  # confirms this WOULD pass a score_eMLG-only gate

  eMLG <- cbind(A = A, B = B)
  n_loci <- c(A = 1, B = 1)

  res <- dynamic_cut_eMLG(eMLG, n_loci, threshold = 0.80, min_r2 = 0.2)
  expect_length(res, 2)  # stayed separate

  ## with the min_r2 gate disabled, they DO merge -- demonstrates this test
  ## actually exercises the gate, not some other rejection reason
  res_nogate <- dynamic_cut_eMLG(eMLG, n_loci, threshold = 0.80, min_r2 = 0)
  expect_length(res_nogate, 1)
})

test_that("no clusters are silently lost across a mix of accepted and rejected merges (bug #1 regression)", {
  ## a frozen sibling's active partner must be finalized immediately, not
  ## left dangling -- missing this dropped ~90% of markers historically.
  ## Mixing a 3-cluster block, a 2-cluster block, and independent
  ## singletons in one call forces the dendrogram walk through both
  ## accepted and rejected merges, which is what originally caught this.
  set.seed(7)
  n <- 80
  mk_block <- function(latent, k, noise = 0.05) {
    sapply(seq_len(k), function(i) {
      pmin(pmax(latent + rbinom(n, 1, noise) - rbinom(n, 1, noise), 0), 2)
    })
  }
  latent1 <- rbinom(n, 2, 0.3)
  latent2 <- rbinom(n, 2, 0.5)
  latent3 <- rbinom(n, 2, 0.2)

  blockA <- mk_block(latent1, 3)
  blockB <- mk_block(latent2, 2)
  indepC <- mk_block(latent3, 1)
  indepD <- matrix(rbinom(n * 2, 2, 0.3), ncol = 2)

  ## name the block columns BEFORE cbinding -- blockA/blockB come from
  ## sapply() and have no colnames of their own otherwise
  colnames(blockA) <- paste0("a", seq_len(ncol(blockA)))
  colnames(blockB) <- paste0("b", seq_len(ncol(blockB)))
  colnames(indepC) <- "c1"
  colnames(indepD) <- paste0("d", seq_len(ncol(indepD)))

  eMLG <- cbind(blockA, blockB, indepC, indepD)
  n_loci <- setNames(rep(1, ncol(eMLG)), colnames(eMLG))

  res <- dynamic_cut_eMLG(eMLG, n_loci, threshold = 0.80, min_r2 = 0.2)

  all_members <- unlist(lapply(res, `[[`, "members"))
  expect_length(all_members, ncol(eMLG))
  expect_setequal(all_members, colnames(eMLG))

  ## the two real blocks should each end up as one group
  expect_true(any(vapply(res, function(g) setequal(g$members, colnames(blockA)), logical(1))))
  expect_true(any(vapply(res, function(g) setequal(g$members, colnames(blockB)), logical(1))))
})

test_that("merging weights by n_loci, not equally, so a single-marker cluster can't outweigh a large one", {
  set.seed(1)
  n <- 60
  latent <- rbinom(n, 2, 0.4)
  A <- pmin(pmax(latent + rbinom(n, 1, 0.03) - rbinom(n, 1, 0.03), 0), 2)  # weight 100
  B <- pmin(pmax(latent + rbinom(n, 1, 0.15) - rbinom(n, 1, 0.15), 0), 2)  # weight 1, noisier

  eMLG <- cbind(A = A, B = B)
  n_loci <- c(A = 100, B = 1)

  res <- dynamic_cut_eMLG(eMLG, n_loci, threshold = 0.80, min_r2 = 0.2)
  expect_length(res, 1)

  merged <- res[[1]]$emlg
  expected_mean <- (100 * mean(A) + 1 * mean(B)) / 101
  expect_equal(mean(merged), expected_mean, tolerance = 1e-8)
  expect_gt(cor(merged, A), 0.99)
})

test_that("every input cluster ID must have an n_loci entry", {
  eMLG <- cbind(A = rbinom(20, 2, 0.3), B = rbinom(20, 2, 0.3))
  expect_error(dynamic_cut_eMLG(eMLG, c(A = 1), threshold = 0.80), "n_loci")
})

test_that("eMLG column names must be unique", {
  eMLG <- cbind(rbinom(20, 2, 0.3), rbinom(20, 2, 0.3))
  colnames(eMLG) <- c("A", "A")
  expect_error(dynamic_cut_eMLG(eMLG, c(A = 1), threshold = 0.80), "unique")
})
