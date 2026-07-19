# Regression tests for plot_pruning_comparison(), in particular the
# zero-matching-clusters edge case (issue: data.table errors out of
# `groups_chr[, .(marker = unlist(members)), by = group_id]` when
# `groups_chr` has zero rows, because `unlist(members)` on an empty
# by-group input returns NULL and data.table can't infer a column type
# from it).

## Small hand-built stage1/result pair with a single unflagged cluster on
## Chr1, so that ld_w_threshold = 0.5 flags nothing on that chromosome --
## the "high" (flagged) direction then has zero matching clusters/groups.
build_stage1_single_cluster <- function() {
  set.seed(1)
  n_ind <- 60
  GTs <- matrix(rbinom(n_ind * 3, 2, 0.3), ncol = 3)
  colnames(GTs) <- paste0("cl1_", 1:3)

  map_snp <- data.table::data.table(
    marker = colnames(GTs), Chr = "Chr1",
    Pos = 1000:1002,
    CL_id = 1L,
    n_loci = 3L,
    is_core = c(TRUE, FALSE, FALSE),
    median_ld = NA_real_,
    ld_w_095 = c(0.05, 0.05, 0.05)
  )

  clusters <- data.table::data.table(
    Chr = "Chr1", CL_id = 1L,
    core_snp = "cl1_1",
    median_ld = 0.9,
    n_snps = 3L,
    members = list(colnames(GTs))
  )

  list(GTs = GTs, stage1 = list(map_snp = map_snp, clusters = clusters))
}

test_that("plot_pruning_comparison() handles zero matching clusters without erroring", {
  d <- build_stage1_single_cluster()

  res <- ld_prune_and_eMLG(
    GTs = d$GTs, stage1 = d$stage1, ld_w_col = "ld_w_095", ld_w_threshold = 0.5,
    score_threshold = 0.80, min_r2 = 0.2, distance_threshold = 100, cores = 1
  )

  ## nothing on Chr1 exceeds ld_w_threshold, so the "high" (flagged) side
  ## has zero Stage-1 clusters and zero combined groups
  expect_equal(nrow(res$groups[startsWith(group_id, "F")]), 0)

  out_dir <- tempfile("plot_pruning_comparison_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  p <- plot_pruning_comparison(
    "Chr1", d$stage1, res, d$stage1$map_snp,
    ld_w_col = "ld_w_095", ld_w_threshold = 0.5,
    direction = "high", out_folder = paste0(out_dir, "/")
  )

  expect_s3_class(p, "patchwork")
  expect_true(file.exists(file.path(out_dir, "Chr1_stage1_vs_combined_high.png")))
})

test_that("plot_pruning_comparison() defaults ld_w_col/ld_w_threshold/min_n_loci_flag from result$params", {
  d <- build_stage1_single_cluster()

  res <- ld_prune_and_eMLG(
    GTs = d$GTs, stage1 = d$stage1, ld_w_col = "ld_w_095", ld_w_threshold = 0.5,
    score_threshold = 0.80, min_r2 = 0.2, distance_threshold = 100, cores = 1
  )
  expect_identical(res$params, list(ld_w_col = "ld_w_095", ld_w_threshold = 0.5, min_n_loci_flag = Inf))

  out_dir <- tempfile("plot_pruning_comparison_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  ## no ld_w_col/ld_w_threshold/min_n_loci_flag passed at all -- should
  ## behave identically to passing them explicitly at their result$params
  ## values (0.5 flags nothing on Chr1, same as the zero-clusters test above)
  p <- plot_pruning_comparison(
    "Chr1", d$stage1, res, d$stage1$map_snp,
    direction = "high", out_folder = paste0(out_dir, "/")
  )

  expect_s3_class(p, "patchwork")
  expect_true(file.exists(file.path(out_dir, "Chr1_stage1_vs_combined_high.png")))
})

test_that("plot_pruning_comparison() warns when ld_w_threshold disagrees with result$params", {
  d <- build_stage1_single_cluster()

  res <- ld_prune_and_eMLG(
    GTs = d$GTs, stage1 = d$stage1, ld_w_col = "ld_w_095", ld_w_threshold = 0.5,
    score_threshold = 0.80, min_r2 = 0.2, distance_threshold = 100, cores = 1
  )

  out_dir <- tempfile("plot_pruning_comparison_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  ## 0.01 flags every marker on Chr1 (all have ld_w_095 = 0.05), unlike the
  ## 0.5 result was actually built with -- Stage 1 and Combined then reflect
  ## different flagging, which should warn
  expect_warning(
    plot_pruning_comparison(
      "Chr1", d$stage1, res, d$stage1$map_snp,
      ld_w_threshold = 0.01, direction = "high", out_folder = paste0(out_dir, "/")
    ),
    "does not match"
  )
})

test_that("plot_pruning_comparison() falls back to historical defaults when result has no params", {
  d <- build_stage1_single_cluster()

  res <- ld_prune_and_eMLG(
    GTs = d$GTs, stage1 = d$stage1, ld_w_col = "ld_w_095", ld_w_threshold = 0.5,
    score_threshold = 0.80, min_r2 = 0.2, distance_threshold = 100, cores = 1
  )
  res$params <- NULL  # simulate a result predating the params field

  out_dir <- tempfile("plot_pruning_comparison_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  ## no explicit args, no result$params -- should use the historical
  ## hardcoded defaults (ld_w_col = "ld_w_095", ld_w_threshold = 0.2) without
  ## erroring or warning
  expect_no_warning(
    p <- plot_pruning_comparison(
      "Chr1", d$stage1, res, d$stage1$map_snp,
      direction = "low", out_folder = paste0(out_dir, "/")
    )
  )
  expect_s3_class(p, "patchwork")
})

test_that("plot_pruning_comparison() still works normally when clusters match (low direction)", {
  d <- build_stage1_single_cluster()

  res <- ld_prune_and_eMLG(
    GTs = d$GTs, stage1 = d$stage1, ld_w_col = "ld_w_095", ld_w_threshold = 0.5,
    score_threshold = 0.80, min_r2 = 0.2, distance_threshold = 100, cores = 1
  )

  out_dir <- tempfile("plot_pruning_comparison_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  p <- plot_pruning_comparison(
    "Chr1", d$stage1, res, d$stage1$map_snp,
    ld_w_col = "ld_w_095", ld_w_threshold = 0.5,
    direction = "low", out_folder = paste0(out_dir, "/")
  )

  expect_s3_class(p, "patchwork")
  expect_true(file.exists(file.path(out_dir, "Chr1_stage1_vs_combined_low.png")))
})
