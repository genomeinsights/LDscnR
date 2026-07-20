# Regression tests for plot_genotype_heatmap().

build_gt_matrix <- function() {
  set.seed(1)
  n_ind <- 20
  n_mk <- 6
  GTs <- matrix(rbinom(n_ind * n_mk, 2, 0.4), nrow = n_ind, ncol = n_mk)
  rownames(GTs) <- paste0("ind", seq_len(n_ind))
  colnames(GTs) <- paste0("mk", seq_len(n_mk))
  GTs
}

test_that("plot_genotype_heatmap() returns a Heatmap object with no annotations", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  GTs <- build_gt_matrix()
  ht <- plot_genotype_heatmap(GTs)

  expect_s4_class(ht, "Heatmap")
})

test_that("plot_genotype_heatmap() accepts col/row annotations and an explicit row_order", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  GTs <- build_gt_matrix()
  col_anno <- setNames(rep(c("A", "B"), each = 3), colnames(GTs))
  row_anno <- setNames(rep(c("pop1", "pop2"), each = 10), rownames(GTs))
  custom_order <- rownames(GTs)[c(11:20, 1:10)]

  ht <- plot_genotype_heatmap(
    GTs,
    row_order = custom_order,
    col_annotation = col_anno, col_annotation_name = "cluster",
    row_annotation = row_anno, row_annotation_name = "population"
  )

  expect_s4_class(ht, "Heatmap")
  expect_identical(rownames(ht@matrix), custom_order)
})

test_that("plot_genotype_heatmap() errors informatively when annotation is missing ids", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  GTs <- build_gt_matrix()
  col_anno <- setNames(rep("A", 3), colnames(GTs)[1:3])  # missing 3 of 6 markers

  expect_error(
    plot_genotype_heatmap(GTs, col_annotation = col_anno),
    "col_annotation is missing"
  )
})

test_that("plot_genotype_heatmap() saves a PNG when out_file is given", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  GTs <- build_gt_matrix()
  out_dir <- tempfile("plot_genotype_heatmap_")
  out_file <- file.path(out_dir, "test_heatmap.png")

  ht <- plot_genotype_heatmap(GTs, out_file = out_file)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  expect_s4_class(ht, "Heatmap")
  expect_true(file.exists(out_file))
})

test_that("plot_genotype_heatmap() polarize=TRUE flips negatively-correlated markers for display only", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  GTs <- build_gt_matrix()
  GTs[, 2] <- 2 - GTs[, 1]  # perfectly anti-correlated with marker 1 (the most-complete reference)

  ht_polarized <- plot_genotype_heatmap(GTs, polarize = TRUE)
  ht_raw <- plot_genotype_heatmap(GTs, polarize = FALSE)

  expect_equal(ht_polarized@matrix[, 2], GTs[, 1])       # flipped to match marker 1
  expect_equal(ht_raw@matrix[, 2], 2 - GTs[, 1])          # unflipped, i.e. original GTs[,2]
})
