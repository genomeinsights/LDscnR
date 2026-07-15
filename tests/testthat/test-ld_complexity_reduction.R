# Regression tests for ld_complexity_reduction(), exercised end-to-end on
# the bundled sim_ex dataset: create_gds_from_geno() -> compute_LD_decay()
# -> ld_complexity_reduction().

build_decay <- function(keep_el = TRUE) {
  skip_if_not_installed("SNPRelate")
  data(sim_ex, package = "LDscnR")

  gds_path <- tempfile(fileext = ".gds")
  gds <- create_gds_from_geno(sim_ex$GTs, sim_ex$map, gds_path)
  on.exit({
    SNPRelate::snpgdsClose(gds)
    unlink(gds_path)
  })

  compute_LD_decay(
    gds,
    n_win_decay    = 5,
    max_SNPs_decay = 2000,
    slide          = 1000,
    keep_el        = keep_el,
    cores          = 1
  )
}

test_that("ld_complexity_reduction returns a well-formed, self-consistent result", {
  data(sim_ex, package = "LDscnR")
  map <- data.table::as.data.table(sim_ex$map)[, .(Chr, Pos, marker)]
  ld  <- build_decay()

  res <- ld_complexity_reduction(sim_ex$GTs, map, ld, rho = 0.5, cores = 1)

  expect_named(res, c("map_snp", "clusters", "pruned"))

  # every marker retained exactly once, every marker assigned a cluster
  expect_equal(nrow(res$map_snp), nrow(map))
  expect_true(all(!is.na(res$map_snp$CL_id)))

  # exactly one core marker per cluster
  n_core_by_cluster <- res$map_snp[, sum(is_core), by = CL_id]$V1
  expect_true(all(n_core_by_cluster == 1))

  # cluster-level table is internally consistent with the per-marker table
  expect_equal(res$clusters$n_snps, lengths(res$clusters$members))
  expect_equal(sort(res$pruned), sort(res$clusters$core_snp))
  expect_true(all(res$pruned %in% res$map_snp[is_core == TRUE, marker]))

  # singletons have NA median_ld; non-singletons don't
  singleton_markers <- res$map_snp[n_loci == 1, marker]
  if (length(singleton_markers)) {
    expect_true(all(is.na(res$map_snp[marker %in% singleton_markers, median_ld])))
  }
  non_singleton <- res$map_snp[n_loci > 1]
  if (nrow(non_singleton)) {
    expect_true(all(!is.na(non_singleton$median_ld)))
  }
})

test_that("clusters are true complete linkage: every member pair clears r2_th", {
  data(sim_ex, package = "LDscnR")
  map <- data.table::as.data.table(sim_ex$map)[, .(Chr, Pos, marker)]
  ld  <- build_decay()

  # create_gds_from_geno() only relies on positional alignment between GTs
  # columns and map rows, not name matching -- sim_ex$GTs's own colnames
  # (e.g. "1:5570") don't carry the "Chr" prefix map$marker does (e.g.
  # "Chr1:5570"), so name-based indexing below needs them aligned first.
  GTs <- sim_ex$GTs
  colnames(GTs) <- map$marker

  rho <- 0.5
  res <- ld_complexity_reduction(GTs, map, ld, rho = rho, cores = 1)

  big_clusters <- res$clusters[n_snps >= 3][order(-n_snps)][seq_len(min(5, .N))]

  for (i in seq_len(nrow(big_clusters))) {
    ch      <- big_clusters$Chr[i]
    members <- big_clusters$members[[i]]

    ds    <- ld$decay_sum[Chr == ch]
    r2_th <- ld_from_rho(b = ds$b, c = ds$c, rho = rho)

    R2 <- suppressWarnings(cor(GTs[, members, drop = FALSE], use = "pairwise.complete.obs")^2)
    diag(R2) <- NA

    # every pairwise r2 within a final cluster must clear r2_th -- the
    # complete-linkage guarantee that specifically avoids chaining
    expect_true(all(R2 >= r2_th, na.rm = TRUE))
  }
})

test_that("ld_w_col/ld_w_threshold lets low-ld_w markers skip clustering", {
  data(sim_ex, package = "LDscnR")
  map <- data.table::as.data.table(sim_ex$map)[, .(Chr, Pos, marker)]
  ld  <- build_decay()

  map[, ld_w := as.numeric(compute_ld_w(ld, rho = 0.95))]

  res_full     <- ld_complexity_reduction(sim_ex$GTs, map, ld, rho = 0.5, cores = 1)
  res_filtered <- ld_complexity_reduction(
    sim_ex$GTs, map, ld, rho = 0.5, cores = 1,
    ld_w_col = "ld_w", ld_w_threshold = 0.2
  )

  low_ld_w_markers <- map[ld_w < 0.2, marker]

  # every low-ld_w marker is its own singleton cluster when pre-filtered
  expect_true(all(
    res_filtered$map_snp[marker %in% low_ld_w_markers, n_loci] == 1
  ))
  expect_true(all(low_ld_w_markers %in% res_filtered$pruned))

  # same total marker count either way, nothing dropped
  expect_equal(nrow(res_full$map_snp), nrow(res_filtered$map_snp))
})
