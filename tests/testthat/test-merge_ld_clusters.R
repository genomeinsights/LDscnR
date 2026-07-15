# Regression tests for merge_ld_clusters(), the cluster-level second stage
# on top of ld_complexity_reduction().

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

test_that("merge_ld_clusters reunites an artificially fragmented block and leaves an independent cluster alone", {
  set.seed(1)
  n_ind <- 100

  latent <- rbinom(n_ind, 2, 0.3)
  block <- sapply(1:6, function(i) {
    pmin(pmax(latent + rbinom(n_ind, 1, 0.05) - rbinom(n_ind, 1, 0.05), 0), 2)
  })
  colnames(block) <- paste0("m", 1:6)

  indep_latent <- rbinom(n_ind, 2, 0.3)
  indep <- sapply(1:3, function(i) {
    pmin(pmax(indep_latent + rbinom(n_ind, 1, 0.05) - rbinom(n_ind, 1, 0.05), 0), 2)
  })
  colnames(indep) <- paste0("i", 1:3)

  GTs <- cbind(block, indep)

  ## simulate ld_complexity_reduction()'s output: the block artificially
  ## split into two clusters (m1-m3, m4-m6), independent cluster kept apart
  ld_result <- list(
    map_snp = data.table::data.table(
      marker = colnames(GTs), Chr = "Chr1", Pos = 1:9,
      CL_id = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
      n_loci = 3L,
      is_core = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE),
      median_ld = NA_real_
    ),
    clusters = data.table::data.table(
      Chr = "Chr1", CL_id = c(1, 2, 3),
      core_snp = c("m1", "m4", "i1"),
      median_ld = c(0.9, 0.9, 0.9),
      n_snps = c(3L, 3L, 3L),
      members = list(c("m1", "m2", "m3"), c("m4", "m5", "m6"), c("i1", "i2", "i3"))
    )
  )

  LD_decay <- list(decay_sum = data.table::data.table(Chr = "Chr1", b = 0, c = 1))

  res <- merge_ld_clusters(GTs, ld_result, LD_decay, rho = 0.5, cores = 1)

  expect_named(res, c("map_snp", "clusters", "pruned"))
  expect_equal(nrow(res$clusters), 2)
  expect_equal(nrow(res$map_snp), nrow(ld_result$map_snp))

  expect_equal(uniqueN(res$map_snp[marker %in% paste0("m", 1:6), CL_id]), 1)
  expect_false(res$map_snp[marker == "i1", CL_id] == res$map_snp[marker == "m1", CL_id])

  ## exactly one core per cluster, n_snps matches actual member-list length
  expect_true(all(res$map_snp[, sum(is_core), by = CL_id]$V1 == 1))
  expect_equal(res$clusters$n_snps, lengths(res$clusters$members))
})

test_that("density is correctly aligned per cluster, not scrambled across the merged/passthrough branches", {
  set.seed(1)
  n_ind <- 100

  latent <- rbinom(n_ind, 2, 0.3)
  block <- sapply(1:6, function(i) {
    pmin(pmax(latent + rbinom(n_ind, 1, 0.05) - rbinom(n_ind, 1, 0.05), 0), 2)
  })
  colnames(block) <- paste0("m", 1:6)

  indep_latent <- rbinom(n_ind, 2, 0.3)
  indep <- sapply(1:3, function(i) {
    pmin(pmax(indep_latent + rbinom(n_ind, 1, 0.05) - rbinom(n_ind, 1, 0.05), 0), 2)
  })
  colnames(indep) <- paste0("i", 1:3)

  GTs <- cbind(block, indep)

  ld_result <- list(
    map_snp = data.table::data.table(
      marker = colnames(GTs), Chr = "Chr1", Pos = 1:9,
      CL_id = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
      n_loci = 3L,
      is_core = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE),
      median_ld = NA_real_
    ),
    clusters = data.table::data.table(
      Chr = "Chr1", CL_id = c(1, 2, 3),
      core_snp = c("m1", "m4", "i1"),
      median_ld = c(0.9, 0.9, 0.9),
      n_snps = c(3L, 3L, 3L),
      members = list(c("m1", "m2", "m3"), c("m4", "m5", "m6"), c("i1", "i2", "i3"))
    )
  )
  LD_decay <- list(decay_sum = data.table::data.table(Chr = "Chr1", b = 0, c = 1))

  res <- merge_ld_clusters(GTs, ld_result, LD_decay, rho = 0.5, cores = 1)
  expect_true("density" %in% names(res$clusters))

  ## independently recompute density for each cluster from its OWN members
  ## and check it landed on the right row -- this is the exact scenario
  ## that caught a real rbindlist() column-misalignment bug during
  ## development (merged and passthrough branches build density in
  ## different column positions; rbindlist()'s use.names="check" default
  ## only warns and still binds positionally)
  for (i in seq_len(nrow(res$clusters))) {
    mk <- res$clusters$members[[i]]
    r2 <- suppressWarnings(cor(GTs[, mk, drop = FALSE], use = "pairwise.complete.obs")^2)
    expected <- mean(r2[upper.tri(r2)] >= 0.5)
    expect_equal(res$clusters$density[i], expected,
                 info = paste("cluster", res$clusters$CL_id[i]))
  }
})

test_that("merge_ld_clusters on real data: self-consistent and never increases cluster count", {
  data(sim_ex, package = "LDscnR")
  map <- data.table::as.data.table(sim_ex$map)[, .(Chr, Pos, marker)]
  ld  <- build_decay()

  GTs <- sim_ex$GTs
  colnames(GTs) <- map$marker  ## see test-ld_complexity_reduction.R for why

  res_stage1 <- ld_complexity_reduction(map = map, LD_decay = ld, rho = 0.5, cores = 1)
  res_merged <- merge_ld_clusters(GTs, res_stage1, ld, rho = 0.5, cores = 1)

  expect_named(res_merged, c("map_snp", "clusters", "pruned"))

  ## merging can only reduce (or hold) the number of clusters, never increase it
  expect_lte(nrow(res_merged$clusters), nrow(res_stage1$clusters))

  ## every marker retained exactly once, still fully assigned
  expect_equal(nrow(res_merged$map_snp), nrow(res_stage1$map_snp))
  expect_true(all(!is.na(res_merged$map_snp$CL_id)))

  ## still self-consistent: one core per cluster, sizes match, pruned == core_snp
  expect_true(all(res_merged$map_snp[, sum(is_core), by = CL_id]$V1 == 1))
  expect_equal(res_merged$clusters$n_snps, lengths(res_merged$clusters$members))
  expect_equal(sort(res_merged$pruned), sort(res_merged$clusters$core_snp))

  ## density: present for every cluster, in [0,1] wherever defined (NA only
  ## for true singletons) -- real multi-chromosome data is what actually
  ## caught the CL_id integer-vs-character indexing bug during development,
  ## so this checks on sim_ex specifically, not just the single-chromosome
  ## synthetic case above
  expect_true("density" %in% names(res_merged$clusters))
  singleton <- res_merged$clusters$n_snps == 1
  expect_true(all(is.na(res_merged$clusters$density[singleton])))
  non_singleton_density <- res_merged$clusters$density[!singleton]
  expect_true(all(!is.na(non_singleton_density)))
  expect_true(all(non_singleton_density >= 0 & non_singleton_density <= 1))

  ## every original marker is still accounted for somewhere in the merged clusters
  expect_setequal(
    unlist(res_merged$clusters$members),
    unlist(res_stage1$clusters$members)
  )
})

test_that("merge_ld_clusters only merges clusters on the same chromosome", {
  data(sim_ex, package = "LDscnR")
  map <- data.table::as.data.table(sim_ex$map)[, .(Chr, Pos, marker)]
  ld  <- build_decay()

  GTs <- sim_ex$GTs
  colnames(GTs) <- map$marker

  res_stage1 <- ld_complexity_reduction(map = map, LD_decay = ld, rho = 0.5, cores = 1)
  res_merged <- merge_ld_clusters(GTs, res_stage1, ld, rho = 0.5, cores = 1)

  ## a merged cluster's members must all still belong to the SAME chromosome
  ## in map_snp -- merging across chromosomes would be physically meaningless
  chr_by_marker <- res_merged$map_snp[, .(marker, Chr)]
  per_cluster_chr <- merge(
    data.table::data.table(marker = unlist(res_merged$clusters$members),
                            CL_id = rep(res_merged$clusters$CL_id, lengths(res_merged$clusters$members))),
    chr_by_marker, by = "marker"
  )
  n_chr_by_cluster <- per_cluster_chr[, uniqueN(Chr), by = CL_id]$V1
  expect_true(all(n_chr_by_cluster == 1))
})
