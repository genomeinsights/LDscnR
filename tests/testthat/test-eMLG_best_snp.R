# Tests for eMLG_best_snp(): best consensus-correlated SNP per block + fill.

# Build a small eMLG result once from the bundled sim_ex data.
build_eMLG_result <- function() {
  skip_if_not_installed("SNPRelate")
  data(sim_ex, package = "LDscnR")
  map <- sim_ex$map
  GTs <- sim_ex$GTs
  colnames(GTs) <- map$marker

  gds_path <- tempfile(fileext = ".gds")
  gds <- create_gds_from_geno(GTs, map, gds_path)
  on.exit({ SNPRelate::snpgdsClose(gds); unlink(gds_path) })

  ld_decay <- compute_LD_decay(gds, n_win_decay = 5, slide = 200,
                               keep_el = TRUE, cores = 1)
  ld_w <- compute_ld_w(ld_decay, rho = 0.95, cores = 1)
  ids  <- unlist(lapply(ld_decay$by_chr, function(o) o$snp_ids), use.names = FALSE)
  map[, ld_w_095 := ld_w[match(marker, ids)]]

  stage1 <- ld_complexity_reduction(map = map, LD_decay = ld_decay,
                                    rho = 0.5, cores = 1)
  res <- ld_prune_and_eMLG(GTs = GTs, stage1 = stage1, ld_w_col = "ld_w_095",
                           ld_w_threshold = 0, LD_decay = ld_decay, rho = 0.95,
                           score_threshold = 0.80, min_r2 = 0.2,
                           compute_unflagged_eMLG = TRUE, min_n_loci_eMLG = 3,
                           cores = 1)
  list(res = res, GTs = GTs)
}

test_that("eMLG_best_snp returns one stats row per eMLG block", {
  be <- build_eMLG_result()
  skip_if(is.null(be$res$eMLG) || ncol(be$res$eMLG) == 0, "no eMLG blocks in sim_ex")

  out <- eMLG_best_snp(be$res, be$GTs)
  expect_equal(nrow(out$stats), ncol(be$res$eMLG))
  expect_setequal(out$stats$group_id, colnames(be$res$eMLG))
  expect_true(all(out$stats$best_marker %in% colnames(be$GTs)))
})

test_that("best member is at least as correlated as the representative", {
  be <- build_eMLG_result()
  skip_if(is.null(be$res$eMLG) || ncol(be$res$eMLG) == 0, "no eMLG blocks in sim_ex")

  out <- eMLG_best_snp(be$res, be$GTs)
  ok <- is.finite(out$stats$best_abs_r) & is.finite(out$stats$rep_abs_r)
  # allow tiny numerical slack
  expect_true(all(out$stats$best_abs_r[ok] >= out$stats$rep_abs_r[ok] - 1e-8))
})

test_that("fill preserves observed calls and reduces missingness", {
  be <- build_eMLG_result()
  skip_if(is.null(be$res$eMLG) || ncol(be$res$eMLG) == 0, "no eMLG blocks in sim_ex")

  out <- eMLG_best_snp(be$res, be$GTs, fill = TRUE, round_fill = TRUE)

  # geno aligns with the eMLG matrix (same shape, same column names)
  expect_equal(dim(out$geno), dim(be$res$eMLG))
  expect_equal(colnames(out$geno), colnames(be$res$eMLG))
  expect_true(all(out$geno %in% c(0L, 1L, 2L, NA_integer_)))

  # observed calls of the chosen SNP are preserved exactly
  G <- if (!is.null(rownames(be$res$eMLG))) be$GTs[rownames(be$res$eMLG), ] else be$GTs
  snp_raw <- G[, out$stats$best_marker]
  obs <- !is.na(snp_raw)
  expect_true(all(out$geno[obs] == snp_raw[obs]))

  # filling never increases missingness
  expect_lte(sum(is.na(out$geno)), sum(is.na(snp_raw)))

  # accounting adds up
  expect_equal(out$stats$n_obs + out$stats$n_filled + out$stats$n_resid_na,
               rep(nrow(be$res$eMLG), nrow(out$stats)))
})

test_that("fill = FALSE returns no geno matrix", {
  be <- build_eMLG_result()
  skip_if(is.null(be$res$eMLG) || ncol(be$res$eMLG) == 0, "no eMLG blocks in sim_ex")

  out <- eMLG_best_snp(be$res, be$GTs, fill = FALSE)
  expect_null(out$geno)
  expect_false("n_obs" %in% names(out$stats))
})
