# Tests for the p-value-first API: ld_null_from_p() -> ld_gate() ->
# ld_region_scan() -> ld_region_c2(), and the ld_scan() wrapper.
#
# A small synthetic panel with real LD blocks and one planted signal region.
# The pipeline must recover the planted region, call nothing else significant,
# and stay silent when the "observed" data is itself a surrogate.

make_case <- function(seed = 42, n_ind = 60, n_blk = 40, blk_size = 10) {
  set.seed(seed)
  n_mk <- n_blk * blk_size
  map <- data.table::data.table(
    Chr = rep(c("Chr1", "Chr2"), each = n_mk / 2),
    Pos = rep(seq_len(n_mk / 2) * 1e4, 2))
  map[, marker := paste0(Chr, ":", Pos)]

  # genotypes in LD blocks, so ld_edges() finds genuine edges
  blk <- rep(seq_len(n_blk), each = blk_size)
  GTs <- matrix(0, n_ind, n_mk, dimnames = list(NULL, map$marker))
  for (b in unique(blk)) {
    h <- stats::rbinom(n_ind, 2, 0.3)
    for (j in which(blk == b))
      GTs[, j] <- pmin(2, pmax(0, h + stats::rbinom(n_ind, 1, 0.08) -
                                    stats::rbinom(n_ind, 1, 0.08)))
  }
  keep <- apply(GTs, 2, stats::sd) > 0
  GTs <- GTs[, keep]; map <- map[keep]
  n_mk <- ncol(GTs)

  ld_ws <- matrix(stats::runif(n_mk * 3, 0.1, 0.6), n_mk, 3,
                  dimnames = list(colnames(GTs), c("rho_0.9", "rho_0.95", "rho_0.99")))
  decay_sum <- data.table::data.table(Chr = c("Chr1", "Chr2"), a = c(0.4, 0.4),
                                      b = c(1e-5, 1e-5), c = c(1, 1))
  sig_idx <- 51:60                              # one whole LD block on Chr1
  mk_p <- function(signal) {
    p <- stats::runif(n_mk)
    if (signal) p[sig_idx] <- stats::runif(length(sig_idx), 0, 1e-8)
    p
  }
  list(map = map, GTs = GTs, ld_ws = ld_ws, decay_sum = decay_sum, sig_idx = sig_idx,
       p_obs = mk_p(TRUE), p_perm = replicate(50, mk_p(FALSE), simplify = FALSE))
}

edges_for <- function(cs, null)
  ld_edges(null$universe, cs$GTs, cs$map[, .(marker, Chr, Pos)],
           cs$decay_sum, rho_ld = 0.6, dcap = 5e5)

test_that("ld_null_from_p accepts matrix, list and function identically", {
  cs <- make_case()
  n_l <- ld_null_from_p(cs$p_obs, cs$p_perm, cs$ld_ws, verbose = FALSE)
  n_m <- ld_null_from_p(cs$p_obs, do.call(cbind, cs$p_perm), cs$ld_ws, verbose = FALSE)
  n_f <- ld_null_from_p(cs$p_obs, function(b) cs$p_perm[[b]], cs$ld_ws, B = 50, verbose = FALSE)

  expect_s3_class(n_l, "ld_null")
  expect_identical(n_l$C_surr, n_m$C_surr)
  expect_identical(n_l$C_surr, n_f$C_surr)
  expect_identical(n_l$B, 50L)
  # surrogate C-vectors are stored sparsely (C > 0 only)
  expect_true(all(vapply(n_l$C_surr, function(s) all(s > 0), logical(1))))
})

test_that("ld_null_from_p validates its inputs", {
  cs <- make_case()
  expect_error(ld_null_from_p(cs$p_obs[-1], cs$p_perm, cs$ld_ws, verbose = FALSE), "aligned")
  expect_error(ld_null_from_p(cs$p_obs, function(b) cs$p_perm[[b]], cs$ld_ws, verbose = FALSE), "`B` is required")
  expect_error(ld_null_from_p(cs$p_obs, "nonsense", cs$ld_ws, verbose = FALSE), "must be a matrix")
  expect_warning(ld_null_from_p(cs$p_obs, cs$p_perm[1:3], cs$ld_ws, verbose = FALSE), "Only 3 surrogate")
})

test_that("ld_region_scan recovers a planted region and nothing else", {
  cs <- make_case()
  null <- ld_null_from_p(cs$p_obs, cs$p_perm, cs$ld_ws, verbose = FALSE)
  s <- ld_region_scan(null, edges_for(cs, null), tau = 0.05, l_min = 3L)

  expect_s3_class(s, "ld_region_scan")
  expect_equal(sum(s$regions$sig), 1L)
  hit <- s$regions[sig == TRUE]
  expect_identical(hit$chr, "Chr1")
  # its span must cover the planted block
  expect_lte(hit$lo, min(cs$map$Pos[cs$sig_idx]))
  expect_gte(hit$hi, max(cs$map$Pos[cs$sig_idx]))
  # no surrogate reaches it, so p sits at the floor
  expect_equal(hit$p_R, 1 / (1 + null$B))
  expect_equal(s$p_floor, 1 / (1 + null$B))
})

test_that("a surrogate used as the observed data yields no significant region", {
  cs <- make_case()
  null0 <- ld_null_from_p(cs$p_perm[[1]], cs$p_perm[-1], cs$ld_ws, verbose = FALSE)
  s0 <- ld_region_scan(null0, edges_for(cs, null0), tau = 0.05, l_min = 3L)
  expect_equal(sum(s0$regions$sig), 0L)
})

test_that("ld_gate reports medians per surrogate and flags a bad basis", {
  cs <- make_case()
  null <- ld_null_from_p(cs$p_obs, cs$p_perm, cs$ld_ws, verbose = FALSE)
  g <- ld_gate(null, edges_for(cs, null), tau = 0.05, l_min = 3L)

  expect_s3_class(g, "ld_gate")
  expect_equal(nrow(g), 1L)
  expect_true(g$pass)
  expect_equal(g$obs_regions, 1L)
  expect_true(g$med_regions <= g$obs_regions)
  expect_true(all(c("q3_regions", "max_regions", "frac_surr_any_region") %in% names(g)))

  # a null whose "surrogates" are the observed data must fail the gate
  bad <- null
  bad$C_surr <- rep(list(null$C_obs[null$C_obs > 0]), 10L)
  bad$B <- 10L
  expect_warning(gb <- ld_gate(bad, edges_for(cs, null), tau = 0.05, l_min = 3L), "Gate FAILED")
  expect_false(gb$pass)
})

test_that("ld_gate accepts a named list of nulls, one row each", {
  cs <- make_case()
  null <- ld_null_from_p(cs$p_obs, cs$p_perm, cs$ld_ws, verbose = FALSE)
  g <- ld_gate(list(a = null, b = null), edges_for(cs, null), tau = 0.05, l_min = 3L)
  expect_equal(nrow(g), 2L)
  expect_identical(g$basis, c("a", "b"))
})

test_that("ld_region_c2 returns scores in (0, 1] normalised by the usable grid", {
  cs <- make_case()
  null <- ld_null_from_p(cs$p_obs, cs$p_perm, cs$ld_ws, verbose = FALSE)
  c2 <- ld_region_c2(null, edges_for(cs, null), tau_grid = seq(0.05, 0.3, by = 0.05),
                     lmin_grid = c(2L, 3L, 5L), verbose = FALSE)

  expect_s3_class(c2, "ld_region_c2")
  expect_equal(c2$n_cells, 18L)
  expect_true(c2$n_usable > 0 && c2$n_usable <= c2$n_cells)
  expect_true(all(c2$regions$c2 >= 0 & c2$regions$c2 <= 1))
  expect_equal(nrow(c2$landscape), 18L)
  # the planted region is significant across the grid, so it scores at the top
  expect_equal(c2$regions$c2[1], 1)
})

test_that("ld_scan reproduces the step-by-step result", {
  cs <- make_case()
  fit <- ld_scan(cs$p_obs, cs$p_perm, cs$ld_ws, cs$map, cs$GTs, cs$decay_sum,
                 tau = 0.05, l_min = 3L, rho_ld = 0.6, c2 = TRUE,
                 tau_grid = seq(0.05, 0.3, by = 0.05), lmin_grid = c(2L, 3L, 5L),
                 verbose = FALSE)
  null <- ld_null_from_p(cs$p_obs, cs$p_perm, cs$ld_ws, verbose = FALSE)
  s <- ld_region_scan(null, edges_for(cs, null), tau = 0.05, l_min = 3L)

  expect_s3_class(fit, "ld_scan")
  expect_equal(fit$regions$regions$s_R, s$regions$s_R)
  expect_equal(fit$regions$regions$p_R, s$regions$p_R)
  expect_true(fit$gate$pass)
  expect_false(is.null(fit$c2))

  # c2 = FALSE skips the expensive grid
  quick <- ld_scan(cs$p_obs, cs$p_perm, cs$ld_ws, cs$map, cs$GTs, cs$decay_sum,
                   tau = 0.05, l_min = 3L, rho_ld = 0.6, c2 = FALSE, verbose = FALSE)
  expect_null(quick$c2)
})

test_that("ld_scan requires a well-formed map", {
  cs <- make_case()
  bad_map <- data.table::copy(cs$map)[, marker := NULL]
  expect_error(ld_scan(cs$p_obs, cs$p_perm, cs$ld_ws, bad_map, cs$GTs, cs$decay_sum,
                       verbose = FALSE), "needs columns")
})

test_that("ld_cscore tolerates NA p-values", {
  set.seed(1)
  n <- 200L
  ld_ws <- matrix(runif(n * 3), n, 3,
                  dimnames = list(paste0("m", seq_len(n)), c("0.5", "0.7", "0.9")))
  p <- runif(n); p[c(3L, 17L, 200L)] <- NA_real_
  C <- expect_silent(ld_cscore(p, ld_ws))
  expect_length(C, n)
  expect_true(all(C[c(3L, 17L, 200L)] == 0))
  expect_true(all(C >= 0 & C <= 1))
})
