test_that("el_floor filters during construction and defaults to no-op", {
  skip_if_not_installed("SNPRelate")
  set.seed(11)
  n_ind <- 60; n_snp <- 120
  ## a genotype matrix with real local LD: blocks of correlated markers
  block <- rep(seq_len(n_snp / 6), each = 6)
  base  <- matrix(rbinom(n_ind * max(block), 2, 0.4), n_ind)
  G <- base[, block] 
  flip <- sample(seq_len(n_snp), n_snp %/% 3)
  G[, flip] <- pmin(2, pmax(0, G[, flip] + sample(c(-1, 0, 1), length(flip) * n_ind, TRUE)))
  colnames(G) <- paste0("s", seq_len(n_snp))
  rownames(G) <- paste0("i", seq_len(n_ind))
  map <- data.table::data.table(marker = colnames(G), Chr = "chr1",
                                Pos = seq_len(n_snp) * 1000L)
  f <- tempfile(fileext = ".gds"); on.exit(unlink(f), add = TRUE)
  gds <- create_gds_from_geno(G, map, f); on.exit(SNPRelate::snpgdsClose(gds), add = TRUE)

  full <- get_el(gds, slide_win_ld = 30, method = "corr")
  expect_gt(nrow(full), 0)

  ## default is exactly the old behaviour
  expect_identical(get_el(gds, slide_win_ld = 30, method = "corr", el_floor = 0), full)

  ## a floor returns precisely the rows the caller would have kept anyway
  for (thr in c(0.05, 0.2, 0.5)) {
    got  <- get_el(gds, slide_win_ld = 30, method = "corr", el_floor = thr)
    want <- full[full$r2 >= thr, ]
    expect_equal(nrow(got), nrow(want))
    expect_equal(sort(got$r2), sort(want$r2), tolerance = 1e-12)
    expect_true(all(got$r2 >= thr))
  }

  ## and it really does shrink the object, which is the point
  expect_lt(nrow(get_el(gds, slide_win_ld = 30, method = "corr", el_floor = 0.5)), nrow(full))

  ## by_chr path honours it too
  expect_true(all(get_el(gds, slide_win_ld = 30, method = "corr",
                         by_chr = TRUE, el_floor = 0.3)$r2 >= 0.3))
})

test_that("el_floor rejects nonsense", {
  expect_error(get_el(NULL, el_floor = -0.1), "el_floor")
  expect_error(get_el(NULL, el_floor = 1.5), "el_floor")
  expect_error(get_el(NULL, el_floor = c(0.1, 0.2)), "el_floor")
  expect_error(get_el(NULL, el_floor = NA_real_), "el_floor")
})
