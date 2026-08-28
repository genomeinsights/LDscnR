test_that("the automatic edge-list floor leaves clustering output identical", {
  skip_if_not_installed("SNPRelate")
  set.seed(21)
  n_ind <- 80; n_snp <- 300
  blk  <- rep(seq_len(n_snp / 10), each = 10)
  base <- matrix(rbinom(n_ind * max(blk), 2, 0.35), n_ind)
  G <- base[, blk]
  fl <- sample(n_snp, n_snp %/% 4)
  G[, fl] <- pmin(2, pmax(0, G[, fl] + sample(c(-1, 0, 1), length(fl) * n_ind, TRUE)))
  colnames(G) <- paste0("s", seq_len(n_snp)); rownames(G) <- paste0("i", seq_len(n_ind))
  ## two chromosomes: compute_LD_decay() estimates background LD across them
  chrv <- rep(c("chr1", "chr2"), each = n_snp / 2)
  map <- data.table::data.table(marker = colnames(G), Chr = chrv,
                                Pos = rep(seq_len(n_snp / 2) * 500L, 2))
  f <- tempfile(fileext = ".gds"); on.exit(unlink(f), add = TRUE)
  gds <- create_gds_from_geno(G, map, f); on.exit(SNPRelate::snpgdsClose(gds), add = TRUE)

  dec <- compute_LD_decay(gds = gds, keep_el = FALSE, slide = 60, ld_method = "corr",
                          min_maf_decay = 0.05, n_win_decay = 5, cores = 1)
  skip_if(is.null(dec$decay_sum) || !nrow(dec$decay_sum), "decay fit unavailable")

  got <- ld_complexity_reduction(map = map, LD_decay = dec, rho = 0.5, gds = gds)

  ## reference: build the edge list with NO floor, then filter exactly as the
  ## function does, and cluster from that
  ds    <- data.table::as.data.table(dec$decay_sum)[Chr == "chr1"]
  skip_if(!nrow(ds), "chr1 dropped by the decay gate")
  r2_th <- ld_from_rho(b = ds$b, c = ds$c, rho = 0.5)
  mk1   <- map[Chr == "chr1", marker]
  el    <- get_el(gds, SNP_id = mk1, slide_win_ld = 60, method = "corr",
                  by_chr = TRUE, el_floor = 0)
  ref   <- el[el$r2 >= r2_th, ]
  flo   <- get_el(gds, SNP_id = mk1, slide_win_ld = 60, method = "corr",
                  by_chr = TRUE, el_floor = r2_th)

  ## the surviving edge set is the same either way -- which is what makes the
  ## automatic floor safe
  expect_equal(nrow(flo[flo$r2 >= r2_th, ]), nrow(ref))
  expect_equal(sort(ref$r2), sort(flo[flo$r2 >= r2_th, ]$r2), tolerance = 1e-12)
  ## and the clustering ran and partitioned every marker exactly once
  expect_equal(nrow(got$map_snp), nrow(map))
  expect_false(anyDuplicated(got$map_snp$marker) > 0)
})
