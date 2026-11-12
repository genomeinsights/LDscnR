# Extracted from test-ld_unit_matrix.R:97

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "LDscnR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
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
  GTs <- cbind(CL1, CL2, CL3, CL4)
  rownames(GTs) <- paste0("ind", seq_len(n_ind))

  map_snp <- data.table::data.table(
    marker = colnames(GTs), Chr = "Chr1",
    Pos = c(1000:1002, 1010:1012, 5000:5002, 2000:2001),
    CL_id = rep(1:4, c(3, 3, 3, 2)),
    n_loci = rep(c(3, 3, 3, 2), c(3, 3, 3, 2)),
    is_core = unlist(list(c(TRUE, FALSE, FALSE), c(TRUE, FALSE, FALSE),
                          c(TRUE, FALSE, FALSE), c(TRUE, FALSE))),
    median_ld = NA_real_,
    ld_w_095 = c(0.9, 0.05, 0.05,  0.9, 0.9, 0.9,  0.9, 0.9, 0.9,  0.05, 0.05)
  )
  clusters <- data.table::data.table(
    Chr = "Chr1", CL_id = 1:4,
    core_snp = c("cl1_1", "cl2_1", "cl3_1", "cl4_1"),
    median_ld = c(0.9, 0.9, 0.9, NA),
    n_snps = c(3, 3, 3, 2),
    members = list(colnames(CL1), colnames(CL2), colnames(CL3), colnames(CL4))
  )
  map <- map_snp[, .(marker, Chr, Pos)]
  list(GTs = GTs, map = map, stage1 = list(map_snp = map_snp, clusters = clusters))
}

# test -------------------------------------------------------------------------
d <- build_stage1()
set.seed(2)
G <- d$GTs[, sample(ncol(d$GTs))]
m <- ld_unit_matrix(G, d$stage1, d$map, size_floor = 2L, repr = "consensus_dosage")
units <- attr(m, "units")
mem <- d$stage1$clusters$members[match(units$unit_id, d$stage1$clusters$CL_id)]
ref <- vapply(mem, function(mk) consensus_dosage(G, mk), numeric(nrow(G)))
expect_identical(unname(m), unname(ref))
