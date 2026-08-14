## A causal LD block (associated) and a decoy LD block (correlated but not
## associated) both pass the ld_w candidate filter; the outlier caller should
## flag the causal block and NOT the decoy, and the background null should agree.

make_toy <- function(seed = 1) {
  set.seed(seed)
  n <- 140
  ## two correlated blocks (high LD) + independent background (low LD)
  latent_causal <- sample(0:2, n, replace = TRUE, prob = c(.25, .5, .25))
  ## phenotype driven by the causal latent only
  Y <- latent_causal + rnorm(n)
  ## decoy latent constructed to be genuinely null (orthogonal to Y), so the
  ## specificity test does not depend on a chance block-phenotype correlation
  repeat {
    latent_decoy <- sample(0:2, n, replace = TRUE, prob = c(.25, .5, .25))
    if (abs(stats::cor(latent_decoy, Y)) < 0.03) break
  }
  block <- function(latent, m) {              # m tightly-correlated dosages
    vapply(seq_len(m), function(j) {
      x <- latent; flip <- sample(n, round(0.05 * n))
      x[flip] <- sample(0:2, length(flip), replace = TRUE); x
    }, numeric(n))
  }
  causal <- block(latent_causal, 25)          # chr A, associated
  decoy  <- block(latent_decoy, 25)           # chr B, NOT associated
  bg     <- matrix(sample(0:2, n * 150, replace = TRUE), n, 150)  # low-LD background
  GTs <- cbind(causal, decoy, bg)
  colnames(GTs) <- paste0("m", seq_len(ncol(GTs)))
  map <- data.frame(
    marker = colnames(GTs),
    Chr = c(rep("A", 25), rep("B", 25), rep(c("A", "B"), each = 75)),
    Pos = c(seq_len(25) * 1e4, seq_len(25) * 1e4,
            25e4 + seq_len(75) * 1e4, 25e4 + seq_len(75) * 1e4),
    stringsAsFactors = FALSE
  )
  ## per-SNP association p-values (simple linear model; method-agnostic input)
  pval <- apply(GTs, 2, function(x) {
    if (stats::sd(x) == 0) return(1)
    stats::cor.test(Y, x)$p.value
  })
  ## ld_w: high for the two blocks, low for background
  ld_w <- c(runif(50, 0.6, 0.9), runif(150, 0, 0.2))
  names(pval) <- names(ld_w) <- colnames(GTs)
  list(GTs = GTs, map = map, pval = pval, ld_w = ld_w, Y = Y)
}

test_that("rmsc_threshold returns a valid selection curve", {
  d <- make_toy()
  s <- rmsc_threshold(d$pval, d$ld_w, fdr = 0.1)
  expect_named(s, c("grid", "rejections", "q_star", "no_elbow"))
  expect_length(s$rejections, length(s$grid))
  expect_true(all(s$rejections >= 0))
  expect_true(s$q_star >= 0 && s$q_star <= max(s$grid))
  ## when ld_w carries no information (constant), filtering cannot help -> no elbow
  flat <- rmsc_threshold(d$pval, rep(1, length(d$pval)), fdr = 0.1)
  expect_true(flat$no_elbow)
})

test_that("background null flags the causal block and not the decoy", {
  d <- make_toy()
  set.seed(42)
  r <- ld_outlier_clusters(d$pval, d$ld_w, d$map, d$GTs,
                           null = "background", fdr = 0.1, B = 499, verbose = FALSE)
  expect_s3_class(r, "ld_outlier_clusters")
  sigc <- r$clusters[significant == TRUE]
  ## at least one significant cluster, all on chr A (the causal block)
  expect_gt(nrow(sigc), 0)
  expect_true(all(sigc$Chr == "A"))
  ## the decoy chromosome B is never flagged
  expect_false(any(r$clusters[Chr == "B", significant]))
  ## the flagged cluster is built from causal markers (m1..m25)
  causal_markers <- paste0("m", 1:25)
  expect_true(mean(unlist(sigc$members) %in% causal_markers) > 0.8)
})

test_that("permutation null (fast EMMAX path) flags the causal block", {
  d <- make_toy()
  K <- tcrossprod(scale(d$GTs)) / ncol(d$GTs)          # simple genomic relationship
  set.seed(7)
  r <- ld_outlier_clusters(d$pval, d$ld_w, d$map, d$GTs, null = "permutation",
                           Y = d$Y, K = K, fdr = 0.1, B = 299, verbose = FALSE)
  sigc <- r$clusters[significant == TRUE]
  expect_gt(nrow(sigc), 0)
  expect_true(all(sigc$Chr == "A"))
  expect_false(any(r$clusters[Chr == "B", significant]))
})

test_that("no significant SNPs yields no flagged clusters without error", {
  d <- make_toy()
  flat <- stats::setNames(rep(0.9, length(d$pval)), names(d$pval))  # nothing significant
  expect_warning(
    r <- ld_outlier_clusters(flat, d$ld_w, d$map, d$GTs,
                             null = "background", B = 99, verbose = FALSE),
    "No SNP significant"
  )
  expect_equal(sum(r$clusters$significant, na.rm = TRUE), 0)
})
