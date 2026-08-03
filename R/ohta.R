## =========================================================
## Fast Ohta D stats for unphased 0/1/2 diploid data
## =========================================================
## Moved into LDscnR (2026-08) from the formica_hybrid pipeline's local
## `Ohta.R`, so the LD/Ohta statistics live in one exported place. Function
## bodies are unchanged from the validated pipeline version.

#' Per-locus among-population R-squared
#'
#' Proportion of dosage variance explained by population membership
#' (`lm(g ~ pop)`), for a single locus.
#'
#' @param g Numeric vector of 0/1/2 allele dosages, one per individual.
#' @param pop Population labels, same length as `g`.
#' @return The R-squared, or `NA` when there is too little variation.
#' @importFrom stats lm
#' @export
calc_R2_pop_one <- function(g, pop) {
  pop <- as.factor(pop)
  ok <- !is.na(g) & !is.na(pop)

  if (sum(ok) < 3L || length(unique(g[ok])) < 2L || nlevels(droplevels(pop[ok])) < 2L) {
    return(NA_real_)
  }

  summary(lm(g[ok] ~ pop[ok]))$r.squared
}

#' Per-locus F_ST from 0/1/2 dosages
#'
#' Nei-style \eqn{F_{ST} = (H_T - H_S)/H_T} for a single locus, from allele
#' dosages, weighting subpopulations by sample size.
#'
#' @param g Numeric vector of 0/1/2 allele dosages, one per individual.
#' @param pop Population labels, same length as `g`.
#' @return F_ST clamped to \[0, 1\], or `NA` when undefined.
#' @export
calc_Fst_one <- function(g, pop) {
  pop <- as.factor(pop)
  ok <- !is.na(g) & !is.na(pop)

  g <- g[ok]
  pop <- droplevels(pop[ok])

  if (length(g) < 3L || length(unique(g)) < 2L || nlevels(pop) < 2L) {
    return(NA_real_)
  }

  ## dosage allele frequency
  p_total <- mean(g) / 2
  Ht <- 2 * p_total * (1 - p_total)

  if (!is.finite(Ht) || Ht <= 0) return(NA_real_)

  pop_stats <- tapply(g, pop, function(x) {
    p <- mean(x, na.rm = TRUE) / 2
    n <- sum(!is.na(x))
    c(n = n, Hs = 2 * p * (1 - p))
  })

  pop_stats <- do.call(rbind, pop_stats)
  w <- pop_stats[, "n"] / sum(pop_stats[, "n"])
  Hs <- sum(w * pop_stats[, "Hs"], na.rm = TRUE)

  Fst <- (Ht - Hs) / Ht
  max(0, min(1, Fst))
}

#' Prepare a genotype matrix for fast Ohta D-statistics
#'
#' Precomputes the per-population and total per-locus mean dosages (and keeps the
#' per-population submatrices) so that many locus-pair D-statistics can be scanned
#' without repeating the population split.
#'
#' @param data_set Individuals \eqn{\times} loci matrix of 0/1/2 dosages.
#' @param pops Population label per individual (length `nrow(data_set)`).
#' @return A `prep` list consumed by [dstat_unphased_fast()] /
#'   [dstat_unphased_scan()].
#' @export
ohta_fast_prepare <- function(data_set, pops) {
  G <- as.matrix(data_set)
  storage.mode(G) <- "numeric"

  if (length(pops) != nrow(G)) {
    stop("length(pops) must equal nrow(data_set)")
  }

  pop_levels <- unique(pops)
  pop_index <- split(seq_len(nrow(G)), pops)

  pop_mats <- lapply(pop_index, function(idx) G[idx, , drop = FALSE])
  pop_n <- vapply(pop_mats, nrow, integer(1))

  ## per-pop per-locus mean dosage
  pop_means <- do.call(rbind, lapply(pop_mats, function(M) {
    colMeans(M, na.rm = TRUE)
  }))

  ## total per-locus mean dosage
  total_means <- colMeans(G, na.rm = TRUE)

  list(
    G = G,
    pops = pops,
    pop_levels = pop_levels,
    pop_index = pop_index,
    pop_mats = pop_mats,
    pop_n = pop_n,
    pop_means = pop_means,
    total_means = total_means
  )
}

## joint genotype counts for one population and one pair of loci
.get_joint_counts_fast <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]

  n <- length(x)
  if (n == 0) return(list(n = 0L, s = rep(0L, 9)))

  s <- c(
    sum(x == 0 & y == 0),  # s00
    sum(x == 0 & y == 1),  # s01
    sum(x == 0 & y == 2),  # s02
    sum(x == 1 & y == 0),  # s10
    sum(x == 1 & y == 1),  # s11
    sum(x == 1 & y == 2),  # s12
    sum(x == 2 & y == 0),  # s20
    sum(x == 2 & y == 1),  # s21
    sum(x == 2 & y == 2)   # s22
  )

  list(n = n, s = s)
}

#' Fast Ohta two-locus D-statistics for one locus pair (unphased 0/1/2)
#'
#' Computes the Ohta (1982) among/within-population two-locus disequilibrium
#' decomposition (`D2it`, `D2is`, `D2st`, `Dp2st`, `Dp2is`) plus a global dosage
#' `r2`, per-population mean `R2_pop_mean` and `Fst_mean`, from a `prep` object.
#'
#' @param index Length-2 integer vector: the two locus (column) indices.
#' @param prep Output of [ohta_fast_prepare()].
#' @param tot_maf Total-sample minor-allele-frequency floor (default 0.1).
#' @param pop_maf Per-population minor-allele-frequency floor (default 0.05).
#' @return A named numeric vector of the statistics (or `NA`s when the pair fails
#'   the MAF filters).
#' @importFrom stats sd cor complete.cases
#' @export
dstat_unphased_fast <- function(index, prep, tot_maf = 0.1, pop_maf = 0.05) {
  if (length(index) != 2) stop("index must have length 2")
  i <- index[1]
  j <- index[2]

  na_out <- c(
    nPops = NA,
    r2 = NA,
    R2_pop_mean = NA,
    Fst_mean = NA,
    D2it = NA,
    D2is = NA,
    D2st = NA,
    Dp2st = NA,
    Dp2is = NA
  )

  tot_maf2 <- tot_maf * 2
  pop_maf2 <- pop_maf * 2
  tot_max <- 2 - tot_maf2
  pop_max <- 2 - pop_maf2

  m1 <- prep$total_means[i]
  m2 <- prep$total_means[j]

  if (!is.finite(m1) || !is.finite(m2) ||
      m1 < tot_maf2 || m1 > tot_max ||
      m2 < tot_maf2 || m2 > tot_max) {
    return(na_out)
  }

  pm1 <- prep$pop_means[, i]
  pm2 <- prep$pop_means[, j]

  keep <- is.finite(pm1) & is.finite(pm2) &
    pm1 >= pop_maf2 & pm1 <= pop_max &
    pm2 >= pop_maf2 & pm2 <= pop_max

  if (!any(keep)) return(na_out)

  mats <- prep$pop_mats[keep]

  kept_pops <- prep$pop_levels[keep]
  keep_ind <- prep$pops %in% kept_pops

  g1 <- prep$G[keep_ind, i]
  g2 <- prep$G[keep_ind, j]
  pops_pair <- prep$pops[keep_ind]

  R2_pop_1 <- calc_R2_pop_one(g1, pops_pair)
  R2_pop_2 <- calc_R2_pop_one(g2, pops_pair)
  R2_pop_mean <- mean(c(R2_pop_1, R2_pop_2), na.rm = TRUE)
  if (!is.finite(R2_pop_mean)) R2_pop_mean <- NA_real_

  Fst_1 <- calc_Fst_one(g1, pops_pair)
  Fst_2 <- calc_Fst_one(g2, pops_pair)
  Fst_mean <- mean(c(Fst_1, Fst_2), na.rm = TRUE)
  if (!is.finite(Fst_mean)) Fst_mean <- NA_real_

  nPops <- length(mats)

  ## global r2 after same population filtering
  geno_filt <- do.call(rbind, lapply(mats, function(M) M[, c(i, j), drop = FALSE]))
  geno_cc <- geno_filt[complete.cases(geno_filt), , drop = FALSE]

  if (nrow(geno_cc) >= 3 &&
      sd(geno_cc[, 1]) > 0 &&
      sd(geno_cc[, 2]) > 0) {
    r2 <- cor(geno_cc[, 1], geno_cc[, 2])^2
  } else {
    r2 <- NA_real_
  }

  Tijs <- matrix(NA_real_, nrow = nPops, ncol = 4)
  Pf   <- matrix(NA_real_, nrow = nPops, ncol = 4)

  for (k in seq_len(nPops)) {
    M <- mats[[k]]
    x <- M[, i]
    y <- M[, j]

    ok <- !is.na(x) & !is.na(y)
    x2 <- x[ok]
    y2 <- y[ok]
    n <- length(x2)

    if (n == 0) next

    p1.0 <- 1 - sum(x2) / (2 * n)
    p1.2 <- 1 - p1.0
    p2.0 <- 1 - sum(y2) / (2 * n)
    p2.2 <- 1 - p2.0

    Pf[k, ] <- c(p1.0, p1.2, p2.0, p2.2)

    s00 <- sum(x2 == 0 & y2 == 0)
    s01 <- sum(x2 == 0 & y2 == 1)
    s02 <- sum(x2 == 0 & y2 == 2)
    s10 <- sum(x2 == 1 & y2 == 0)
    s11 <- sum(x2 == 1 & y2 == 1)
    s12 <- sum(x2 == 1 & y2 == 2)
    s20 <- sum(x2 == 2 & y2 == 0)
    s21 <- sum(x2 == 2 & y2 == 1)
    s22 <- sum(x2 == 2 & y2 == 2)

    Tijs[k, ] <- c(
      (2*s00 + s01 + s10 + 0.5*s11) / n,
      (2*s02 + s01 + s12 + 0.5*s11) / n,
      (2*s20 + s10 + s21 + 0.5*s11) / n,
      (2*s22 + s21 + s12 + 0.5*s11) / n
    )
  }

  good <- rowSums(is.na(Tijs)) == 0 & rowSums(is.na(Pf)) == 0
  Tijs <- Tijs[good, , drop = FALSE]
  Pf   <- Pf[good, , drop = FALSE]
  nPops <- nrow(Tijs)

  if (nPops == 0) return(na_out)

  P1m <- cbind(
    2 * Pf[,1] * Pf[,3],
    2 * Pf[,1] * Pf[,4],
    2 * Pf[,2] * Pf[,3],
    2 * Pf[,2] * Pf[,4]
  )

  P2m <- cbind(
    Pf[,1] * Pf[,3],
    Pf[,1] * Pf[,4],
    Pf[,2] * Pf[,3],
    Pf[,2] * Pf[,4]
  )

  Pf_mean <- colMeans(Pf)

  Pf_total_manip1 <- matrix(
    c(
      2*Pf_mean[1]*Pf_mean[3],
      2*Pf_mean[1]*Pf_mean[4],
      2*Pf_mean[2]*Pf_mean[3],
      2*Pf_mean[2]*Pf_mean[4]
    ),
    nrow = nPops, ncol = 4, byrow = TRUE
  )

  Pf_total_manip2 <- matrix(
    c(
      Pf_mean[1]*Pf_mean[3],
      Pf_mean[1]*Pf_mean[4],
      Pf_mean[2]*Pf_mean[3],
      Pf_mean[2]*Pf_mean[4]
    ),
    nrow = nPops, ncol = 4, byrow = TRUE
  )

  Tbar <- matrix(colMeans(Tijs), nrow = nPops, ncol = 4, byrow = TRUE)

  D2is  <- sum((Tijs - P1m)^2) / nPops
  D2it  <- sum((Tijs - Pf_total_manip1)^2) / nPops
  D2st  <- sum((P2m - Pf_total_manip2)^2) / nPops
  Dp2st <- sum((colMeans(Tijs) - Pf_total_manip1[1, ])^2) / 4
  Dp2is <- sum((Tijs - Tbar)^2) / nPops / 4

  round(c(
    nPops = nPops,
    r2 = r2,
    R2_pop_mean = R2_pop_mean,
    Fst_mean = Fst_mean,
    D2it = D2it,
    D2is = D2is,
    D2st = D2st,
    Dp2st = Dp2st,
    Dp2is = Dp2is
  ), 6)
}

#' Scan many locus pairs for Ohta D-statistics
#'
#' Applies [dstat_unphased_fast()] over many locus pairs, optionally in parallel.
#'
#' @param pairs Two-column matrix/data.frame of locus index pairs.
#' @param prep Output of [ohta_fast_prepare()].
#' @param tot_maf Total-sample MAF floor (default 0.1).
#' @param pop_maf Per-population MAF floor (default 0.05).
#' @param cores Number of cores for [parallel::mclapply()] (default 1).
#' @return A `data.table`, one row of statistics per input pair.
#' @importFrom parallel mclapply
#' @export
dstat_unphased_scan <- function(pairs, prep, tot_maf = 0.1, pop_maf = 0.05, cores = 1) {
  pairs <- as.matrix(pairs)
  if (ncol(pairs) != 2) stop("pairs must have 2 columns")

  out <- do.call(rbind, mclapply(seq_len(nrow(pairs)), function(ii) {
    data.table(t(dstat_unphased_fast(
      index = pairs[ii, ],
      prep = prep,
      tot_maf = tot_maf,
      pop_maf = pop_maf
    )))
  }, mc.cores = cores))

  out[]
}

#' Reference Ohta D-statistics for one locus pair
#'
#' The original (slower, name-keyed) Ohta two-locus decomposition, kept as a
#' cross-check for [dstat_unphased_fast()]. Populations are taken from
#' `rownames(data_set)`.
#'
#' @param index Length-2 integer vector: the two locus (column) indices.
#' @param data_set Individuals \eqn{\times} loci dosage matrix, row names = pops.
#' @param tot_maf Total-sample MAF floor (default 0.1).
#' @param pop_maf Per-population MAF floor (default 0.05).
#' @return A rounded named numeric vector `nPops, r2, D2it, D2is, D2st, Dp2st, Dp2is`.
#' @export
ohtadstats_dstat <- function (index, data_set, tot_maf = 0.1, pop_maf = 0.05)
{
  tot_maf = tot_maf * 2
  pop_maf = pop_maf * 2
  tot_max_thresh = 2 - tot_maf
  pop_max_thresh = 2 - pop_maf
  geno <- data_set[, c(index[1], index[2])]
  if (mean(geno[, 1], na.rm = TRUE) >= tot_maf & mean(geno[, 2],
                                                      na.rm = TRUE) >= tot_maf & mean(geno[, 1], na.rm = TRUE) <=
      tot_max_thresh & mean(geno[, 2], na.rm = TRUE) <= tot_max_thresh) {
    freqs1 <- unlist(by(geno[, 1], rownames(geno), mean,
                        na.rm = TRUE))
    freqs2 <- unlist(by(geno[, 2], rownames(geno), mean,
                        na.rm = TRUE))
    rm <- c(names(freqs1)[which(freqs1 < pop_maf | freqs1 >
                                  pop_max_thresh)], names(freqs2)[which(freqs2 < pop_maf |
                                                                          freqs2 > pop_max_thresh)])
    if (length(rm) > 0)
      geno <- geno[-which(rownames(geno) %in% rm), ]
    if (nrow(geno) > 0) {
      nPops <- length(table(rownames(geno)))
      T <- function(geno) {
        geno <- geno[which(is.na(geno[, 1]) == F & is.na(geno[, 2]) == FALSE), ]
        length <- nrow(geno)
        s.00 <- length(which(geno[, 1] == 0 & geno[,2] == 0))
        s.01 <- length(which(geno[, 1] == 0 & geno[,2] == 1))
        s.02 <- length(which(geno[, 1] == 0 & geno[,2] == 2))
        s.10 <- length(which(geno[, 1] == 1 & geno[,2] == 0))
        s.11 <- length(which(geno[, 1] == 1 & geno[,2] == 1))
        s.12 <- length(which(geno[, 1] == 1 & geno[,2] == 2))
        s.20 <- length(which(geno[, 1] == 2 & geno[,2] == 0))
        s.21 <- length(which(geno[, 1] == 2 & geno[,2] == 1))
        s.22 <- length(which(geno[, 1] == 2 & geno[,2] == 2))
        T00 <- (2 * s.00 + s.01 + s.10 + 0.5 * s.11)/length
        T02 <- (2 * s.02 + s.01 + s.12 + 0.5 * s.11)/length
        T20 <- (2 * s.20 + s.10 + s.21 + 0.5 * s.11)/length
        T22 <- (2 * s.22 + s.21 + s.12 + 0.5 * s.11)/length
        return(c(T00, T02, T20, T22))
      }
      P <- function(geno, out) {
        geno <- geno[which(is.na(geno[, 1]) == F & is.na(geno[,
                                                              2]) == FALSE), ]
        length <- nrow(geno)
        p1.0 <- 1 - sum(geno[, 1], na.rm = TRUE)/(2 * length)
        p1.2 <- 1 - p1.0
        p2.0 <- 1 - sum(geno[, 2], na.rm = TRUE)/(2 * length)
        p2.2 <- 1 - p2.0
        if (out == "freq")
          return(c(p1.0, p1.2, p2.0, p2.2))
        if (out == "manip1")
          return(c(2 * p1.0 * p2.0, 2 * p1.0 * p2.2,
                   2 * p1.2 * p2.0, 2 * p1.2 * p2.2))
        if (out == "manip2")
          return(c(p1.0 * p2.0, p1.0 * p2.2, p1.2 * p2.0,
                   p1.2 * p2.2))
      }
      Tijs <- by(geno, rownames(geno), T)
      Tijs.n <- as.numeric(unlist(Tijs))
      P1m <- by(geno, rownames(geno), P, "manip1")
      P1m.n <- as.numeric(unlist(P1m))
      P2m <- by(geno, rownames(geno), P, "manip2")
      P2m.n <- as.numeric(unlist(P2m))
      Pf <- by(geno, rownames(geno), P, "freq")
      Pf.n <- as.numeric(unlist(Pf))
      n.pops <- length(levels(as.factor(rownames(geno))))
      D2is <- sum((Tijs.n - P1m.n)^2)/n.pops
      p1.0 <- mean(Pf.n[seq(1, length(Pf.n), 4)])
      p1.2 <- mean(Pf.n[seq(2, length(Pf.n), 4)])
      p2.0 <- mean(Pf.n[seq(3, length(Pf.n), 4)])
      p2.2 <- mean(Pf.n[seq(4, length(Pf.n), 4)])
      Pf.v <- c()
      Pf.v[seq(1, length(Pf.n), 4)] <- 2 * p1.0 * p2.0
      Pf.v[seq(2, length(Pf.n), 4)] <- 2 * p1.0 * p2.2
      Pf.v[seq(3, length(Pf.n), 4)] <- 2 * p1.2 * p2.0
      Pf.v[seq(4, length(Pf.n), 4)] <- 2 * p1.2 * p2.2
      T00 <- mean(Tijs.n[seq(1, length(Tijs.n), 4)])
      T02 <- mean(Tijs.n[seq(2, length(Tijs.n), 4)])
      T20 <- mean(Tijs.n[seq(3, length(Tijs.n), 4)])
      T22 <- mean(Tijs.n[seq(4, length(Tijs.n), 4)])
      Tijs.n2 <- c()
      Tijs.n2[seq(1, length(Tijs.n), 4)] <- T00
      Tijs.n2[seq(2, length(Tijs.n), 4)] <- T02
      Tijs.n2[seq(3, length(Tijs.n), 4)] <- T20
      Tijs.n2[seq(4, length(Tijs.n), 4)] <- T22
      D2it <- sum((Tijs.n - Pf.v)^2)/n.pops
      D2st <- sum((P2m.n - Pf.v/2)^2)/n.pops
      p0p0 <- 2 * p1.0 * p2.0
      p0p2 <- 2 * p1.0 * p2.2
      p2p0 <- 2 * p1.2 * p2.0
      p2p2 <- 2 * p1.2 * p2.2
      Dp2st <- (T00 - p0p0)^2 + (T02 - p0p2)^2 + (T20 - p2p0)^2 + (T22 - p2p2)^2
      Dp2is <- sum((Tijs.n - Tijs.n2)^2)/n.pops
    }
    else {
      D2it <- NA
      D2is <- NA
      D2st <- NA
      Dp2st <- NA
      Dp2is <- NA
      nPops <- NA
      r2 <- NA
    }
  }
  else {
    D2it <- NA
    D2is <- NA
    D2st <- NA
    Dp2st <- NA
    Dp2is <- NA
    nPops <- NA
    r2 <- NA
  }

  if (nrow(geno) > 0) {
    nPops <- length(table(rownames(geno)))

    geno_cc <- geno[complete.cases(geno[, 1:2]), , drop = FALSE]

    if (nrow(geno_cc) >= 3 &&
        sd(geno_cc[, 1], na.rm = TRUE) > 0 &&
        sd(geno_cc[, 2], na.rm = TRUE) > 0) {
      r2 <- cor(geno_cc[, 1], geno_cc[, 2], use = "complete.obs")^2
      } else {
      r2 <- NA
      }
  }


  return(round(c(nPops, r2,D2it, D2is, D2st, Dp2st, Dp2is), 6))
}

#' Prepare and compute Ohta D-statistics for one locus pair
#'
#' Convenience wrapper: builds a `prep` with [ohta_fast_prepare()] and returns
#' [dstat_unphased_fast()] for a single pair.
#'
#' @param index Length-2 integer vector: the two locus (column) indices.
#' @param geno Individuals \eqn{\times} loci dosage matrix.
#' @param pops Population label per individual.
#' @param tot_maf Total-sample MAF floor (default 0.1).
#' @param pop_maf Per-population MAF floor (default 0.05).
#' @return A named numeric vector of the statistics.
#' @export
run_ohtadstats <- function(index, geno, pops, tot_maf = 0.1, pop_maf = 0.05) {
  X <- geno
  rownames(X) <- pops

  prep <- ohta_fast_prepare(data_set = X, pops = pops)

  out <- dstat_unphased_fast(
    index = index,
    prep = prep,
    tot_maf = tot_maf,
    pop_maf = pop_maf
  )

  out <- as.numeric(out)
  names(out) <- c(
    "nPops", "r2", "R2_pop_mean", "Fst_mean",
    "D2it", "D2is", "D2st", "Dp2st", "Dp2is"
  )
  out
}
