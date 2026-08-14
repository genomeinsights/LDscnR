## Candidate-first LD-aware outlier detection.
##
## An inversion of the two-stage complexity-reduction pipeline for outlier /
## association scans: instead of clustering the whole genome and testing every
## cluster, we (1) pick a data-driven local-LD threshold that maximises
## discoveries, (2) cluster ONLY the high-ld_w candidates, (3) test single SNPs,
## and (4) let each significant SNP recruit its LD cluster, which becomes the
## reporting / test unit. A size-aware cluster-level null (permutation or
## below-threshold background) decides which clusters are real.

#' Rejection-maximising selection of a local-LD threshold
#'
#' Independent-filtering (Bourgon et al. 2010) choice of an `ld_w` quantile
#' threshold: over a grid of quantiles, retain the SNPs above each, FDR-correct
#' their association p-values, and return the threshold that maximises the number
#' of discoveries. Because `ld_w` is computed from LD structure alone it is
#' approximately independent of the association p-value under the null, which is
#' the condition independent filtering requires.
#'
#' @param pval Numeric vector of per-SNP association p-values (from any method).
#' @param ld_w Numeric vector of per-SNP local LD support ([compute_ld_w()]),
#'   aligned to `pval`.
#' @param fdr SNP-level false-discovery rate used when counting discoveries.
#' @param grid Grid of `ld_w` quantile thresholds to scan.
#'
#' @return A list with `grid`, `rejections` (discoveries at each threshold),
#'   `q_star` (the rejection-maximising quantile) and `no_elbow` (TRUE if
#'   filtering never beats keeping all SNPs, i.e. the prioritisation does not
#'   apply to these data).
#' @references Bourgon R, Gentleman R, Huber W (2010). Independent filtering
#'   increases detection power for high-throughput experiments. \emph{PNAS}
#'   107:9546-9551.
#' @export
rmsc_threshold <- function(pval, ld_w, fdr = 0.05, grid = seq(0, 0.98, by = 0.01)) {
  ok <- is.finite(pval) & is.finite(ld_w)
  rej <- vapply(grid, function(q) {
    thr <- stats::quantile(ld_w[ok], q, na.rm = TRUE)
    k <- ok & ld_w >= thr
    sum(stats::p.adjust(pval[k], "BH") < fdr, na.rm = TRUE)
  }, integer(1))
  qs <- grid[which.max(rej)]
  list(grid = grid, rejections = rej, q_star = qs,
       no_elbow = qs <= min(grid) + 1e-9 || max(rej) <= rej[which.min(abs(grid))])
}

#' Candidate-first LD-aware outlier clusters
#'
#' Detect outlier / association regions by clustering only the high local-LD
#' candidates and testing each cluster against a size-aware null. The steps are:
#'
#' 1. **Threshold.** [rmsc_threshold()] picks a data-driven `ld_w` quantile
#'    `q*` (or supply `ld_w_threshold`). Candidates are SNPs with
#'    `ld_w >= quantile(ld_w, q*)` -- only these can be outliers.
#' 2. **Cluster.** Single-linkage clustering of the candidates on `1 - r^2`
#'    within each chromosome (connected components at `r2_link`); a SNP joins a
#'    cluster if it is in LD with any member.
#' 3. **Test SNPs.** FDR within candidates flags significant SNPs.
#' 4. **Recruit & test clusters.** Each cluster (significant SNPs plus the LD
#'    neighbours recruited in step 2) is the reporting unit. A cluster with at
#'    least `min_cluster_sig` significant SNPs is tested against a size-aware
#'    null of its significant-SNP count.
#'
#' Two nulls are available. `"background"` draws the null count from the
#' below-threshold (low-`ld_w`) SNPs -- the genomic background -- size-matched to
#' each cluster; it needs only `pval`, so it works with **any** association
#' method (LFMM, pcadapt, BayPass XtX, ...), but it treats a cluster's SNPs as
#' independent draws and is therefore anti-conservative in magnitude for large
#' correlated blocks (e.g. inversions). `"permutation"` re-runs an EMMAX scan on
#' the candidates under permuted phenotypes (requires `Y` and `K`); it captures
#' within-cluster LD and structure confounding but only applies when the
#' relatedness model `K` is independent of the phenotype (mixed models, PC/Q
#' covariates, BayPass with a fixed \eqn{\Omega}).
#'
#' @param pval Per-SNP association p-values (any method). If named, reordered to
#'   `map$marker`; otherwise assumed aligned to `map` rows.
#' @param ld_w Per-SNP local LD support ([compute_ld_w()]), aligned as `pval`.
#' @param map Data frame / data.table with `marker`, `Chr`, `Pos`.
#' @param GTs Genotype dosage matrix, individuals (rows) by markers (columns);
#'   `colnames(GTs)` must contain every `map$marker`. Used for candidate LD.
#' @param ld_w_threshold Optional fixed `ld_w` quantile in `(0, 1)`. `NULL`
#'   (default) selects it with [rmsc_threshold()].
#' @param fdr SNP-level FDR within candidates.
#' @param r2_link Single-linkage \eqn{r^2} threshold for candidate clustering.
#'   Ignored on a chromosome for which `rho_ld` supplies a value.
#' @param distance_threshold Optional maximum consecutive-position gap (bp) within
#'   a cluster: after single-linkage, a cluster is split wherever adjacent members
#'   are farther apart than this (a distance-restricted single linkage, as in
#'   [ld_prune_and_eMLG()]'s `split_by_distance`). Stops long-range LD from
#'   chaining physically distant signals through a coherent block. `NULL` (default)
#'   applies no distance cap unless `rho_d` is given.
#' @param rho_ld,rho_d Optional relative-LD thresholds that derive the linkage and
#'   distance cuts *per chromosome* from the LD-decay fit (requires `LD_decay`):
#'   `r2_link = ld_from_rho(b, c, rho_ld)` (higher `rho_ld` = looser \eqn{r^2}),
#'   `distance_threshold = d_from_rho(a_pred, rho_d)` (higher `rho_d` = wider gap).
#'   Either may be given alone; each overrides the corresponding direct argument
#'   on chromosomes present in the decay fit.
#' @param LD_decay Object from [compute_LD_decay()] (its `$decay_sum` supplies the
#'   per-chromosome `a_pred`, `b`, `c`); required only when `rho_ld`/`rho_d` are set.
#' @param min_cluster_sig Minimum significant SNPs for a cluster to be tested.
#' @param null Cluster-level null: `"background"` (default, method-agnostic) or
#'   `"permutation"` (EMMAX, needs `Y`/`K`).
#' @param cluster_fdr Cluster-level FDR (on the permutation/background `emp_p`).
#' @param B Number of null replicates.
#' @param Y,K Phenotype vector and kinship for `null = "permutation"`. `K` may be
#'   a single `n x n` matrix or a named list of per-chromosome (LOCO) matrices.
#' @param rmsc_grid Grid passed to [rmsc_threshold()].
#' @param cores Cores for the permutation loop.
#' @param verbose Emit diagnostic warnings: no RMSC elbow; every chromosome
#'   flagged; `ld_w` tracking genome-wide structure (positive `ld_w`-vs-signal
#'   correlation on most chromosomes); significant clusters spanning whole
#'   chromosomes. The last two flag structure-dominated runs whose calls are
#'   unreliable (see `$diagnostics$ldw_tracks_structure`, `$whole_chr_clusters`).
#'
#' @return An object of class `ld_outlier_clusters`: a list with `clusters`
#'   (per-cluster `Chr`, size `n`, `n_sig`, `start`/`end`, `members`, `emp_p`,
#'   `emp_q`, `significant`), `candidates` (per-candidate marker table),
#'   `ld_w_threshold`/`ld_w_value`, `pcut` (SNP significance cutoff), `rmsc`
#'   (the selection curve), and `diagnostics`.
#' @seealso [rmsc_threshold()], [compute_ld_w()], [emmax()]
#' @importFrom stats as.dist cor cutree hclust p.adjust quantile
#' @export
ld_outlier_clusters <- function(pval, ld_w, map, GTs,
                                ld_w_threshold = NULL, fdr = 0.05, r2_link = 0.5,
                                distance_threshold = NULL,
                                rho_ld = NULL, rho_d = NULL, LD_decay = NULL,
                                min_cluster_sig = 1L,
                                null = c("background", "permutation"),
                                cluster_fdr = 0.05, B = 1000L,
                                Y = NULL, K = NULL,
                                rmsc_grid = seq(0, 0.98, by = 0.01),
                                cores = 1L, verbose = TRUE) {
  null <- match.arg(null)
  map <- data.table::copy(data.table::as.data.table(map))
  if (!all(c("marker", "Chr", "Pos") %in% names(map)))
    stop("`map` must have columns marker, Chr, Pos.")
  if (!is.null(names(pval))) pval <- pval[map$marker]
  if (!is.null(names(ld_w))) ld_w <- ld_w[map$marker]
  if (length(pval) != nrow(map) || length(ld_w) != nrow(map))
    stop("`pval` and `ld_w` must align to `map` (length nrow(map) or named by marker).")
  if (is.null(colnames(GTs)) || !all(map$marker %in% colnames(GTs)))
    stop("`GTs` must have colnames covering every map$marker.")
  map[, `:=`(pv0 = as.numeric(pval), ldw0 = as.numeric(ld_w),
             Chr = as.character(Chr))]

  ## 1. ld_w threshold (data-driven unless supplied) --------------------------
  sel <- rmsc_threshold(map$pv0, map$ldw0, fdr = fdr, grid = rmsc_grid)
  qstar <- if (is.null(ld_w_threshold)) sel$q_star else ld_w_threshold
  no_elbow <- is.null(ld_w_threshold) && sel$no_elbow
  if (verbose && no_elbow)
    warning("RMSC has no elbow: ld_w filtering does not raise discoveries, so ",
            "the local-LD prioritisation may not apply to these data.", call. = FALSE)
  ok <- is.finite(map$ldw0) & is.finite(map$pv0)
  thr_ldw <- stats::quantile(map$ldw0[ok], qstar, na.rm = TRUE)
  map[, cand := ok & ldw0 >= thr_ldw]
  if (!any(map$cand)) stop("No candidate SNPs at the chosen ld_w threshold.")

  ## 2. SNP-level significance within candidates ------------------------------
  ci <- which(map$cand)
  map[, q_snp := NA_real_][ci, q_snp := stats::p.adjust(pv0, "BH")]
  map[, sig := is.finite(q_snp) & q_snp < fdr]
  pcut <- if (any(map$sig)) max(map$pv0[map$sig]) else NA_real_
  if (is.na(pcut))
    warning("No SNP significant within candidates at fdr = ", fdr,
            "; no clusters can be flagged.", call. = FALSE)

  ## 3. distance-restricted single-linkage LD clusters (per chromosome) -------
  ## Per-chromosome r2_link and distance_threshold may be derived from the LD-
  ## decay fit via rho: r2_link = ld_from_rho(b, c, rho_ld) (higher rho_ld = looser
  ## r2), distance_threshold = d_from_rho(a_pred, rho_d) (higher rho_d = wider gap).
  r2_by_chr <- dist_by_chr <- NULL
  if (!is.null(rho_ld) || !is.null(rho_d)) {
    if (is.null(LD_decay) || is.null(LD_decay$decay_sum))
      stop("`rho_ld`/`rho_d` require `LD_decay` (with $decay_sum).")
    ds <- data.table::as.data.table(LD_decay$decay_sum)
    cc <- if ("c" %in% names(ds)) ds$c else rep(1, nrow(ds))
    if (!is.null(rho_ld)) r2_by_chr   <- stats::setNames(ld_from_rho(ds$b, cc, rho_ld), ds$Chr)
    if (!is.null(rho_d))  dist_by_chr <- stats::setNames(d_from_rho(ds$a_pred, rho_d), ds$Chr)
  }
  ## split one single-linkage cluster into physically-contiguous runs
  split_gap <- function(cl, pos, dthr) {
    if (!is.finite(dthr)) return(cl)
    out <- integer(length(cl)); nxt <- 0L
    for (g in unique(cl)) {
      gi <- which(cl == g); o <- order(pos[gi])
      run <- integer(length(gi)); run[o] <- cumsum(c(TRUE, diff(pos[gi][o]) > dthr))
      out[gi] <- nxt + run; nxt <- nxt + max(run)
    }
    out
  }
  map[, cluster := NA_character_]
  for (ch in unique(map$Chr[map$cand])) {
    idx <- which(map$cand & map$Chr == ch)
    mk <- map$marker[idx]; pos <- map$Pos[idx]
    r2c  <- if (!is.null(r2_by_chr) && !is.na(r2_by_chr[ch])) r2_by_chr[[ch]] else r2_link
    dthr <- if (!is.null(dist_by_chr) && !is.na(dist_by_chr[ch])) dist_by_chr[[ch]]
            else if (!is.null(distance_threshold)) distance_threshold else Inf
    if (length(mk) == 1L) { map$cluster[idx] <- paste0(ch, "_c1"); next }
    R <- suppressWarnings(stats::cor(GTs[, mk], use = "pairwise.complete.obs")^2)
    R[!is.finite(R)] <- 0
    cl <- stats::cutree(stats::hclust(stats::as.dist(1 - R), method = "single"), h = 1 - r2c)
    cl <- split_gap(cl, pos, dthr)                    # distance cap (split_by_distance)
    map$cluster[idx] <- paste0(ch, "_c", cl)
  }

  ## 4. cluster table + size-aware null ---------------------------------------
  cd <- map[cand == TRUE]
  cl <- cd[, .(Chr = Chr[1L], n = .N, n_sig = sum(sig),
               start = min(Pos), end = max(Pos), members = list(marker)),
           by = cluster]
  tested <- cl$n_sig >= min_cluster_sig & is.finite(pcut)
  cl[, emp_p := NA_real_]
  if (any(tested)) {
    nd <- switch(null,
      background  = .bg_null(map, cl, pcut, B),
      permutation = {
        if (is.null(Y) || is.null(K)) stop("null = 'permutation' requires Y and K.")
        .perm_null(map, cl, pcut, B, Y, K, GTs, cores)
      })
    for (i in which(tested))
      cl$emp_p[i] <- (1 + sum(nd[[cl$cluster[i]]] >= cl$n_sig[i])) / (B + 1)
  }
  cl[, emp_q := NA_real_][tested, emp_q := stats::p.adjust(emp_p, "BH")]
  cl[, significant := is.finite(emp_q) & emp_q < cluster_fdr]
  data.table::setorder(cl, -significant, emp_p, -n_sig)

  ## diagnostics --------------------------------------------------------------
  chr_sig <- cl[significant == TRUE, unique(Chr)]
  n_chr <- length(unique(map$Chr))
  all_chr <- n_chr > 1L && length(chr_sig) == n_chr
  if (verbose && all_chr)
    warning("Every chromosome contains a significant cluster; genome-wide ",
            "outliers are unusual -- check for residual structure or an ",
            "over-permissive threshold.", call. = FALSE)

  ## structure flag 1: does the association statistic rise with ld_w on EVERY
  ## chromosome?  On a localized-signal dataset the ld_w vs -log10(p) rank
  ## correlation is positive only on the few chromosomes that carry outliers;
  ## when it is positive genome-wide (including nominally neutral chromosomes)
  ## ld_w is tracking residual structure, not signal.
  map[, sig_stat := -log10(pmax(pv0, 1e-300))]
  ldw_cor <- map[is.finite(ldw0) & is.finite(sig_stat),
                 .(cor = suppressWarnings(stats::cor(ldw0, sig_stat,
                        method = "spearman"))),
                 by = Chr]
  frac_chr_pos <- mean(ldw_cor$cor > 0.1, na.rm = TRUE)
  ldw_tracks_structure <- n_chr > 2L && isTRUE(frac_chr_pos > 0.8)
  if (verbose && ldw_tracks_structure)
    warning(sprintf(paste0("ld_w correlates with the association signal on ",
            "%.0f%% of chromosomes -- ld_w is likely tracking genome-wide ",
            "structure rather than localized outliers; calls are unreliable."),
            100 * frac_chr_pos), call. = FALSE)

  ## structure flag 2: do significant clusters span whole chromosomes?  A span
  ## fraction near 1 means a chromosome collapsed into a single cluster, so
  ## per-cluster TP/FP is meaningless (at most one call per chromosome).
  chrlen <- map[, .(len = max(Pos) - min(Pos)), by = Chr]
  sigcl <- cl[significant == TRUE]
  span_frac <- if (nrow(sigcl))
    (sigcl$end - sigcl$start) / chrlen$len[match(sigcl$Chr, chrlen$Chr)] else numeric(0)
  whole_chr_clusters <- length(span_frac) > 0L &&
    stats::median(span_frac, na.rm = TRUE) > 0.5
  if (verbose && whole_chr_clusters)
    warning("Significant clusters span whole chromosomes (median span ",
            "fraction > 0.5); no localizable signal -- structure-dominated.",
            call. = FALSE)

  structure(list(
    clusters = cl[],
    candidates = map[cand == TRUE,
                     .(marker, Chr, Pos, ld_w = ldw0, pval = pv0, q_snp, sig, cluster)],
    ld_w_threshold = qstar, ld_w_value = unname(thr_ldw), pcut = pcut,
    rmsc = data.table::data.table(q = sel$grid, rejections = sel$rejections),
    null = null, B = B,
    diagnostics = list(
      no_elbow = no_elbow, peak_rejections = max(sel$rejections),
      background_rate = if (is.finite(pcut))
        map[cand == FALSE & is.finite(pv0), mean(pv0 <= pcut)] else NA_real_,
      chr_with_signal = chr_sig, all_chr_flagged = all_chr,
      ldw_signal_cor = stats::setNames(ldw_cor$cor, ldw_cor$Chr),
      frac_chr_ldw_positive = frac_chr_pos,
      ldw_tracks_structure = ldw_tracks_structure,
      cluster_span_fraction = span_frac,
      median_span_fraction = if (length(span_frac))
        stats::median(span_frac, na.rm = TRUE) else NA_real_,
      whole_chr_clusters = whole_chr_clusters)
  ), class = "ld_outlier_clusters")
}

## below-threshold background null: draw each cluster's null n_sig from the
## low-ld_w SNPs, size-matched and rotated to preserve background autocorrelation
.bg_null <- function(map, cl, pcut, B) {
  out <- stats::setNames(vector("list", nrow(cl)), cl$cluster)
  for (ch in unique(cl$Chr)) {
    sbg <- as.integer(map[cand == FALSE & Chr == ch & is.finite(pv0)][order(Pos), pv0] <= pcut)
    nb <- length(sbg); if (nb < 2L) next
    for (cc in cl[Chr == ch, cluster]) {
      ni <- cl[cluster == cc, n]
      out[[cc]] <- vapply(seq_len(B), function(b) {
        st <- sample.int(nb, 1L); sum(sbg[((st - 1L + seq_len(ni) - 1L) %% nb) + 1L])
      }, numeric(1))
    }
  }
  out
}

## EMMAX permutation null, fast path. K is fixed across permutations, so each
## chromosome's kinship is eigendecomposed and the candidate genotypes rotated
## ONCE; every permutation then only re-estimates the variance component (REML
## with a pre-computed eig.R) and does an O(nm) whitened single-SNP F-test. The
## F is numerically identical to emmax() -- RSS is invariant to the choice of
## whitening square-root -- so it is comparable to the observed emmax() p-values.
.perm_null <- function(map, cl, pcut, B, Y, K, GTs, cores) {
  n <- length(Y); df1 <- 1L; df2 <- n - 2L
  chrs <- unique(cl$Chr)
  one_n <- matrix(1, n, n) / n
  pre <- lapply(chrs, function(ch) {
    mk <- map[cand == TRUE & Chr == ch, marker]
    Kc <- if (is.list(K)) K[[ch]] else K
    K_norm <- (n - 1) / sum((diag(n) - one_n) * Kc) * Kc
    Xo <- matrix(1, n, 1L)
    ev <- eigen(K_norm, symmetric = TRUE)                    # for whitening
    list(mk = mk, clus = map[match(mk, marker), cluster],
         K_norm = K_norm, Xo = Xo,
         eigR = emma.eigen.R.wo.Z(K_norm, Xo),               # for REML
         V = ev$vectors, lam = ev$values,
         Xt = crossprod(ev$vectors, GTs[, mk, drop = FALSE]),
         xot = as.numeric(crossprod(ev$vectors, Xo)))
  })
  one <- function(b) {
    yb <- sample(Y); counts <- stats::setNames(integer(nrow(cl)), cl$cluster)
    for (p in pre) {
      re <- emma.REMLE(yb, p$Xo, p$K_norm, eig.R = p$eigR)
      wv <- 1 / sqrt(re$vg * p$lam + re$ve)                  # whitening weights
      yt <- as.numeric(crossprod(p$V, yb)) * wv
      xo <- p$xot * wv; Xt <- p$Xt * wv
      a2 <- sum(xo * xo); ra <- yt - xo * (sum(xo * yt) / a2)
      RSS0 <- sum(ra * ra)
      Rb   <- Xt - outer(xo, as.numeric(crossprod(xo, Xt)) / a2)  # residualise on intercept
      RSSf <- RSS0 - as.numeric(crossprod(ra, Rb))^2 / colSums(Rb * Rb)
      Fj   <- (RSS0 / RSSf - 1) * df2 / df1
      pj   <- stats::pf(Fj, df1, df2, lower.tail = FALSE)
      hit  <- p$mk[is.finite(pj) & pj <= pcut]
      if (length(hit)) { tt <- table(p$clus[match(hit, p$mk)]); counts[names(tt)] <- counts[names(tt)] + as.integer(tt) }
    }
    counts
  }
  reps <- if (cores > 1L) parallel::mclapply(seq_len(B), one, mc.cores = cores) else lapply(seq_len(B), one)
  mat <- do.call(rbind, reps)
  stats::setNames(lapply(cl$cluster, function(cc) mat[, cc]), cl$cluster)
}

#' @export
print.ld_outlier_clusters <- function(x, ...) {
  cl <- x$clusters
  cat(sprintf("<ld_outlier_clusters>  ld_w q* = %.2f (value %.3f) | null = %s | B = %d\n",
              x$ld_w_threshold, x$ld_w_value, x$null, x$B))
  cat(sprintf("  %d candidates in %d LD clusters | SNP p-cutoff = %.3g | background rate = %.3f\n",
              nrow(x$candidates), nrow(cl), x$pcut,
              x$diagnostics$background_rate))
  ns <- sum(cl$significant, na.rm = TRUE)
  cat(sprintf("  significant clusters (emp_q < cutoff): %d%s\n", ns,
              if (x$diagnostics$all_chr_flagged) "  [every chromosome flagged!]" else ""))
  if (ns) print(cl[significant == TRUE,
                   .(cluster, Chr, n, n_sig, start, end,
                     emp_p = round(emp_p, 4), emp_q = round(emp_q, 4))])
  invisible(x)
}
