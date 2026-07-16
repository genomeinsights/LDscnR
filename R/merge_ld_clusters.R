#' Merge Fragmented LD Clusters via Cluster-Level Correlation
#'
#' Second-stage clustering on top of [ld_complexity_reduction()]'s output.
#' As documented there ("Known limitation"), a genuinely large,
#' low-recombination block (a centromere, an inversion) can get fragmented
#' into several separate clusters purely because pairwise LD beyond the
#' `LD_decay` edge list's sliding window was never evaluated -- not because
#' the block isn't really one unit. This reunites such fragments.
#'
#' The trick is doing the merge at the CLUSTER level, not the SNP level:
#' [ld_complexity_reduction()] already reduces the marker set to a much
#' smaller number of clusters, so an unrestricted, fully all-pairwise
#' comparison across every cluster on a chromosome -- something far too
#' costly for the raw markers -- becomes cheap. Each cluster is first
#' reduced to one consensus genotype ([expected_gt_dosage()] on its
#' polarized member genotypes), those consensus signals are correlated
#' against every other cluster on the same chromosome with no distance
#' restriction at all, and complete linkage at the same rho-derived
#' threshold [ld_complexity_reduction()] uses merges clusters whose
#' consensus signals are still highly correlated despite having been split
#' apart. Comparisons are restricted to within the same chromosome --
#' there's no physical basis for merging clusters from different
#' chromosomes, since the fragmentation this fixes is specifically a
#' sliding-window artifact along one chromosome's marker order.
#'
#' @param GTs Numeric matrix of genotype dosages, individuals x markers,
#'   with column names matching the markers in `ld_result`.
#' @param ld_result Output of [ld_complexity_reduction()].
#' @param LD_decay The same `ld_decay` object passed to
#'   [ld_complexity_reduction()], used for the per-chromosome threshold.
#' @param rho Numeric in (0, 1), same meaning as in
#'   [ld_complexity_reduction()]. Default `0.5`.
#' @param cores Number of cores. Used both to process chromosomes in
#'   parallel and, within each chromosome, to compute cluster consensus
#'   genotypes ([expected_gt_dosage()]) in parallel -- the latter can be a
#'   real bottleneck when a chromosome has many clusters, since it's one
#'   function call per cluster. If `cores > 1` and multiple chromosomes are
#'   processed at once, the two loops nest (each chromosome-level worker
#'   spawns its own consensus-level workers), so total processes can reach
#'   `cores^2`; this is only a concern for genome-wide, many-chromosome,
#'   high-`cores` runs -- the common case here is one chromosome (or
#'   `cores = 1`), where there's no nesting.
#' @param min_n_snps Integer. Only clusters with `n_snps >= min_n_snps` enter
#'   the all-pairwise cluster comparison; smaller clusters pass straight
#'   through unmerged (`density` is still computed for them). The dominant
#'   cost here is the cross-cluster correlation matrix and `hclust()`, both
#'   roughly quadratic in the NUMBER of clusters compared, not their size --
#'   so excluding a chromosome's small clusters (usually the large majority)
#'   can be a large speedup. This is a heuristic, not a guarantee: a
#'   genuinely large true block could in principle fragment into many small
#'   pieces rather than a few large ones, in which case a small piece could
#'   still need merging and would be missed here. Default `1` (no
#'   filtering, i.e. identical to previous behaviour) since whether that
#'   trade-off is worth it depends on the data -- e.g. check
#'   `cor(ld_result$map_snp$ld_w_09x, ld_result$map_snp$n_loci)` first:
#'   a real positive correlation (fragmentation-prone/high-LD regions
#'   already coming out of Stage 1 with larger clusters) is what makes a
#'   higher threshold safe to use.
#'
#' @return An object of class `"merge_ld_clusters"` (a list with a
#'   [print.merge_ld_clusters()] method -- printed automatically once at the
#'   end of the call, and again any time the returned object is printed
#'   later), with the same shape as [ld_complexity_reduction()]'s output:
#'   `map_snp`, `clusters`, `pruned`. Clusters that get merged share a new,
#'   combined `members` list and a representative chosen the same way as
#'   the first stage (highest median r2 -- here, to the other merged
#'   clusters' consensus signals rather than to raw SNPs); clusters that
#'   don't merge with anything pass through with their original `CL_id`,
#'   `core_snp`, and `median_ld` unchanged. `clusters` also gains a
#'   `density` column (for every cluster, not just merged ones): the
#'   fraction of ALL pairwise raw-marker comparisons among a cluster's
#'   members -- computed directly from `GTs`, with no distance limit -- that
#'   clear `r2_th`. `median_ld` says how well one representative correlates
#'   with the rest; `density` says how uniformly correlated the whole
#'   cluster actually is, which is the more direct read on whether a merge
#'   reunited a genuine block (high density) or just chained loosely
#'   related clusters together (low density) -- e.g. a real large block on
#'   Chr26 checked this way came out at ~0.99, versus 0.19-0.32 for
#'   demonstrably chained components found earlier by the same check.
#'   `NA` for singleton (one-marker) clusters, where there's no pair to
#'   compare.
#'
#' @examples
#' \dontrun{
#' ## Recommended workflow when only some loci need the merge step -- e.g.
#' ## flagging by ld_w (see compute_ld_w()): loci with unusually high local
#' ## LD support are where the sliding window behind LD_decay's edge list is
#' ## most likely to have fragmented one true block into several Stage-1
#' ## clusters. Run Stage 1 on the WHOLE marker set together FIRST (cheap
#' ## regardless of scale, and it uses the real window-covered edges to
#' ## group markers correctly), then flag CLUSTERS -- not markers -- and
#' ## only pay merge_ld_clusters()'s all-pairwise cost on the flagged ones.
#' ##
#' ## Pre-splitting markers into "high"/"low" ld_w BEFORE clustering, then
#' ## running Stage 1 separately on each half, looks equivalent but isn't:
#' ## it can sever a real block right at the threshold boundary, since two
#' ## adjacent markers can straddle the cutoff even though `el` connects
#' ## them directly. Checked on one real chromosome: 226 of 15,524 Stage-1
#' ## clusters from a combined run mixed low- and high-ld_w members, and a
#' ## pre-split would have cut 975 low-ld_w markers loose from their real
#' ## cluster. Flagging clusters (post hoc, using Stage 1's own boundaries)
#' ## avoids that entirely.
#' stage1 <- ld_complexity_reduction(map = map, LD_decay = ld_decay, rho = 0.5)
#'
#' ## classify via a separate vector, not a mutated column on $clusters --
#' ## merge_ld_clusters() expects ld_result$clusters in exactly the shape
#' ## ld_complexity_reduction() returns it; a stray extra column carried
#' ## through wouldn't survive every internal branch identically (merged
#' ## clusters get freshly built rows without it, passed-through ones keep
#' ## it), producing inconsistent columns and a rbindlist() error
#' ld_w_threshold <- 0.2
#' marker_ldw <- setNames(map$ld_w_095, map$marker)
#' needs_merge <- vapply(
#'   stage1$clusters$members, function(mk) any(marker_ldw[mk] > ld_w_threshold), logical(1)
#' )
#'
#' flagged   <- stage1$clusters[needs_merge]
#' unflagged <- stage1$clusters[!needs_merge]
#'
#' ## a ld_result-shaped list restricted to the flagged clusters -- the only
#' ## piece merge_ld_clusters() actually reads is $clusters; $map_snp is
#' ## trimmed to match purely so the two stay consistent, and $pruned isn't
#' ## used internally at all
#' ld_result_flagged <- list(
#'   map_snp  = stage1$map_snp[marker %in% unlist(flagged$members)],
#'   clusters = flagged,
#'   pruned   = flagged$core_snp
#' )
#'
#' merged <- merge_ld_clusters(
#'   GTs = GTs[, unlist(flagged$members)],
#'   ld_result = ld_result_flagged, LD_decay = ld_decay, rho = 0.5
#' )
#'
#' ## unflagged clusters' representatives pass straight through untouched
#' pruned_markers <- c(unflagged$core_snp, merged$pruned)
#' }
#'
#' @export
merge_ld_clusters <- function(GTs, ld_result, LD_decay, rho = 0.5, cores = 1, min_n_snps = 1) {

  clusters <- data.table::copy(ld_result$clusters)
  chr_levels <- unique(clusters$Chr)

  merged_by_chr <- parallel_apply(chr_levels, function(ch) {

    chr_clusters <- clusters[Chr == ch]

    ds <- LD_decay$decay_sum[Chr == ch]
    if (nrow(ds) != 1) {
      stop("Expected exactly one decay_sum row for chr '", ch, "', found ", nrow(ds),
           ". Check LD_decay$decay_sum$Chr labels.")
    }
    r2_th <- ld_from_rho(b = ds$b, c = ds$c, rho = rho)

    raw_density <- function(mk) {
      if (length(mk) < 2) return(NA_real_)
      r2 <- suppressWarnings(stats::cor(GTs[, mk, drop = FALSE], use = "pairwise.complete.obs")^2)
      r2[!is.finite(r2)] <- 0
      mean(r2[upper.tri(r2)] >= r2_th)
    }

    ## clusters below min_n_snps skip the expensive all-pairs comparison
    ## entirely and pass through unmerged -- see @param min_n_snps for the
    ## reasoning. Appended back after everything below. vapply() over an
    ## empty list correctly returns numeric(0), so this is safe (and still
    ## creates the `density` column) even when below_threshold has 0 rows --
    ## the default min_n_snps=1 means it always does, since n_snps is never
    ## < 1; skipping the assignment there entirely used to leave `density`
    ## missing and break the final rbindlist() downstream.
    below_threshold <- chr_clusters[n_snps < min_n_snps]
    chr_clusters     <- chr_clusters[n_snps >= min_n_snps]
    below_threshold[, density := vapply(members, raw_density, numeric(1))]

    if (nrow(chr_clusters) < 2) {
      chr_clusters[, density := vapply(members, raw_density, numeric(1))]
      return(data.table::rbindlist(list(chr_clusters, below_threshold), use.names = TRUE))
    }

    ## one consensus genotype per cluster (rows = individuals, cols = clusters).
    ## Singletons skip expected_gt_dosage() entirely: polarize_genotypes()
    ## already short-circuits on a single column, so for n_snps == 1 the
    ## function is provably the identity on that column (rowMeans(na.rm=TRUE)
    ## over one column returns that column, NaN->NA conversion included) --
    ## just cheaper to grab the column directly than pay ~12k redundant calls.
    ## The non-singleton remainder is parallelized with the same `cores` used
    ## for the outer per-chromosome split -- fine in practice since this
    ## function is almost always called with a single chromosome's worth of
    ## clusters at a time (or cores = 1 for genome-wide runs), so the two
    ## loops don't actually run nested.
    is_single <- chr_clusters$n_snps == 1L
    consensus <- matrix(NA_real_, nrow = nrow(GTs), ncol = nrow(chr_clusters))
    if (any(is_single)) {
      single_markers <- vapply(chr_clusters$members[is_single], `[[`, character(1), 1)
      consensus[, is_single] <- GTs[, single_markers, drop = FALSE]
    }
    if (any(!is_single)) {
      consensus[, !is_single] <- do.call(cbind, parallel_apply(chr_clusters$members[!is_single], function(mk) {
        expected_gt_dosage(GTs[, mk, drop = FALSE])
      }, cores = cores))
    }
    colnames(consensus) <- chr_clusters$CL_id

    R2 <- suppressWarnings(stats::cor(consensus, use = "pairwise.complete.obs")^2)
    R2[!is.finite(R2)] <- 0  ## e.g. a constant consensus signal -- treat as unlinked
    diag(R2) <- 1

    hc  <- stats::hclust(stats::as.dist(1 - R2), method = "complete")
    grp <- stats::cutree(hc, h = 1 - r2_th)

    groups <- split(chr_clusters$CL_id, grp)
    n_grp  <- length(groups)

    message(
      ch, ": ", nrow(chr_clusters), " clusters >= ", min_n_snps, " SNPs -> ", n_grp, " groups",
      if (nrow(below_threshold)) paste0(" (+", nrow(below_threshold), " smaller, passed through)") else ""
    )
    pb <- utils::txtProgressBar(min = 0, max = n_grp, style = 3)
    on.exit(close(pb))

    out <- vector("list", n_grp)

    for (i in seq_len(n_grp)) {
      ids <- groups[[i]]
      sub <- chr_clusters[CL_id %in% ids]

      out[[i]] <- if (nrow(sub) == 1) {
        sub[, density := raw_density(members[[1]])]
        sub
      } else {
        ## representative: the original cluster whose consensus has the
        ## highest median r2 to the OTHER clusters being merged with it --
        ## same "highest median r2" principle Stage 1 uses, just one level up.
        ## R2's dimnames are character (colnames() always coerces), but CL_id
        ## itself is an integer (ld_complexity_reduction()'s .GRP, assigned
        ## globally across chromosomes) -- indexing a matrix with a bare
        ## integer vector is POSITIONAL, not a name lookup, so this must be
        ## character or it silently grabs the wrong (or out-of-bounds) cells.
        sub_r2 <- R2[as.character(sub$CL_id), as.character(sub$CL_id), drop = FALSE]
        diag(sub_r2) <- NA
        med_r2 <- apply(sub_r2, 1, stats::median, na.rm = TRUE)
        best <- sub$CL_id[which.max(med_r2)]

        merged_members <- unlist(sub$members, use.names = FALSE)

        data.table::data.table(
          Chr       = ch,
          CL_id     = best,  ## reuse an original CL_id so the merge stays traceable
          core_snp  = sub[CL_id == best, core_snp],
          median_ld = max(med_r2),
          density   = raw_density(merged_members),
          n_snps    = sum(sub$n_snps),
          members   = list(merged_members)
        )
      }

      utils::setTxtProgressBar(pb, i)
    }

    ## column ORDER differs between branches (density lands in a different
    ## position); rbindlist()'s default use.names="check" only warns on that
    ## and still binds positionally, which would silently put density values
    ## in the wrong rows.
    data.table::rbindlist(c(out, list(below_threshold)), use.names = TRUE)
  }, cores = cores)

  new_clusters <- data.table::rbindlist(merged_by_chr, use.names = TRUE)

  member_dt <- new_clusters[, .(marker = unlist(members)),
                             by = .(CL_id, n_loci = n_snps, core = core_snp, median_ld)]
  member_dt[, is_core := marker == core]
  member_dt[, core := NULL]

  ## same X[i, on=] pattern ld_complexity_reduction() uses -- preserves
  ## map_snp's row order and every other original column, only the cluster
  ## assignment columns are replaced
  old_map_snp <- data.table::copy(ld_result$map_snp)
  old_map_snp[, c("CL_id", "n_loci", "is_core", "median_ld") := NULL]

  map_snp <- member_dt[, .(marker, CL_id, n_loci, is_core, median_ld)][
    old_map_snp, on = "marker"
  ]

  out <- list(
    map_snp  = map_snp,
    clusters = new_clusters,
    pruned   = new_clusters$core_snp
  )
  attr(out, "n_input_clusters") <- nrow(ld_result$clusters)
  class(out) <- "merge_ld_clusters"
  print(out)
  out
}

#' Print a merge_ld_clusters Object
#'
#' Summarizes the result of [merge_ld_clusters()]: how many clusters were
#' merged, the resulting cluster-size distribution, representative markers'
#' median r2, and each cluster's raw-marker-pair density (how much of a
#' merge reunited a genuine block versus chained loosely related clusters).
#'
#' @param x A `"merge_ld_clusters"` object.
#' @param digits Number of significant digits for printed statistics.
#' @param ... Unused; present for S3 consistency.
#'
#' @return `x`, invisibly.
#'
#' @export
print.merge_ld_clusters <- function(x, digits = 3, ...) {

  cat("<merge_ld_clusters>\n")

  n_before <- attr(x, "n_input_clusters")
  n_after  <- nrow(x$clusters)

  cat("\nClusters:\n")
  if (!is.null(n_before)) {
    cat(
      "  input clusters:", format(n_before, big.mark = ","),
      " -> merged clusters:", format(n_after, big.mark = ","),
      " (", signif(100 * (1 - n_after / n_before), digits), "% reduction)\n",
      sep = ""
    )
  } else {
    cat("  merged clusters:", format(n_after, big.mark = ","), "\n")
  }
  cat("  chromosomes:", data.table::uniqueN(x$clusters$Chr), "\n")

  n_snps <- x$clusters$n_snps

  cat("\nCluster sizes:\n")
  cat(
    "  median =", stats::median(n_snps),
    " range = [", min(n_snps), ", ", max(n_snps), "]\n",
    sep = ""
  )
  cat(
    "  singleton clusters:", sum(n_snps == 1L),
    " (", signif(100 * mean(n_snps == 1L), digits), "%)\n",
    sep = ""
  )

  ld <- x$clusters$median_ld
  ld <- ld[!is.na(ld)]
  if (length(ld)) {
    cat("\nRepresentative median r2 (non-singleton clusters):\n")
    cat(
      "  median =", signif(stats::median(ld), digits),
      " range = [", signif(min(ld), digits), ", ", signif(max(ld), digits), "]\n",
      sep = ""
    )
  }

  dens <- x$clusters$density
  dens <- dens[!is.na(dens)]
  if (length(dens)) {
    cat("\nRaw-marker-pair density (non-singleton clusters):\n")
    cat(
      "  median =", signif(stats::median(dens), digits),
      " range = [", signif(min(dens), digits), ", ", signif(max(dens), digits), "]\n",
      sep = ""
    )
    cat("  (fraction of all pairwise raw-marker comparisons that clear r2_th;\n")
    cat("   low values on a merged cluster mean it was chained, not reunited)\n")
  }

  cat("\nStored components:\n")
  cat("  map_snp, clusters, pruned\n")

  invisible(x)
}
