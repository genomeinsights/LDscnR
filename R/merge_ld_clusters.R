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
#' @param cores Number of chromosomes to process in parallel.
#'
#' @return A list with the same shape as [ld_complexity_reduction()]:
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
#' @export
merge_ld_clusters <- function(GTs, ld_result, LD_decay, rho = 0.5, cores = 1) {

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

    if (nrow(chr_clusters) < 2) {
      chr_clusters[, density := vapply(members, raw_density, numeric(1))]
      return(chr_clusters)
    }

    ## one consensus genotype per cluster (rows = individuals, cols = clusters)
    consensus <- vapply(chr_clusters$members, function(mk) {
      expected_gt_dosage(GTs[, mk, drop = FALSE])
    }, numeric(nrow(GTs)))
    colnames(consensus) <- chr_clusters$CL_id

    R2 <- suppressWarnings(stats::cor(consensus, use = "pairwise.complete.obs")^2)
    R2[!is.finite(R2)] <- 0  ## e.g. a constant consensus signal -- treat as unlinked
    diag(R2) <- 1

    hc  <- stats::hclust(stats::as.dist(1 - R2), method = "complete")
    grp <- stats::cutree(hc, h = 1 - r2_th)

    data.table::rbindlist(lapply(split(chr_clusters$CL_id, grp), function(ids) {
      sub <- chr_clusters[CL_id %in% ids]
      if (nrow(sub) == 1) {
        sub[, density := raw_density(members[[1]])]
        return(sub)
      }

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
    }), use.names = TRUE)  ## column ORDER differs between the two branches
    ## above (density lands in a different position); rbindlist()'s default
    ## use.names="check" only warns on that and still binds positionally,
    ## which would silently put density values in the wrong rows.
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

  list(
    map_snp  = map_snp,
    clusters = new_clusters,
    pruned   = new_clusters$core_snp
  )
}
