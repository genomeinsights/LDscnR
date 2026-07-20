# Combined LD-pruning + eMLG generation, replacing merge_ld_clusters() for
# clusters flagged by high local LD support (ld_w).
#
# Context: a common pattern is to run ld_complexity_reduction() once
# genome-wide, flag clusters containing a marker above an ld_w threshold,
# then run merge_ld_clusters() (complete linkage) on just the flagged
# clusters to reunite sliding-window-fragmented blocks -- and SEPARATELY
# build eMLGs from that output and dynamic-cut them AGAIN (average linkage,
# score_eMLG-gated) for Ohta-statistic/long-range-LD work. That's two
# independent all-pairs correlation passes over overlapping/related data.
#
# This does ONE pass instead: Stage 1 -> flag -> distance-restricted
# dynamic cut directly on the flagged Stage-1 clusters, producing BOTH a
# pruned representative marker set (for LD-pruning/GRM) AND an eMLG matrix
# (for downstream long-range LD analysis) from the same clustering.
#
# Why this is defensible for LD-pruning too, not just eMLG summarization:
#   - complete linkage's strict "every pairwise r2 must clear the
#     threshold" guarantee exists to avoid pruning away real, independent
#     genetic signal -- appropriate when markers COULD be genuinely
#     unlinked. But in a young population's low-recombination region, that
#     concern doesn't really apply: with very few actual recombination
#     events having occurred there, the region genuinely behaves as a
#     small number of effectively-independent haplotype blocks, not many
#     independent markers being falsely conflated. Representing it with
#     fewer units isn't discarding signal, it's accurately reflecting the
#     true effective degrees of freedom.
#   - the real risk was merging PHYSICALLY DISTANT clusters whose
#     correlation might reflect intra-chromosomal epistatic selection
#     rather than physical linkage -- exactly the signal later analysis
#     (Ohta statistics) wants to detect, not erase by premature merging.
#     Fixed by restricting merging to physically-contiguous runs
#     (single-linkage BY DISTANCE: consecutive cluster gap <=
#     distance_threshold). Note this deliberately does NOT cap a run's
#     total physical span -- a genuine non-recombining block is expected
#     to be one long, closely-spaced, contiguous run that can legitimately
#     span megabases; what should break it into separate runs is a real
#     gap between neighbours, not the block's overall size.
#   - distance_threshold itself defaults to a per-chromosome value derived
#     from rho (d_from_rho(a_pred, rho), reusing the same rho already
#     implied by ld_w_col's naming, e.g. "ld_w_095") rather than one fixed
#     bp constant -- and specifically from a_pred (chromosome-size-
#     corrected predicted decay rate), not each chromosome's own observed
#     a. Checked directly on real data: a chromosome with an ~100x slower
#     observed decay rate than typical (plausibly a large non-recombining
#     supergene-type region) and a typical chromosome of similar size get
#     a near-identical, well-calibrated ~11kb window at rho=0.95 this way,
#     capturing ~98-99% of real local gaps between flagged clusters while
#     still treating the rare large gaps (the top 1-2%, up to several
#     hundred kb) as presumptive breaks. Using the chromosome's own
#     observed a instead would have given the anomalous chromosome a
#     ~1.5Mb window -- over 2x its largest actually-observed gap, i.e.
#     license to bridge nearly the whole chromosome into one run
#     regardless of real local structure. That's exactly backwards: a
#     chromosome's own anomalous decay is often the reason its internal
#     structure is worth NOT smoothing over by construction, not a
#     justification for a more permissive merge threshold.
#
# See dynamic_cut_eMLG.R for the average-vs-single-vs-complete-linkage
# comparison and the score_eMLG/pair_r2 gating rationale/bugs already
# fixed; this file reuses that function unchanged as the within-run
# merging step and adds: flagging, distance-restricted run splitting
# across chromosomes, unflagged cluster passthrough, and representative
# marker selection.
#
# Cost: the dominant cost (all-pairs cor()) scales quadratically with
# cluster count (confirmed directly: ~0.01s at 292 clusters, ~31s at
# 15,000) -- restricting to ld_w-flagged clusters is load-bearing for
# keeping this tractable, not a workaround to relax. At full-chromosome,
# unfiltered scale this would hit the same wall merge_ld_clusters() did
# before min_n_snps filtering.
#
# eMLG computation is itself opt-out/opt-in at two levels, since callers
# doing LD-pruning only don't need eMLGs at all, and long-range-LD analysis
# will apply its own n_loci threshold downstream anyway (too few loci isn't
# informative for detecting extended LD patterns regardless):
#   - compute_unflagged_eMLG: skip eMLG computation for unflagged clusters
#     entirely (usually the large majority of clusters) when only the
#     pruned marker set is needed.
#   - min_n_loci_eMLG: skip storing an eMLG for any group (flagged or
#     unflagged) below this many raw loci -- the pruned representative
#     marker is unaffected either way, this only trims the returned eMLG
#     matrix.

## expected_gt_dosage() is provably the identity for a single-marker
## cluster (polarize_genotypes() short-circuits at ncol<=1, rowMeans() over
## one column returns that column unchanged) -- skip the redundant
## function-call/polarization overhead and use the raw column directly.
## Always correct, no output change, so applied unconditionally rather than
## gated by a parameter.
consensus_dosage <- function(GTs, mk) {
  if (length(mk) == 1) return(GTs[, mk])
  expected_gt_dosage(GTs[, mk, drop = FALSE])
}

## Splits clusters into physically-contiguous runs: single linkage BY
## DISTANCE, i.e. a new run starts whenever the gap to the previous
## (position-sorted) cluster exceeds distance_threshold. A run's TOTAL span
## is deliberately unbounded -- see file header.
split_by_distance <- function(pos_min, pos_max, distance_threshold = 5e5) {
  ord <- order(pos_min)
  gap <- pos_min[ord] - c(NA, pos_max[ord][-length(ord)])
  new_run <- is.na(gap) | gap > distance_threshold
  run <- cumsum(new_run)
  run_id <- integer(length(ord))
  run_id[ord] <- run
  names(run_id) <- names(pos_min)[ord]
  run_id[names(pos_min)]
}

## Representative marker for a final group: among its constituent Stage-1
## clusters, the one whose consensus signal has the highest median r2 to
## the OTHERS in the group (same "highest median r2" principle
## ld_complexity_reduction()/merge_ld_clusters() use), using that cluster's
## own core_snp. Cheap even though computed per-group: groups are small
## (tens of constituent clusters at most, not raw markers).
pick_representative <- function(cl_ids, eMLG, stage1_clusters) {
  if (length(cl_ids) == 1) {
    return(stage1_clusters[CL_id == cl_ids, core_snp])
  }
  sub_r2 <- suppressWarnings(stats::cor(eMLG[, cl_ids, drop = FALSE], use = "pairwise.complete.obs")^2)
  sub_r2[!is.finite(sub_r2)] <- 0
  diag(sub_r2) <- NA
  med_r2 <- apply(sub_r2, 1, stats::median, na.rm = TRUE)
  best_cl <- cl_ids[which.max(med_r2)]
  stage1_clusters[CL_id == best_cl, core_snp]
}

#' Combined LD-pruning + eMLG generation via distance-restricted dynamic cut
#'
#' Second-stage clustering on top of [ld_complexity_reduction()]'s output,
#' combining what a prior two-pass approach (complete-linkage re-merging of
#' flagged clusters, then a separate [dynamic_cut_eMLG()] pass) used to do
#' as two separate passes into one: Stage-1 clusters flagged by
#' local LD support (`ld_w`) are consolidated via a single distance-
#' restricted, quality-gated average-linkage dynamic cut (see
#' `dynamic_cut_eMLG.R`'s file header for the linkage comparison and the two
#' historical bugs it's gated against), producing both a pruned
#' representative marker set and an eMLG matrix from the same clustering.
#' Unflagged clusters pass straight through unchanged. See this file's
#' header comment for the full rationale, including why this is defensible
#' for LD-pruning itself (not just eMLG summarization) in a young,
#' low-recombination population, and why merging is restricted to
#' physically-contiguous runs rather than allowed genome-wide.
#'
#' @param GTs Numeric matrix, individuals x markers.
#' @param stage1 Output of [ld_complexity_reduction()].
#' @param ld_w_col Name of the ld_w column in `stage1$map_snp` used to flag
#'   which clusters need the expensive treatment (a cluster is flagged if
#'   ANY member exceeds `ld_w_threshold`).
#' @param ld_w_threshold Flagging threshold.
#' @param LD_decay Object from [compute_LD_decay()], used to derive a
#'   per-chromosome `distance_threshold` from `rho` when `distance_threshold`
#'   is not supplied directly (via `d_from_rho(a_pred, rho)`, using each
#'   chromosome's size-corrected predicted decay rate `a_pred`, not its own
#'   observed `a` -- see `distance_threshold` below for why). Only required
#'   when `distance_threshold` is `NULL` (the default).
#' @param rho Relative LD threshold used to derive `distance_threshold` when
#'   it isn't supplied directly. Reuses the same `rho` already implied by
#'   `ld_w_col`'s naming convention (e.g. `"ld_w_095"`) by default, so one
#'   `rho` governs both what counts as "locally elevated" (upstream, via
#'   [compute_ld_w()]) and what counts as "still physically contiguous"
#'   (here). Default `0.95`.
#' @param score_threshold Minimum `score_eMLG` a candidate merge's result
#'   must clear (see [dynamic_cut_eMLG()]).
#' @param min_r2 Minimum r2 required directly between the two sides being
#'   merged (see [dynamic_cut_eMLG()]).
#' @param distance_threshold Max consecutive-gap in bp allowed within one
#'   physically-contiguous run (see `split_by_distance()`). Clusters more
#'   than this apart are never merged, regardless of correlation. Default
#'   `NULL`: derived per chromosome from `rho` and `LD_decay` instead of a
#'   single fixed value (via `d_from_rho(a_pred, rho)`) -- deliberately
#'   using `a_pred` (chromosome-size-corrected decay rate), not each
#'   chromosome's own observed `a`. Checked directly on real data: an
#'   anomalous, ~100x-slower-decay chromosome and a typical one (similar
#'   size) get a near-identical, well-calibrated ~11kb window at
#'   `rho = 0.95` this way, capturing ~98-99% of real local gaps while
#'   still treating the rare large gaps as presumptive breaks. Using the
#'   observed `a` instead would have given the anomalous chromosome a
#'   ~1.5Mb window -- more than double the largest gap actually observed on
#'   it, i.e. license to bridge essentially the whole chromosome into one
#'   run regardless of its real local structure. Supplying a single numeric
#'   value here still works as before (e.g. for hand-built test fixtures
#'   with no matching decay object) and skips the `rho`-based derivation
#'   entirely. Ignored entirely if `genetic_map` is supplied.
#' @param genetic_map If supplied (together with `cM_threshold`), run
#'   splitting is done in genetic (cM) distance instead of physical (bp)
#'   distance -- `distance_threshold`/`rho`/`LD_decay` are not consulted at
#'   all in that case. A data.table/data.frame with `Chr`, `Pos`, `cM`
#'   columns (see [interpolate_cM()]; any chromosome-naming differences
#'   from `stage1`/`map_snp` must be reconciled by the caller first).
#'   Physical distance is a poor, sometimes wildly misleading proxy for
#'   recombination distance: checked directly on real data, some flagged
#'   cluster pairs over 1 Mb apart physically sat at cM distance 0 (fully
#'   linked -- no measurable recombination between them at all), while
#'   others only a few hundred kb apart already spanned several cM. Use
#'   this whenever a genetic map for the study system is available.
#'   Default `NULL` (bp-based `distance_threshold` behaviour).
#' @param cM_threshold Max consecutive-gap in cM allowed within one
#'   genetically-contiguous run, analogous to `distance_threshold` but in
#'   cM units -- required if `genetic_map` is supplied. Unlike
#'   `distance_threshold`'s per-chromosome bp derivation, a single
#'   `cM_threshold` is used for every chromosome: recombination distance is
#'   already the thing that varies by chromosome/region, so there is
#'   nothing left for a physical unit to further correct for once
#'   distances are measured in cM.
#' @param compute_unflagged_eMLG If `FALSE`, skip eMLG computation for
#'   unflagged clusters entirely (usually the large majority) -- for
#'   callers who only need the pruned marker set, not eMLGs. Default `TRUE`.
#' @param min_n_loci_eMLG Skip storing an eMLG for any group (flagged or
#'   unflagged) with fewer than this many raw loci -- the pruned
#'   representative marker is returned regardless, this only trims the
#'   returned `eMLG` matrix. Default `1` (no filtering).
#' @param min_n_loci_flag Additionally flag (send into the distance-
#'   restricted dynamic cut, alongside the `ld_w_threshold` criterion) any
#'   Stage-1 cluster with at least this many raw loci, even if its ld_w
#'   never exceeds `ld_w_threshold` -- giving already-substantial low-ld_w
#'   clusters a chance to merge with a physically-nearby, correlated
#'   neighbour instead of just passing through unchanged. Default `Inf`
#'   (disabled -- flagging is by `ld_w_threshold` alone, existing
#'   behaviour). NOTE: this only usefully restricts the flagged set when
#'   `ld_w_threshold` is a real, discriminating cutoff -- at
#'   `ld_w_threshold = 0`, virtually every cluster may already exceed it
#'   (checked on one real dataset: 99.998% of markers had a positive local
#'   LD-support value), so `min_n_loci_flag` can end up adding nothing on
#'   top of that; the intended use is a real threshold (e.g. 0.025) plus a
#'   `min_n_loci_flag` that pulls in the small number of additional
#'   substantial-but-low-ld_w clusters on top of it.
#' @param cores Number of cores for the unflagged-cluster eMLG loop (the
#'   large, embarrassingly-parallel one -- flagged clusters and the
#'   distance-restricted dynamic cut are usually far fewer and stay
#'   serial). Uses `parallel::mclapply()` (fork-based copy-on-write) when
#'   `cores > 1` and not on Windows, so `GTs` is NOT duplicated per worker
#'   as long as workers only read it (they do) -- no need to pre-split GTs
#'   into a list yourself; that would cost more memory, not less, since it
#'   forces eager materialization of many separate objects in the parent
#'   before any forking happens. Default `1` (serial, with a progress bar;
#'   parallel execution drops the live bar).
#'
#' @return A list: `eMLG` (matrix, individuals x groups that have one --
#'   see `compute_unflagged_eMLG`/`min_n_loci_eMLG`), `groups` (data.table:
#'   group_id, Chr, representative, n_loci, score, has_eMLG, members --
#'   always lists EVERY group regardless of eMLG filtering), `pruned`
#'   (character vector of representative markers, one per group -- for
#'   LD-pruning use, unaffected by eMLG filtering), `params` (the
#'   `ld_w_col`/`ld_w_threshold`/`min_n_loci_flag`/`rho`/`distance_threshold`/
#'   `use_cM`/`cM_threshold` this call actually used -- [plot_pruning_comparison()] defaults to
#'   these so a Stage 1 vs Combined comparison can't silently use a
#'   different threshold on each side).
#'
#' @examples
#' \dontrun{
#' stage1 <- ld_complexity_reduction(map = map, LD_decay = ld_decay, rho = 0.5)
#' result <- ld_prune_and_eMLG(
#'   GTs = GTs, stage1 = stage1, ld_w_col = "ld_w_095",
#'   ld_w_threshold = 0.2, LD_decay = ld_decay, rho = 0.95,
#'   score_threshold = 0.80, min_r2 = 0.2
#' )
#' pruned_markers <- result$pruned
#' eMLG <- result$eMLG
#' }
#'
#' @export
ld_prune_and_eMLG <- function(GTs, stage1, ld_w_col, ld_w_threshold,
                               LD_decay = NULL, rho = 0.95,
                               score_threshold = 0.80, min_r2 = 0.2,
                               distance_threshold = NULL,
                               genetic_map = NULL, cM_threshold = NULL,
                               compute_unflagged_eMLG = TRUE,
                               min_n_loci_eMLG = 1,
                               min_n_loci_flag = Inf,
                               cores = 1) {

  map_snp  <- stage1$map_snp
  clusters <- stage1$clusters

  ## genetic_map/cM_threshold, when supplied, take over run-splitting
  ## entirely (bp-based distance_threshold/rho/LD_decay is not consulted at
  ## all in that case) -- physical distance is a poor, sometimes wildly
  ## misleading proxy for recombination distance: checked directly on real
  ## data, some marker-cluster pairs over 1 Mb apart physically sat at
  ## cM distance 0 (fully linked -- no measurable recombination between
  ## them at all), while others only a few hundred kb apart already spanned
  ## several cM. A single cM_threshold (not per-chromosome, unlike the
  ## bp-based default) is appropriate here precisely because recombination
  ## distance is already the thing that varies by chromosome/region --
  ## physical units have nothing further to correct for once you're
  ## measuring in cM.
  use_cM <- !is.null(genetic_map)
  if (use_cM && is.null(cM_threshold)) {
    stop("`cM_threshold` must be supplied when `genetic_map` is supplied.")
  }

  ## distance_threshold defaults to a per-chromosome, rho-derived value
  ## (d_from_rho(a_pred, rho)) rather than one fixed bp constant, reusing
  ## the same rho already implied by ld_w_col's naming convention (e.g.
  ## "ld_w_095") -- one rho value then governs both what counts as
  ## "locally elevated" and what counts as "still physically contiguous".
  ## Deliberately uses a_pred (chromosome-size-corrected decay rate), not
  ## each chromosome's own observed a: checked directly on real data (an
  ## anomalous ~100x-slower-decay chromosome vs. a typical one, similar
  ## size), using observed a would blow the anomalous chromosome's window
  ## up to ~1.5Mb -- more than double the largest gap actually observed on
  ## it, i.e. license to bridge essentially the whole chromosome into one
  ## run regardless of its real local structure. a_pred instead gives both
  ## chromosomes a near-identical, well-calibrated ~11kb window at
  ## rho=0.95, capturing ~98-99% of real local gaps while still treating
  ## the rare large gaps as presumptive breaks -- exactly the chromosome
  ## whose internal structure is most worth not smoothing over by
  ## construction doesn't get a free pass just because its own decay is
  ## anomalous.
  dist_threshold_by_chr <- NULL
  if (!use_cM && is.null(distance_threshold)) {
    if (is.null(LD_decay)) {
      stop(
        "Either `distance_threshold`, `LD_decay` (to derive a per-chromosome ",
        "threshold from `rho` via d_from_rho(a_pred, rho)), or `genetic_map` + ",
        "`cM_threshold` must be supplied."
      )
    }
    dsum <- LD_decay$decay_sum
    dist_threshold_by_chr <- stats::setNames(d_from_rho(dsum$a_pred, rho), dsum$Chr)
  }

  needs_merge_ids <- map_snp[get(ld_w_col) > ld_w_threshold, unique(CL_id)]
  if (is.finite(min_n_loci_flag)) {
    needs_merge_ids <- union(needs_merge_ids, clusters[n_snps >= min_n_loci_flag, CL_id])
  }
  flagged   <- clusters[CL_id %in% needs_merge_ids]
  unflagged <- clusters[!CL_id %in% needs_merge_ids]

  ## unflagged clusters pass straight through unchanged -- their own
  ## core_snp is already the representative regardless of whether an eMLG
  ## is computed for them at all.
  if (nrow(unflagged)) {
    eligible <- compute_unflagged_eMLG & unflagged$n_snps >= min_n_loci_eMLG
    unflagged_for_eMLG <- unflagged[eligible]

    unflagged_eMLG <- if (nrow(unflagged_for_eMLG)) {
      n_u <- nrow(unflagged_for_eMLG)
      message("Computing eMLGs for ", n_u, " unflagged clusters...")
      m_list <- if (cores > 1 && .Platform$OS.type != "windows") {
        ## fork-based: GTs is shared via copy-on-write, not duplicated per
        ## worker, as long as each worker only reads it (consensus_dosage
        ## never modifies GTs) -- no progress bar under mclapply.
        parallel::mclapply(
          seq_len(n_u),
          function(i) consensus_dosage(GTs, unflagged_for_eMLG$members[[i]]),
          mc.cores = cores
        )
      } else {
        pb <- utils::txtProgressBar(min = 0, max = n_u, style = 3)
        out <- vector("list", n_u)
        for (i in seq_len(n_u)) {
          out[[i]] <- consensus_dosage(GTs, unflagged_for_eMLG$members[[i]])
          utils::setTxtProgressBar(pb, i)
        }
        close(pb)
        out
      }
      m <- do.call(cbind, m_list)
      colnames(m) <- paste0("U", unflagged_for_eMLG$CL_id)
      m
    } else {
      matrix(numeric(0), nrow = nrow(GTs), ncol = 0)
    }

    scores <- setNames(rep(NA_real_, nrow(unflagged)), unflagged$CL_id)
    if (nrow(unflagged_for_eMLG)) {
      scores[as.character(unflagged_for_eMLG$CL_id)] <- apply(unflagged_eMLG, 2, score_eMLG)
    }

    unflagged_groups <- data.table::data.table(
      group_id = paste0("U", unflagged$CL_id),
      Chr = unflagged$Chr,
      representative = unflagged$core_snp,
      n_loci = unflagged$n_snps,
      score = scores[as.character(unflagged$CL_id)],
      has_eMLG = eligible,
      members = unflagged$members
    )
  } else {
    unflagged_eMLG <- matrix(numeric(0), nrow = nrow(GTs), ncol = 0)
    unflagged_groups <- data.table::data.table(
      group_id = character(0), Chr = character(0), representative = character(0),
      n_loci = integer(0), score = numeric(0), has_eMLG = logical(0), members = list()
    )
  }

  params <- list(
    ld_w_col = ld_w_col, ld_w_threshold = ld_w_threshold, min_n_loci_flag = min_n_loci_flag,
    rho = rho, distance_threshold = distance_threshold,
    use_cM = use_cM, cM_threshold = cM_threshold
  )

  if (nrow(flagged) == 0) {
    return(list(
      eMLG = unflagged_eMLG,
      groups = unflagged_groups,
      pruned = unflagged_groups$representative,
      params = params
    ))
  }

  ## flagged clusters' own eMLGs are needed as input to the dynamic cut
  ## regardless of size or min_n_loci_eMLG -- even a single-marker flagged
  ## cluster might merge into something larger, so this can't be skipped
  ## the way unflagged eMLGs can be. consensus_dosage() still applies the
  ## n_loci==1 fast path for the ones that stay singletons.
  n_f <- nrow(flagged)
  message("Computing eMLGs for ", n_f, " flagged clusters...")
  pb <- utils::txtProgressBar(min = 0, max = n_f, style = 3)
  flagged_eMLG_list <- vector("list", n_f)
  for (i in seq_len(n_f)) {
    flagged_eMLG_list[[i]] <- consensus_dosage(GTs, flagged$members[[i]])
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  flagged_eMLG <- do.call(cbind, flagged_eMLG_list)
  colnames(flagged_eMLG) <- flagged$CL_id

  n_loci_flagged <- setNames(flagged$n_snps, flagged$CL_id)

  pos_min <- setNames(
    vapply(flagged$members, function(mk) min(map_snp[marker %in% mk, Pos]), numeric(1)),
    flagged$CL_id
  )
  pos_max <- setNames(
    vapply(flagged$members, function(mk) max(map_snp[marker %in% mk, Pos]), numeric(1)),
    flagged$CL_id
  )

  ## when using a genetic map, run-splitting works in cM space instead of
  ## bp: convert each flagged cluster's physical span to genetic position
  ## once here, so split_by_distance() below (unchanged, unit-agnostic)
  ## just compares cM gaps against cM_threshold. A cluster on a chromosome
  ## genetic_map doesn't cover (or covers with <2 points) gets NA cM
  ## positions, which split_by_distance() already treats as "start a new
  ## run" -- i.e. it conservatively never merges across an unmapped region
  ## rather than erroring.
  if (use_cM) {
    chr_by_cl <- setNames(flagged$Chr, flagged$CL_id)
    pos_min <- setNames(interpolate_cM(genetic_map, chr_by_cl[names(pos_min)], pos_min), names(pos_min))
    pos_max <- setNames(interpolate_cM(genetic_map, chr_by_cl[names(pos_max)], pos_max), names(pos_max))
  }

  flagged_group_list <- list()

  ## restricted to same chromosome AND same physically-contiguous run --
  ## never merge across chromosomes (no physical basis) or across a real
  ## position gap (possible epistatic-selection signal, not redundancy).
  ## Progress is tracked per chromosome, not per run: run count per
  ## chromosome is a poor proxy for time remaining (cost is dominated by
  ## the largest run's O(n^2) cor(), and a chromosome can have very few but
  ## very large runs), so a per-chromosome message reports run count/sizes
  ## for context instead of trying to bar-track them.
  chr_levels <- unique(flagged$Chr)
  message("Distance-restricted dynamic cut: ", length(chr_levels), " chromosomes")
  pb <- utils::txtProgressBar(min = 0, max = length(chr_levels), style = 3)

  for (chr_i in seq_along(chr_levels)) {
    ch <- chr_levels[chr_i]
    ids_ch <- as.character(flagged$CL_id[flagged$Chr == ch])
    dt_ch <- if (use_cM) {
      cM_threshold
    } else if (is.null(dist_threshold_by_chr)) {
      distance_threshold
    } else {
      dist_threshold_by_chr[[ch]]
    }
    runs <- split_by_distance(pos_min[ids_ch], pos_max[ids_ch], dt_ch)
    run_sizes <- table(runs)
    message(
      ch, ": ", length(ids_ch), " flagged clusters -> ", length(run_sizes), " run(s) ",
      "(largest ", max(run_sizes), ")"
    )

    for (r in unique(runs)) {
      run_ids <- names(runs)[runs == r]

      if (length(run_ids) == 1) {
        cid <- run_ids
        flagged_group_list[[length(flagged_group_list) + 1]] <- list(
          Chr = ch, cl_ids = cid, emlg = flagged_eMLG[, cid],
          n_loci = n_loci_flagged[[cid]]
        )
        next
      }

      sub_groups <- dynamic_cut_eMLG(
        flagged_eMLG[, run_ids, drop = FALSE], n_loci_flagged,
        threshold = score_threshold, min_r2 = min_r2
      )
      for (g in sub_groups) {
        flagged_group_list[[length(flagged_group_list) + 1]] <- list(
          Chr = ch, cl_ids = g$members, emlg = g$emlg,
          n_loci = sum(n_loci_flagged[g$members])
        )
      }
    }
    utils::setTxtProgressBar(pb, chr_i)
  }
  close(pb)

  flagged_groups <- data.table::rbindlist(lapply(seq_along(flagged_group_list), function(i) {
    g <- flagged_group_list[[i]]
    data.table::data.table(
      group_id = paste0("F", i),
      Chr = g$Chr,
      representative = pick_representative(g$cl_ids, flagged_eMLG, flagged),
      n_loci = g$n_loci,
      score = score_eMLG(g$emlg),
      has_eMLG = g$n_loci >= min_n_loci_eMLG,
      members = list(unlist(flagged$members[flagged$CL_id %in% g$cl_ids]))
    )
  }))

  ## the merged eMLG signal is already computed as a byproduct of the
  ## dynamic cut regardless of size, so min_n_loci_eMLG only trims which
  ## columns are kept in the returned matrix here, not what gets computed
  keep <- flagged_groups$has_eMLG
  flagged_eMLG_final <- do.call(cbind, lapply(flagged_group_list[keep], `[[`, "emlg"))
  colnames(flagged_eMLG_final) <- flagged_groups$group_id[keep]

  groups <- data.table::rbindlist(list(unflagged_groups, flagged_groups))

  list(
    eMLG = cbind(unflagged_eMLG, flagged_eMLG_final),
    groups = groups,
    pruned = groups$representative,
    params = params
  )
}
