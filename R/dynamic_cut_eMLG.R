# Dynamically consolidate cluster-level eMLGs into coarser groups.
#
# ld_complexity_reduction() reduces a chromosome's raw markers to a much
# smaller set of within-chromosome LD clusters, using complete linkage
# throughout -- appropriate for LD-pruning, where a strict "every pairwise
# r2 must clear the threshold" guarantee matters. But complete linkage is
# the wrong tool for the different question this file answers: which of
# THOSE clusters' consensus genotypes (eMLGs) still belong
# together as one coherent, hard-callable unit, even if their internal
# correlation is diffuse/moderate rather than uniformly high (the signature
# of a polygenic or epistatically-linked complex, not a single strongly-
# linked block).
#
# Validated on real data (a 292-cluster chromosome subset) before writing
# this:
#   - single linkage: catastrophic chaining (292 clusters -> 1 giant group
#     at r2>0.10; still an 81%-of-everything mega-group at r2>0.35)
#   - complete linkage: fragments genuinely coherent diffuse blocks (one
#     45-cluster block with mean internal r2=0.42 split into 31 separate
#     final clusters at the LD-pruning-calibrated rho=0.5 threshold)
#   - average linkage: recovers that same block as ONE cluster, and does so
#     across a much wider/more robust threshold range than complete linkage
#     (which only works in a narrow ~0.10-0.20 window before fragmenting
#     again by 0.25)
#   - a single global r2 threshold (cutree at one height) isn't the right
#     lever anyway, once the actual quality metric that matters is
#     score_eMLG = cor(round(eMLG), eMLG)^2 -- i.e. does the cluster's
#     consensus genotype survive the hard-calling downstream LD/Ohta
#     statistics require. A single height cut doesn't target that directly.
#   - so: walk the average-linkage dendrogram bottom-up (tips -> root, one
#     hclust merge event at a time) and accept each candidate merge only if
#     BOTH (a) the two sides being merged are themselves correlated
#     (min_r2), and (b) the resulting consensus still scores >= threshold on
#     score_eMLG. This lets the stopping point adapt locally per branch,
#     rather than assuming one global r2 threshold is right everywhere.
#   - single linkage was also tried as the TREE structure for this walk
#     (not as the accept/reject criterion, which is score_eMLG/min_r2
#     either way) on the theory that its finer-grained, one-at-a-time
#     candidate steps might give the quality gate more precise control.
#     Empirically this was worse, not better: single linkage's tree
#     fragmented the same 45-cluster block into 7 pieces (vs average
#     linkage's 1), because single linkage proposes candidates by
#     closest-pair-anywhere, which need not fit the group's overall signal
#     well and gets rejected (freezing prematurely) more often than average
#     linkage's already-typical-fit candidates.
#
# Two known bugs fixed during testing, both worth remembering if this is
# ever reimplemented:
#   1. When one side of a candidate merge is already frozen (from an earlier
#      rejected merge) but the OTHER side is still active (a fresh leaf, or
#      an accepted-merge node not yet finalized), that active side must be
#      finalized immediately -- it has nowhere left to grow once its sibling
#      is closed. Missing this silently dropped ~90% of markers (found via a
#      total-member-count check: only 32 of 292 markers survived without
#      this fix).
#   2. score_eMLG alone does not verify the two sides being merged are
#      actually correlated with each other -- two near-monomorphic,
#      unrelated groups can average to something that still rounds cleanly
#      by coincidence (confirmed directly: a merge with min pairwise r2 as
#      low as 1.8e-05 still passed score_eMLG >= 0.80 without the min_r2
#      gate). Hence the min_r2 check alongside score_eMLG.

#' Dynamically cluster eMLG signals via a quality-gated tree walk
#'
#' Consolidates cluster-level consensus genotypes (eMLGs -- e.g. the
#' `eMLG` matrix from [make_eMLGs()]) into coarser groups, by walking an
#' average-linkage dendrogram bottom-up (tips -> root) and accepting each
#' candidate merge only while the result stays hard-callable. This adapts
#' locally per branch, unlike a single global height cut (`cutree()`), which
#' assumes one r2 threshold is right everywhere -- see the file-level
#' comment for why average linkage (not single or complete) was chosen as
#' the tree, and why the accept/reject decision is driven by
#' `score_eMLG`/`min_r2` directly rather than by tree height.
#'
#' [ld_prune_and_eMLG()] is the usual entry point for applying this to a
#' whole [ld_complexity_reduction()] result; call this directly when you
#' already have an eMLG matrix (from [make_eMLGs()] or elsewhere) and just
#' want to consolidate it.
#'
#' @param eMLG Numeric matrix, individuals x clusters (consensus genotype
#'   per cluster). Column names are used as cluster IDs and must be unique.
#' @param n_loci Named numeric vector, cluster ID -> number of raw markers
#'   in that cluster. Used to weight each side when averaging two groups'
#'   consensus signals together (larger, more evidence-backed clusters
#'   count for more). Must have an entry for every `colnames(eMLG)`.
#' @param threshold Minimum `score_eMLG` (`cor(round(x), x)^2` -- how well
#'   the hard-called genotype round-trips against the continuous consensus)
#'   a candidate merge's resulting signal must clear to be accepted.
#'   Validated around 0.80-0.85 on real data; higher is stricter (less
#'   merging).
#' @param min_r2 Minimum r2 required directly between the two sides being
#'   merged, checked in addition to `threshold`. Default `0.2`.
#'
#' @return A list of final groups, each `list(members = <cluster IDs>,
#'   emlg = <consensus genotype vector>, score = <score_eMLG of that
#'   consensus>)`. Every input cluster ID appears in exactly one group.
#'
#' @examples
#' \dontrun{
#' groups <- dynamic_cut_eMLG(eMLGs$eMLG, n_loci = eMLGs$map_eMLG$n_loci, threshold = 0.80)
#' }
#'
#' @export
dynamic_cut_eMLG <- function(eMLG, n_loci, threshold, min_r2 = 0.2) {
  cluster_ids <- colnames(eMLG)
  if (is.null(cluster_ids) || anyDuplicated(cluster_ids)) {
    stop("eMLG must have unique column names (cluster IDs).")
  }
  if (!all(cluster_ids %in% names(n_loci))) {
    stop("n_loci is missing an entry for some column of eMLG.")
  }

  R2 <- suppressWarnings(stats::cor(eMLG, use = "pairwise.complete.obs")^2)
  R2[!is.finite(R2)] <- 0
  diag(R2) <- 1

  hc <- stats::hclust(stats::as.dist(1 - R2), method = "average")
  n  <- length(hc$labels)

  node_members <- vector("list", n - 1)
  node_emlg    <- vector("list", n - 1)
  node_nloci   <- numeric(n - 1)
  node_active  <- rep(NA, n - 1)   # TRUE = still growable, FALSE = frozen/closed
  finalized    <- list()

  get_state <- function(idx) {
    if (idx < 0) {
      lf <- hc$labels[-idx]
      list(members = lf, emlg = eMLG[, lf], nloci = n_loci[lf], active = TRUE)
    } else {
      list(members = node_members[[idx]], emlg = node_emlg[[idx]],
           nloci = node_nloci[idx], active = isTRUE(node_active[idx]))
    }
  }

  finalize_state <- function(st) {
    e <- if (is.matrix(st$emlg)) st$emlg[, 1] else st$emlg
    finalized[[length(finalized) + 1]] <<- list(
      members = st$members, emlg = e, score = score_eMLG(e)
    )
  }

  for (i in seq_len(n - 1)) {
    left  <- get_state(hc$merge[i, 1])
    right <- get_state(hc$merge[i, 2])

    if (!left$active || !right$active) {
      ## A frozen side's members were already finalized when it froze; an
      ## active side has nowhere left to go once its sibling is closed, so
      ## it must be finalized now (see bug #1 in the file-level comment).
      if (left$active)  finalize_state(left)
      if (right$active) finalize_state(right)
      node_active[i] <- FALSE
      next
    }

    ## bug #2: gate on the pair's own r2, not just the resulting average's
    ## score_eMLG, so two uncorrelated groups can't merge just because
    ## their average happens to round cleanly by coincidence.
    pair_r2 <- suppressWarnings(stats::cor(left$emlg, right$emlg, use = "pairwise.complete.obs")^2)

    w <- c(left$nloci, right$nloci)
    cand_emlg <- weighted_row_mean(polarize_genotypes(cbind(left$emlg, right$emlg)), w)
    s <- score_eMLG(cand_emlg)

    if (!is.na(pair_r2) && pair_r2 >= min_r2 && !is.na(s) && s >= threshold) {
      node_members[[i]] <- c(left$members, right$members)
      node_emlg[[i]]    <- cand_emlg
      node_nloci[i]     <- sum(w)
      node_active[i]    <- TRUE
    } else {
      finalize_state(left)
      finalize_state(right)
      node_active[i] <- FALSE
    }
  }

  ## the root: if the very last merge was accepted, that's one final group
  if (isTRUE(node_active[n - 1])) {
    finalize_state(list(members = node_members[[n - 1]], emlg = node_emlg[[n - 1]]))
  }

  finalized
}
