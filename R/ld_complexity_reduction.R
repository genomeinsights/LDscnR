## Edge list for one chromosome: reuse LD_decay's if it has one (in memory, or a
## path in el_data_folder mode), else rebuild it from genotypes on the fly. The
## rebuild is per chromosome and the result is discarded by the caller as soon as
## its clusters are formed, so nothing is written or accumulated.
.chr_edge_list <- function(LD_decay, ch, chr_markers, gds = NULL, reopen_path = NULL,
                           el_floor = 0) {

  el <- LD_decay$by_chr[[ch]]$el
  if (is.character(el)) return(data.table::fread(el, showProgress = FALSE))  ## el_data_folder mode stores a path
  if (!is.null(el)) return(el)

  if (is.null(gds) && is.null(reopen_path)) {
    stop("LD_decay$by_chr[['", ch, "']]$el is NULL -- either run compute_LD_decay() with ",
         "keep_el = TRUE (or a valid el_data_folder), or pass gds = to ",
         "ld_complexity_reduction() so this chromosome's edges can be rebuilt on the fly.")
  }

  g <- gds
  if (!is.null(reopen_path)) {
    g <- SNPRelate::snpgdsOpen(reopen_path, allow.duplicate = TRUE)
    on.exit(SNPRelate::snpgdsClose(g), add = TRUE)
  } else if (is.character(g)) {
    g <- SNPRelate::snpgdsOpen(g, allow.duplicate = TRUE)
    on.exit(SNPRelate::snpgdsClose(g), add = TRUE)
  }

  p <- LD_decay$params
  ## el_floor applies only on this on-the-fly path. A stored edge list (in
  ## memory or on disk) is already materialised, so there is nothing to save.
  get_el(gds          = g,
         SNP_id       = chr_markers,
         slide        = if (!is.null(p$slide))     p$slide     else 1000,
         ld_method    = if (!is.null(p$ld_method)) p$ld_method else "corr",
         cores        = 1,
         by_chr       = TRUE,
         el_floor     = el_floor)
}

#' Reduce a Marker Set to Less-Redundant LD Representatives
#'
#' Clusters markers into local LD blocks and reduces each block to a single
#' representative ("core") marker. Intended as an LD-pruning step ahead of
#' kinship/relatedness matrix estimation (a GRM, EMMAX's `K`, BayPass's
#' `OMEGA`), or as the clustering step feeding [make_eMLGs()] -- it returns
#' `CL_id`/`n_loci` columns on the same per-marker basis `make_eMLGs()`
#' expects from its `map_cl` argument.
#'
#' The clustering threshold is derived per chromosome from that chromosome's
#' own fitted LD-decay curve (`rho`, via [ld_from_rho()]) rather than a fixed
#' r2 cutoff applied genome-wide, so chromosomes with different background LD
#' levels get an appropriately different absolute threshold. This threshold
#' is deliberately NOT allowed to vary further within a chromosome (e.g. by
#' local decay rate): a region with slower local decay -- a centromere or an
#' inversion -- is exactly where pruning needs to be MORE aggressive, since
#' those markers are the most redundant; adapting the threshold to expect
#' high local LD there and tolerate it would defeat the purpose of pruning.
#'
#' Stage-1 clustering has two internal passes, using a per-chromosome edge
#' list of pairwise \eqn{r^2}. These are not the manuscript's Stage 1 and
#' Stage 2: manuscript Stage 2 assembles significant units after testing.
#' The edge list is taken from `LD_decay$by_chr[[ch]]$el` when it
#' is there (the same precomputed edges [compute_ld_w()] uses, built with
#' `keep_el = TRUE`, so no pairwise LD is recomputed); otherwise, if `gds` is
#' supplied, it is rebuilt from the genotypes ONE CHROMOSOME AT A TIME and
#' discarded again, so a caller can run `compute_LD_decay(keep_el = FALSE)` --
#' saving nothing to disk and carrying no edges in the returned object -- and
#' still use this function. The two internal passes are:
#' \enumerate{
#'   \item Connected components on the sparse thresholded graph (single
#'     linkage; cheap).
#'   \item WITHIN each component, refined to true complete-linkage
#'     sub-clusters via `hclust()` on a small dense sub-matrix built from
#'     `el`'s already-available r2 values. This refinement matters:
#'     single linkage alone is vulnerable to "chaining" through a mosaic LD
#'     pattern -- marker A uncorrelated with C, but both bridged into one
#'     component via an intermediate B that's correlated with each -- which
#'     would wrongly treat A and C as redundant. Complete linkage requires
#'     every pairwise r2 within a final cluster to clear the threshold, not
#'     just adjacent ones.
#' }
#'
#' The representative within each final sub-cluster is the marker with the
#' highest median r2 to the rest of its cluster -- not the marker most
#' correlated with the cluster's mean genotype. The mean-genotype approach
#' is fragile to inconsistent allele-coding direction within a cluster (a
#' genuinely central marker coded on the opposite reference allele shows a
#' spuriously negative correlation with the mean and gets passed over for a
#' noisier marker); r2 is sign-agnostic, so median r2 doesn't have that
#' problem, and -- since complete-linkage refinement already builds the pairwise r2 sub-matrix
#' for `hclust()` -- costs nothing extra to compute.
#'
#' Chromosomes are processed with a `txtProgressBar` each, labelled by
#' chromosome, tracking progress through that chromosome's connected
#' components -- some chromosomes can have tens of thousands of them, so a
#' single genome-wide bar would sit at 0% for a long time on the first (or
#' largest) chromosome. With `cores > 1` each chromosome's bar is drawn by
#' its own forked worker and they write to the same stdout, so output can
#' interleave; this is purely cosmetic; use `cores = 1` for clean output.
#'
#' Known limitation: `el` itself is built from a sliding window, so a pair
#' of markers farther apart than that window was never evaluated at all
#' (not "below threshold" -- absent, which behaves as r2 = 0 here). Stage 1
#' tolerates this fine, since single-linkage components only need a chain of
#' overlapping-window edges to link distant markers together. Complete-linkage refinement does
#' not: complete linkage requires every pairwise r2 within a final cluster
#' to individually clear the threshold, so a genuinely one-block region
#' whose extent exceeds the window gets fragmented into more than one
#' representative, purely because the missing pairs default to r2 = 0
#' rather than their real (unknown) value. The window is normally sized off
#' genome/chromosome-average decay to reach a high rho (0.95-0.99), which
#' comfortably covers a lower pruning rho (e.g. 0.5) in typical regions --
#' but centromeres/inversions decay far slower than that average, so this
#' is where fragmentation is most likely. Fixing it would mean recomputing
#' r2 directly from genotypes (unlimited by any window) for markers whose
#' pairwise LD isn't in `el`, which this function does not currently do.
#' [ld_prune_and_eMLG()] fixes it after the fact instead, by comparing
#' clusters (not raw markers) directly from genotypes with no window
#' restriction, restricted to clusters flagged post hoc using this
#' function's own cluster boundaries -- not by pre-splitting markers before
#' clustering, which can sever a real block right at the flagging
#' threshold.
#'
#' @param map A `data.table` with at least `Chr` and `marker` columns.
#' @param LD_decay An `ld_decay` object from [compute_LD_decay()]. Its edge
#'   lists are reused when present (`keep_el = TRUE`, or a valid
#'   `el_data_folder`); otherwise supply `gds` and they are rebuilt on the fly.
#' @param gds Optional. An open GDS object (as returned by
#'   [create_gds_from_geno()]) or a path to a `.gds` file, used to rebuild a
#'   chromosome's edge list when `LD_decay` carries none. `slide` and
#'   `ld_method` are taken from `LD_decay$params`, so the rebuilt edges match
#'   the window the decay curve was fitted on; the one difference is that
#'   rebuilding covers every marker on the chromosome, whereas the original
#'   `el` may have been built on a `max_SNPs_decay` subsample. Ignored when the
#'   edge lists are already present. With `cores > 1` each worker opens its own
#'   handle from the file path; if only an open handle is available and its
#'   path cannot be recovered, `cores` is forced to 1 rather than sharing a GDS
#'   connection across forks.
#' @param r2_rho Numeric in (0, 1). Decay level used to derive each
#'   chromosome's clustering threshold via [ld_from_rho()]. Default `0.5`.
#'   Ignored when `r2_th` is supplied.
#' @param r2_th Optional single number in (0, 1]. A fixed r2 clustering
#'   threshold applied to every chromosome instead of the decay-relative one
#'   derived from `r2_rho`. Default `NULL` keeps the decay-relative threshold.
#' @param cores Number of chromosomes to process in parallel.
#' @param idx Optional integer vector of row indices into `map` (1-based).
#'   If supplied, only these markers are considered at all -- clustered
#'   against each other and returned -- letting the caller decide the marker
#'   set up front (by `ld_w`, MAF, a region of interest, or anything else)
#'   instead of this function taking its own filtering criterion. If `NULL`
#'   (the default), all markers in `map` are used.
#'
#' @return An object of class `"ld_complexity_reduction"` (a list with a
#'   [print.ld_complexity_reduction()] method -- printed automatically once
#'   at the end of the call, and again any time the returned object is
#'   printed later):
#' \describe{
#'   \item{map_snp}{`map` (subset to `idx` if supplied), with `CL_id`,
#'     `n_loci`, `is_core` (logical), and `median_ld` (this marker's own
#'     median r2 to the rest of its cluster, `NA` for singletons) added.
#'     `CL_id`/`n_loci` are named to match what [make_eMLGs()] expects from
#'     `map_cl`.}
#'   \item{clusters}{One row per cluster: `Chr`, `CL_id`, `core_snp`,
#'     `median_ld` (the core marker's), `n_snps`, and `members` (a
#'     list-column of every marker in the cluster, including the core).}
#'   \item{pruned}{Character vector of representative marker names --
#'     `clusters$core_snp`, for direct use as a pruned marker set.}
#'   \item{params}{The resolved `r2_rho`, `r2_th` (`NULL` unless a fixed
#'     threshold was supplied) and `cores` this call used, so a
#'     downstream consumer can read the threshold back rather than have it
#'     passed a second time as a separate argument.}
#' }
#'
#' @param max_dense Largest single-linkage component (in markers) refined with a
#'   dense r2 matrix and [stats::hclust()]. Larger components use an equivalent
#'   sparse complete-linkage pass on the thresholded edges, whose memory scales
#'   with the number of edges rather than the square of the component size (a
#'   dense matrix costs 8 bytes per marker pair: 0.2 GB at 5,000 markers, 29 GB
#'   at 60,000). The sparse pass replays [stats::hclust()]'s own merge order,
#'   so the two give identical clusters, including when r2 values are tied.
#'   Component size grows with `slide` in `LD_decay$params`: raising the slide
#'   (e.g. after an LD-decay check) lets single linkage chain more markers into
#'   one component, so the largest component per chromosome is reported as
#'   each chromosome is clustered. Default `5000`.
#' @param rho Deprecated. Use `r2_rho`.
#' @export
ld_complexity_reduction <- function(map, LD_decay, r2_rho = 0.5, cores = 1, idx = NULL,
                                    gds = NULL, r2_th = NULL, max_dense = 5000L,
                                    rho = NULL) {
  if (!missing(rho)) r2_rho <- .renamed_arg(rho, "rho", "r2_rho")
  if (!(is.numeric(max_dense) && length(max_dense) == 1L && !is.na(max_dense) && max_dense >= 0))
    stop("`max_dense` must be a single non-negative number.")

  if (!is.null(r2_th) &&
      !(is.numeric(r2_th) && length(r2_th) == 1L && is.finite(r2_th) &&
        r2_th > 0 && r2_th <= 1))
    stop("`r2_th` must be NULL or a single number in (0, 1].")
  r2_th_fixed <- r2_th

  if (!is.null(idx)) map <- map[idx]

  ## as.character(): a factor Chr would otherwise be coerced to its integer
  ## level code when used as a list index below (LD_decay$by_chr[[ch]]),
  ## producing a "subscript out of bounds" error on real (factor-Chr) maps.
  chr_levels <- as.character(unique(map$Chr))

  ## Any chromosome without edges has to be rebuilt from `gds`. A forked worker
  ## cannot safely share an open GDS connection, so with cores > 1 each opens its
  ## own from the file path; without a usable path, drop to sequential instead.
  need_el     <- any(vapply(chr_levels, function(ch) is.null(LD_decay$by_chr[[ch]]$el), logical(1)))
  reopen_path <- NULL
  if (need_el && !is.null(gds) && cores > 1) {
    path <- if (is.character(gds)) gds else gds$filename
    if (!is.null(path) && file.exists(path)) {
      reopen_path <- path
    } else {
      message("Rebuilding edge lists from an open GDS handle -- forcing cores = 1 ",
              "(a forked worker cannot share it safely).")
      cores <- 1
    }
  }

  by_chr <- parallel_apply(chr_levels, function(ch) {

    chr_map     <- map[Chr == ch]
    chr_markers <- chr_map$marker

    singletons <- function(markers) {
      if (!length(markers)) {
        return(data.table::data.table(
          marker = character(0), is_core = logical(0),
          median_ld = numeric(0), grp = character(0)
        ))
      }
      data.table::data.table(
        marker = markers, is_core = TRUE,
        median_ld = NA_real_, grp = markers
      )
    }

    if (length(chr_markers) < 2) {
      out <- singletons(chr_markers)
      out[, Chr := ch]
      return(out)
    }

    ## A fixed r2_th, when supplied, replaces the decay-relative threshold on
    ## every chromosome.
    r2_th <- if (!is.null(r2_th_fixed)) r2_th_fixed else {
      ds <- LD_decay$decay_sum[Chr == ch]
      if (nrow(ds) != 1) {
        stop("Expected exactly one decay_sum row for chr '", ch, "', found ", nrow(ds),
             ". Check LD_decay$decay_sum$Chr labels.")
      }
      ld_from_rho(b = ds$b, c = ds$c, rho = r2_rho)
    }

    ## The edge list is consumed on the very next line and nowhere else, and only
    ## to select r2 >= r2_th. Building it with that floor is therefore
    ## OUTPUT-IDENTICAL -- every row it drops is a row the filter below would drop
    ## -- while avoiding the materialisation of pairs that are about to be thrown
    ## away. On a 15,000-marker chromosome at the default sliding window that is
    ## the difference between ~15M rows and the few that clear the threshold.
    ## Guarded: a chromosome whose decay fit failed gives a non-finite threshold,
    ## and then no floor is applied.
    floor_ok <- length(r2_th) == 1L && is.finite(r2_th) && r2_th >= 0 && r2_th <= 1
    el <- .chr_edge_list(LD_decay, ch, chr_markers, gds = gds, reopen_path = reopen_path,
                         el_floor = if (floor_ok) r2_th else 0)
    edges_th <- el[r2 >= r2_th & SNP1 %in% chr_markers & SNP2 %in% chr_markers]

    if (nrow(edges_th) == 0) {
      out <- singletons(chr_markers)
      out[, Chr := ch]
      return(out)
    }

    ## Stage 1 (cheap): connected components from the sparse thresholded
    ## graph -- single-linkage groups
    g <- igraph::graph_from_data_frame(
      edges_th[, .(SNP1, SNP2)], directed = FALSE,
      vertices = data.table::data.table(name = chr_markers)
    )
    comp <- igraph::components(g)

    ## split() groups markers by component in one pass; rescanning
    ## comp$membership once per component (the more obvious
    ## unique(membership) + membership==cl_id pattern) is O(n_markers *
    ## n_components), which dominates runtime once a chromosome has tens of
    ## thousands of components.
    components_list <- split(names(comp$membership), comp$membership)
    n_comp <- length(components_list)

    message(ch, ": ", n_comp, " components (largest ", max(lengths(components_list)),
            " markers)")
    pb <- utils::txtProgressBar(min = 0, max = n_comp, style = 3)
    on.exit(close(pb))

    per_component <- vector("list", n_comp)

    ## Every threshold edge lies inside exactly one component (components are
    ## built from these edges), so split the edges by component ONCE instead of
    ## rescanning all of edges_th for every component (O(n_edges * n_components)).
    ## Output-identical: each component sees exactly the rows the old
    ## `edges_th[SNP1 %in% cm & SNP2 %in% cm]` filter returned, in the same order.
    edge_comp <- comp$membership[edges_th$SNP1]
    edge_rows <- split(seq_len(nrow(edges_th)), factor(edge_comp, levels = names(components_list)))

    for (i in seq_len(n_comp)) {
      cluster_markers <- components_list[[i]]

      per_component[[i]] <- if (length(cluster_markers) == 1) {
        data.table::data.table(
          marker = cluster_markers, is_core = TRUE,
          median_ld = NA_real_, grp = cluster_markers
        )
      } else {
        ## Stage 2: within this (single-linkage) component, refine to TRUE
        ## complete-linkage sub-clusters -- a component can still be a
        ## "chain" where not every pair inside it clears r2_th
        sub_edges <- edges_th[edge_rows[[i]]]
        m <- length(cluster_markers)
        i1 <- match(sub_edges$SNP1, cluster_markers)
        i2 <- match(sub_edges$SNP2, cluster_markers)

        if (m <= max_dense) {
          R2_sub <- matrix(0, m, m, dimnames = list(cluster_markers, cluster_markers))
          diag(R2_sub) <- 1
          R2_sub[cbind(i1, i2)] <- sub_edges$r2
          R2_sub[cbind(i2, i1)] <- sub_edges$r2

          d      <- stats::as.dist(1 - R2_sub)
          hc     <- stats::hclust(d, method = "complete")
          sub_cl <- stats::cutree(hc, h = 1 - r2_th)
          r2_of  <- function(sm) R2_sub[sm, sm, drop = FALSE]
        } else {
          ## Large component: the dense m x m matrix (8 m^2 bytes, plus dist())
          ## is what exhausts memory at full marker density. Same complete-linkage
          ## cut, computed on the sparse threshold graph (see .sparse_complete_linkage).
          sub_cl <- .sparse_complete_linkage(m, i1, i2, 1 - sub_edges$r2, 1 - r2_th)
          ## edges inside each final sub-cluster, split once (same reasoning as edge_rows)
          within <- which(sub_cl[i1] == sub_cl[i2])
          sub_rows <- split(within, sub_cl[i1][within])
          r2_of  <- function(sm) {
            k <- match(sm, cluster_markers)
            R <- matrix(0, length(sm), length(sm), dimnames = list(sm, sm)); diag(R) <- 1
            w <- sub_rows[[as.character(sub_cl[k[1]])]]
            if (length(w)) {
              a1 <- match(i1[w], k); a2 <- match(i2[w], k)
              R[cbind(a1, a2)] <- sub_edges$r2[w]; R[cbind(a2, a1)] <- sub_edges$r2[w]
            }
            R
          }
        }

        ## same split()-not-rescan reasoning as above, one level down
        data.table::rbindlist(lapply(split(cluster_markers, sub_cl), function(sub_markers) {

          if (length(sub_markers) == 1) {
            return(data.table::data.table(
              marker = sub_markers, is_core = TRUE,
              median_ld = NA_real_, grp = sub_markers
            ))
          }

          if (length(sub_markers) == 2) {
            ld <- r2_of(sub_markers)[1, 2]
            return(data.table::data.table(
              marker = sub_markers, is_core = c(TRUE, FALSE),
              median_ld = ld, grp = sub_markers[1]
            ))
          }

          ## representative: highest median r2 to the rest of the sub-cluster
          ## (see roxygen note above for why, and why this reuses R2_sub
          ## instead of recomputing anything from GTs)
          R2_sc <- r2_of(sub_markers)
          diag(R2_sc) <- NA  ## exclude self-similarity from each marker's median
          med_r2 <- apply(R2_sc, 1, stats::median, na.rm = TRUE)
          core <- sub_markers[which.max(med_r2)]

          data.table::data.table(
            marker = sub_markers,
            is_core = sub_markers == core,
            median_ld = as.numeric(med_r2),
            grp = core
          )
        }))
      }

      utils::setTxtProgressBar(pb, i)
    }

    out <- data.table::rbindlist(per_component)
    out[, Chr := ch]
    out
  }, cores = cores)

  long <- data.table::rbindlist(by_chr)
  long[, CL_id := .GRP, by = .(Chr, grp)]
  long[, n_loci := .N, by = CL_id]
  long[, grp := NULL]

  ## X[i, on=] (not merge(), which sorts its output by the join key) so
  ## map_snp keeps the same row order as the input map -- callers relying on
  ## positional alignment with map (e.g. against a GTs matrix of their own)
  ## shouldn't have that silently broken.
  map_snp <- long[, .(marker, CL_id, n_loci, is_core, median_ld)][
    data.table::copy(map), on = "marker"
  ]

  clusters <- long[, .(
    core_snp  = marker[is_core],
    median_ld = median_ld[is_core],
    n_snps    = .N,
    members   = list(marker)
  ), by = .(Chr, CL_id)]

  ## `params` records the resolved call so a downstream consumer -- ld_outlier_test()
  ## in particular -- can read `r2_rho` back rather than have it passed a second time as a
  ## separate argument that could silently disagree with what actually built this object.
  ## Same pattern ld_prune_and_eMLG() already uses for its own resolved thresholds.
  out <- list(
    map_snp  = map_snp,
    clusters = clusters,
    pruned   = clusters$core_snp,
    params   = list(r2_rho = r2_rho, r2_th = r2_th_fixed, cores = cores,
                    ## `rho` is the pre-0.9.2 name for `r2_rho`, kept so callers
                    ## that read `params$rho` (LDscnR-multi up to 0.9.1) still
                    ## get the value instead of NULL.
                    rho = r2_rho)
  )
  class(out) <- "ld_complexity_reduction"
  print(out)
  out
}

## Complete-linkage clustering cut at height h, computed on the SPARSE graph of
## pairs with distance <= h (d = 1 - r2), without a dense m x m matrix, giving the
## same partition as cutree(hclust(as.dist(D), "complete"), h = h) on the dense
## matrix D in which every absent pair is 1 (r2 = 0), INCLUDING under ties.
##
## Why the sparse graph is enough: complete linkage is monotone, so every merge at
## height <= h happens before any merge above h, and cutree(h) keeps exactly the
## merges at height <= h. Two clusters merge at height <= h only if every cross
## pair has d <= h; an absent pair (d = 1 > h, for any r2 threshold > 0) blocks the
## merge for good. The complete-linkage update max(D(I2,K), D(J2,K)) therefore
## leaves only the COMMON neighbours of the two merged clusters.
##
## Why the partitions are identical even with ties (30 diploids give many equal
## r2 values): this replays the merge order of R's own hclust() Fortran routine
## (stats/src/hclust.f, Murtagh's nearest-neighbour-list algorithm), not just any
## valid complete-linkage order: NN(i) is the first j > i with the smallest D(i, j);
## each step merges the pair from the first i with the smallest DISNN(i); the
## merged cluster keeps the lower index I2; NN(I2) is recomputed over K > I2; and
## every i whose NN was I2 or J2 is redetermined. Values above h (absent pairs) only
## break ties among merges above h, which cutree(h) never uses, so they can be
## treated as +Inf. Verified against hclust() on random banded graphs with heavy
## ties and on every component of real full-density chromosomes.
##
## Memory is O(number of edges), not O(m^2): a 60,547-marker component that needs
## 29 GB as a dense matrix (plus 15 GB for dist()) fits in tens of MB here.
##
## @param m Number of markers (vertices are 1..m, in the dense matrix's order).
## @param i1,i2 Integer endpoints of the edges (1-based; self-pairs are ignored and a
##   repeated pair keeps its last distance).
## @param d Distances 1 - r2 of those edges.
## @param h Cut height.
## @return Integer vector of length m: cluster label per marker.
.sparse_complete_linkage <- function(m, i1, i2, d, h) {
  keep <- d <= h & i1 != i2
  i1 <- i1[keep]; i2 <- i2[keep]; d <- d[keep]
  ## A pair listed twice (in either orientation) takes its last value, as the
  ## dense matrix's assignment does; the neighbour lists need each pair once.
  dup <- duplicated(cbind(pmin(i1, i2), pmax(i1, i2)), fromLast = TRUE)
  i1 <- i1[!dup]; i2 <- i2[!dup]; d <- d[!dup]
  a <- c(i1, i2); b <- c(i2, i1); dd <- c(d, d)
  nb <- rep(list(integer(0)), m); nd <- rep(list(numeric(0)), m)
  if (length(a)) {
    o <- order(a, b)
    a <- a[o]; b <- b[o]; dd <- dd[o]
    sp <- split(seq_along(a), factor(a, levels = seq_len(m)))
    nb <- lapply(sp, function(k) b[k]); nd <- lapply(sp, function(k) dd[k])
    names(nb) <- names(nd) <- NULL
  }
  ## neighbour lists are kept sorted by neighbour index, so "first j with the
  ## minimum" is which.min over the j > i part
  nn_after <- function(i) {
    k <- nb[[i]]; v <- nd[[i]]
    s <- k > i
    if (!any(s)) return(c(NA_real_, Inf))
    k <- k[s]; v <- v[s]
    w <- which.min(v)
    c(k[w], v[w])
  }
  nn <- rep(NA_integer_, m); disnn <- rep(Inf, m)
  for (i in seq_len(m)) { r <- nn_after(i); nn[i] <- r[1]; disnn[i] <- r[2] }
  active <- rep(TRUE, m)
  members <- as.list(seq_len(m))
  repeat {
    im <- which.min(disnn)                 # first i with the smallest DISNN (inactive = Inf)
    if (!length(im) || disnn[im] > h) break
    jm <- nn[im]
    i2 <- min(im, jm); j2 <- max(im, jm)
    ## merge J2 into I2 (complete linkage: common neighbours, larger distance)
    p_k <- nb[[i2]]; p_d <- nd[[i2]]; q_k <- nb[[j2]]; q_d <- nd[[j2]]
    sp_ <- p_k != j2; p_k <- p_k[sp_]; p_d <- p_d[sp_]
    sq_ <- q_k != i2; q_k <- q_k[sq_]; q_d <- q_d[sq_]
    common <- p_k[p_k %in% q_k]            # stays sorted
    newd <- pmax(p_d[match(common, p_k)], q_d[match(common, q_k)])
    for (cc in union(p_k, q_k)) {
      kc <- nb[[cc]]; dc <- nd[[cc]]
      rm_ <- kc == i2 | kc == j2
      kc <- kc[!rm_]; dc <- dc[!rm_]
      j <- match(cc, common)
      if (!is.na(j)) {
        pos <- findInterval(i2, kc)
        kc <- append(kc, i2, pos); dc <- append(dc, newd[j], pos)
      }
      nb[[cc]] <- kc; nd[[cc]] <- dc
    }
    nb[[i2]] <- common; nd[[i2]] <- newd
    nb[j2] <- list(integer(0)); nd[j2] <- list(numeric(0))
    active[j2] <- FALSE; disnn[j2] <- Inf; nn[j2] <- NA_integer_
    members[[i2]] <- c(members[[i2]], members[[j2]]); members[j2] <- list(NULL)
    r <- nn_after(i2); nn[i2] <- r[1]; disnn[i2] <- r[2]
    ## redetermine NN of every i whose NN was I2 or J2
    for (i in which(active & !is.na(nn) & (nn == i2 | nn == j2) & seq_len(m) != i2)) {
      r <- nn_after(i); nn[i] <- r[1]; disnn[i] <- r[2]
    }
  }
  label <- integer(m)
  for (k in which(active)) label[members[[k]]] <- k
  label
}

#' Print an ld_complexity_reduction Object
#'
#' Summarizes the result of [ld_complexity_reduction()]: how many markers
#' were kept as cluster representatives, the cluster-size distribution, and
#' the representative markers' median r2 to the rest of their cluster.
#'
#' @param x An `"ld_complexity_reduction"` object.
#' @param digits Number of significant digits for printed statistics.
#' @param ... Unused; present for S3 consistency.
#'
#' @return `x`, invisibly.
#'
#' @export
print.ld_complexity_reduction <- function(x, digits = 3, ...) {

  cat("<ld_complexity_reduction>\n")

  n_markers  <- nrow(x$map_snp)
  n_clusters <- nrow(x$clusters)

  cat("\nMarkers:\n")
  cat("  input markers:", format(n_markers, big.mark = ","), "\n")
  cat(
    "  clusters / pruned markers:", format(n_clusters, big.mark = ","),
    " (", signif(100 * n_clusters / n_markers, digits), "% kept)\n",
    sep = ""
  )
  cat("  chromosomes:", data.table::uniqueN(x$clusters$Chr), "\n")
  if (!is.null(x$params$r2_th)) {
    cat("  threshold: fixed r2 =", signif(x$params$r2_th, digits), "\n")
  } else if (!is.null(.stage1_r2_rho(x))) {
    cat("  threshold: decay-relative, r2_rho =", signif(.stage1_r2_rho(x), digits), "\n")
  }

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

  cat("\nStored components:\n")
  cat("  map_snp, clusters, pruned\n")

  invisible(x)
}
