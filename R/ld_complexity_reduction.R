#' Reduce a Marker Set to LD-Independent Representatives
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
#' Clustering is two-stage, using `LD_decay$by_chr[[ch]]$el` (the same
#' precomputed edge list [compute_ld_w()] uses -- built with
#' `keep_el = TRUE`, so no pairwise LD is recomputed here):
#' \enumerate{
#'   \item Connected components on the sparse thresholded graph (single
#'     linkage; cheap).
#'   \item WITHIN each component, refined to true complete-linkage
#'     sub-clusters via `hclust()` on a small dense sub-matrix built from
#'     `el`'s already-available r2 values. This second stage matters:
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
#' problem, and -- since Stage 2 already builds the pairwise r2 sub-matrix
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
#' overlapping-window edges to link distant markers together. Stage 2 does
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
#'
#' @param map A `data.table` with at least `Chr` and `marker` columns.
#' @param LD_decay An `ld_decay` object from [compute_LD_decay()], built with
#'   `keep_el = TRUE` (or a valid `el_data_folder`).
#' @param rho Numeric in (0, 1). Decay level used to derive each
#'   chromosome's clustering threshold via [ld_from_rho()]. Default `0.5`.
#' @param cores Number of chromosomes to process in parallel.
#' @param idx Optional integer vector of row indices into `map` (1-based).
#'   If supplied, only these markers are considered at all -- clustered
#'   against each other and returned -- letting the caller decide the marker
#'   set up front (by `ld_w`, MAF, a region of interest, or anything else)
#'   instead of this function taking its own filtering criterion. If `NULL`
#'   (the default), all markers in `map` are used.
#'
#' @return A list with:
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
#' }
#'
#' @export
ld_complexity_reduction <- function(map, LD_decay, rho = 0.5, cores = 1, idx = NULL) {

  if (!is.null(idx)) map <- map[idx]

  chr_levels <- unique(map$Chr)

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

    ds <- LD_decay$decay_sum[Chr == ch]
    if (nrow(ds) != 1) {
      stop("Expected exactly one decay_sum row for chr '", ch, "', found ", nrow(ds),
           ". Check LD_decay$decay_sum$Chr labels.")
    }
    r2_th <- ld_from_rho(b = ds$b, c = ds$c, rho = rho)

    el <- LD_decay$by_chr[[ch]]$el
    if (is.null(el)) {
      stop("LD_decay$by_chr[['", ch, "']]$el is NULL -- compute_LD_decay() must be run with ",
           "keep_el = TRUE (or a valid el_data_folder) for ld_complexity_reduction() to reuse its edges.")
    }
    if (is.character(el)) el <- data.table::fread(el, showProgress = FALSE)  ## el_data_folder mode stores a path
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

    message(ch, ": ", n_comp, " components")
    pb <- utils::txtProgressBar(min = 0, max = n_comp, style = 3)
    on.exit(close(pb))

    per_component <- vector("list", n_comp)

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
        sub_edges <- edges_th[SNP1 %in% cluster_markers & SNP2 %in% cluster_markers]
        m <- length(cluster_markers)
        R2_sub <- matrix(0, m, m, dimnames = list(cluster_markers, cluster_markers))
        diag(R2_sub) <- 1
        i1 <- match(sub_edges$SNP1, cluster_markers)
        i2 <- match(sub_edges$SNP2, cluster_markers)
        R2_sub[cbind(i1, i2)] <- sub_edges$r2
        R2_sub[cbind(i2, i1)] <- sub_edges$r2

        d      <- stats::as.dist(1 - R2_sub)
        hc     <- stats::hclust(d, method = "complete")
        sub_cl <- stats::cutree(hc, h = 1 - r2_th)

        ## same split()-not-rescan reasoning as above, one level down
        data.table::rbindlist(lapply(split(cluster_markers, sub_cl), function(sub_markers) {

          if (length(sub_markers) == 1) {
            return(data.table::data.table(
              marker = sub_markers, is_core = TRUE,
              median_ld = NA_real_, grp = sub_markers
            ))
          }

          if (length(sub_markers) == 2) {
            ld <- R2_sub[sub_markers[1], sub_markers[2]]
            return(data.table::data.table(
              marker = sub_markers, is_core = c(TRUE, FALSE),
              median_ld = ld, grp = sub_markers[1]
            ))
          }

          ## representative: highest median r2 to the rest of the sub-cluster
          ## (see roxygen note above for why, and why this reuses R2_sub
          ## instead of recomputing anything from GTs)
          R2_sc <- R2_sub[sub_markers, sub_markers, drop = FALSE]
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

  list(
    map_snp  = map_snp,
    clusters = clusters,
    pruned   = clusters$core_snp
  )
}
