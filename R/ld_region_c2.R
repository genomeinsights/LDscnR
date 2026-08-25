#' Second-tier C-score -- integrate tau_C and l_min away
#'
#' Once the null only assigns per-region significance, `tau` and `l_min` stop
#' being things to calibrate and become nuisance parameters -- exactly like the
#' `ld_w` window `rho` and the stringency `q*` that the per-marker C-score
#' ([ld_cscore()]) already integrates over. This applies the same remedy one level
#' up: run [ld_region_scan()] across a `(tau, l_min)` grid and score each locus by
#' the fraction of the grid in which it comes out significant.
#'
#' \deqn{U = \{g \in G : g \text{ yields} \ge 1 \text{ region with } q_R < \alpha_R\}}
#' \deqn{C^{(2)}_\ell = \frac{1}{|U|} \sum_{g \in U}
#'       \mathbf{1}[\exists R \ni \ell \text{ with } q_R(g) < \alpha_R]}
#'
#' Normalising by the *usable* grid `U` -- the cells that produce any significant
#' region at all -- is what defeats cherry-picking. A locus significant only in the
#' single cell that happens to maximise the region count scores `1/|U|`, near zero,
#' while one significant across most of the grid scores near 1. No cell is
#' privileged, so no cell can be selected after the fact.
#'
#' `C^(2)` is a **ranking**, not a second significance test: the grid is not
#' independent across cells, and no error rate is controlled by it. Its use is to
#' order loci that `q_R` cannot separate -- which is the normal situation, since
#' with a permutation null every region worth reporting tends to sit at the
#' p-floor `1/(1+B)`.
#'
#' @section Relation to the per-SNP C-score:
#' The consistency principle is Fang et al.'s (2021, \doi{10.1093/molbev/msab144});
#' see [ld_cscore()] for the attribution. Applying it a second time, to the
#' region-calling parameters, is what is new here, and it is only possible
#' because [ld_cscore()] keeps the first-level score per-SNP: regions can then be
#' re-formed at every `(tau, l_min)` rather than fixed once.
#'
#' @param null An `ld_null` object ([ld_null_from_p()] or [structured_null()]).
#' @param edges An [ld_edges()] object built over `null$universe`.
#' @param tau_grid C-score thresholds to sweep (default `seq(0.02, 0.5, by = 0.02)`).
#' @param lmin_grid Region sizes to sweep (default `c(1, 2, 3, 5, 10, 15, 20)`).
#' @param fdr BH level defining significance within a cell (default 0.05).
#' @param anchor_tau,anchor_lmin The operating point whose regions are annotated
#'   with `C^(2)` (defaults 0.05 and 3). The anchor supplies clean, LD-supported
#'   region boundaries; `C^(2)` is insensitive to the exact choice, since a real
#'   locus scores high whether its boundaries are drawn here or in a neighbouring
#'   cell.
#' @param verbose Print per-`tau` progress.
#'
#' @return An `ld_region_c2` object: `regions` (the anchor regions with a `c2`
#'   column, ordered by it), `landscape` (per-cell `tau`, `l_min`, `n_obs`,
#'   `n_sig`), `n_usable` (`|U|`), `n_cells`, and `params`.
#' @seealso [ld_region_scan()], [ld_gate()], [ld_scan()]
#' @export
ld_region_c2 <- function(null, edges, tau_grid = seq(0.02, 0.5, by = 0.02),
                         lmin_grid = c(1L, 2L, 3L, 5L, 10L, 15L, 20L), fdr = 0.05,
                         anchor_tau = 0.05, anchor_lmin = 3L, verbose = TRUE) {
  stopifnot(inherits(null, "ld_null"))
  co <- .edge_coords(edges); B <- null$B
  cells <- vector("list", length(tau_grid) * length(lmin_grid)); k <- 0L
  sig_acc <- list()

  for (tau in tau_grid) {
    ## Clustering does not depend on l_min, so cluster once per tau at l_min = 1
    ## and take the size filter afterwards -- the whole reason a grid this size
    ## is affordable at all.
    O_all <- .region_table(null$C_obs, tau, 1L, edges, co)
    data.table::setnames(O_all, "score", "s_R")
    S_all <- .surr_table(null, tau, 1L, edges, co)
    for (lm in lmin_grid) {
      k <- k + 1L
      O <- O_all[get("size") >= lm]
      S <- if (nrow(S_all)) S_all[get("size") >= lm] else S_all
      res <- .region_pq(O, S, B, fdr)
      cells[[k]] <- data.table::data.table(tau = tau, l_min = lm,
                                           n_obs = nrow(res), n_sig = sum(res$sig))
      if (nrow(res) && any(res$sig))
        sig_acc[[length(sig_acc) + 1L]] <- res[get("sig") == TRUE,
          c("chr", "lo", "hi", "s_R", "q_R"), with = FALSE][, c("tau", "l_min") := list(tau, lm)]
    }
    if (verbose) { cat(sprintf("   tau = %.2f done\n", tau)); utils::flush.console() }
  }

  landscape <- data.table::rbindlist(cells)
  n_usable <- landscape[get("n_sig") > 0, .N]
  n_cells <- nrow(landscape)
  if (verbose) cat(sprintf("[ld_region_c2] usable cells |U| = %d of %d\n", n_usable, n_cells))

  ## Anchor regions, annotated by the fraction of usable cells in which a
  ## significant region overlaps their span.
  anchor <- .region_table(null$C_obs, anchor_tau, anchor_lmin, edges, co)
  data.table::setnames(anchor, "score", "s_R")
  anchor[, "c2" := NA_real_]
  if (nrow(anchor)) {
    sig_all <- if (length(sig_acc)) data.table::rbindlist(sig_acc) else NULL
    if (is.null(sig_all) || !n_usable) {
      anchor[, "c2" := 0]
    } else {
      for (i in seq_len(nrow(anchor))) {
        ov <- sig_all[get("chr") == anchor$chr[i] &
                      get("lo") <= anchor$hi[i] & get("hi") >= anchor$lo[i]]
        data.table::set(anchor, i, "c2",
                        data.table::uniqueN(ov[, c("tau", "l_min"), with = FALSE]) / n_usable)
      }
    }
    data.table::setorderv(anchor, c("c2", "s_R"), c(-1L, -1L))
  }
  structure(list(regions = anchor, landscape = landscape, n_usable = n_usable,
                 n_cells = n_cells,
                 params = list(tau_grid = tau_grid, lmin_grid = lmin_grid, fdr = fdr,
                               anchor_tau = anchor_tau, anchor_lmin = anchor_lmin, B = B,
                               basis = null$basis %||% NA_character_,
                               engine = null$engine %||% NA_character_)),
            class = "ld_region_c2")
}

#' @export
print.ld_region_c2 <- function(x, ...) {
  p <- x$params
  cat(sprintf("<ld_region_c2> %d anchor region(s) at (tau = %.2f, l_min = %d) | |U| = %d of %d cells\n",
              nrow(x$regions), p$anchor_tau, p$anchor_lmin, x$n_usable, x$n_cells))
  cat(sprintf("  basis: %s | engine: %s | C2 is a ranking, not a second test\n", p$basis, p$engine))
  if (nrow(x$regions)) print(utils::head(x$regions, 10L))
  invisible(x)
}
