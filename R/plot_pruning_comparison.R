#' Plot Stage 1 vs combined (ld_prune_and_eMLG) cluster diagnostic for one chromosome
#'
#' A diagnostic scatter plot (`ld_w` vs. position, points coloured by
#' cluster) stacking [ld_complexity_reduction()]'s raw Stage-1 clusters on
#' top of [ld_prune_and_eMLG()]'s consolidated groups for one chromosome, so
#' fragmented low-recombination blocks (Stage 1, top) and their reunited
#' counterparts (Combined, bottom) can be compared directly. Saves the
#' stacked figure to disk and returns it invisibly.
#'
#' @param chr Chromosome name (e.g. "Chr26"), used in plot labels/filename.
#' @param pruned_stage1 Output of [ld_complexity_reduction()] (whole-genome).
#' @param result Output of [ld_prune_and_eMLG()] (whole-genome).
#' @param map Map data.table with Chr, Pos, marker, and `ld_w_col`.
#' @param ld_w_col Name of the ld_w column in `map`/`pruned_stage1$map_snp`.
#' @param ld_w_threshold Threshold used to flag clusters for the Stage-1
#'   comparison view -- should match what was used to build `result`.
#' @param min_n_loci_flag Must match whatever `min_n_loci_flag` was passed
#'   to the [ld_prune_and_eMLG()] call that produced `result` (default
#'   `Inf`, i.e. off). Without this, the Stage 1 panel (pure `ld_w_col >
#'   ld_w_threshold`) and the Combined panel (`result$groups`' "F"/"U"
#'   prefixes, which also reflect `min_n_loci_flag`) fall out of sync:
#'   clusters pulled into the flagged/merge pathway purely by size, not
#'   ld_w, show up in the Combined "high" panel and vanish from Combined
#'   "low" with no Stage-1-side counterpart -- checked directly on real
#'   data (ld_w_threshold=0.025, min_n_loci_flag=5): ~1.5% of markers landed
#'   in Combined "F" this way, missing from Combined "U", before this
#'   parameter was added.
#' @param direction "high" (default) shows the FLAGGED clusters/groups
#'   (`ld_w_col > threshold` OR `n_snps >= min_n_loci_flag`, group_id prefix
#'   "F") -- the main diagnostic view. "low" shows the UNFLAGGED side (the
#'   complement, group_id prefix "U") as a sanity check: this view should
#'   never contain markers with `ld_w_col > threshold`, since "unflagged" is
#'   defined as no cluster member exceeding threshold AND the cluster not
#'   being pulled in by `min_n_loci_flag`. The "high" view, by contrast, is
#'   expected to contain some sub-threshold "boundary" markers (cluster-
#'   mates of a flagged marker, or whole clusters flagged only via
#'   `min_n_loci_flag`) since flagging happens at the cluster level, not the
#'   marker level.
#' @param out_folder Folder to write the figure to (created if missing).
#' @param width,height Figure dimensions in inches, passed to
#'   `ggplot2::ggsave()`.
#'
#' @return The stacked comparison plot, invisibly (also saved to disk as
#'   `<out_folder>/<chr>_stage1_vs_combined_<direction>.png`).
#'
#' @examples
#' \dontrun{
#' stage1 <- ld_complexity_reduction(map = map, LD_decay = ld_decay, rho = 0.5)
#' result <- ld_prune_and_eMLG(
#'   GTs = GTs, stage1 = stage1, ld_w_col = "ld_w_095", ld_w_threshold = 0.2
#' )
#' plot_pruning_comparison("Chr1", stage1, result, map, ld_w_threshold = 0.2)
#' }
#'
#' @export
plot_pruning_comparison <- function(chr, pruned_stage1, result, map,
                                     ld_w_col = "ld_w_095", ld_w_threshold = 0.2,
                                     min_n_loci_flag = Inf,
                                     direction = c("high", "low"),
                                     out_folder = "./Figures/", width = 10, height = 9) {
  direction <- match.arg(direction)
  dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)

  pal_cluster <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33",
                    "#A65628","#F781BF","#1B9E77","#D95F02","#7570B3","#66A61E")

  ## cycle a modest, visually-distinct palette across CL_id -- with hundreds
  ## of clusters no palette gives every one a unique color, but neighbouring
  ## clusters (what we care about here) will very likely differ
  plot_clusters <- function(map_snp, title) {
    dt <- data.table::copy(map_snp)
    dt[, cl_rank := match(CL_id, unique(CL_id))]
    dt[, col := pal_cluster[(cl_rank - 1) %% length(pal_cluster) + 1]]
    ggplot2::ggplot(dt, ggplot2::aes(Pos / 1e6, .data[[ld_w_col]], color = col)) +
      ggplot2::geom_point(size = 1.2, alpha = 0.85) +
      ggplot2::scale_color_identity() +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::labs(x = paste(chr, "position (Mbp)"), y = expression(ld["w,"*rho*"=0.95"]), title = title)
  }

  ## "high"/"low" pick the cluster set (Stage 1) and the matching group_id
  ## prefix (Combined) together, so the two panels always stay in sync --
  ## see @param direction above. NOTE: "low" must be the COMPLEMENT of the
  ## "high" (any member > threshold, OR min_n_loci_flag) cluster set, not a
  ## separate "any member < threshold" condition -- those are not the same
  ## thing, since a flagged cluster can (and often does) contain
  ## sub-threshold "boundary" members too. Using "any member < threshold"
  ## directly would pull flagged clusters right back into the "low" view.
  ## Mirrors ld_prune_and_eMLG()'s own flagging logic exactly, so this stays
  ## in sync with whatever produced `result` (see @param min_n_loci_flag).
  group_prefix <- if (direction == "high") "F" else "U"
  flagged_ids <- pruned_stage1$map_snp[get(ld_w_col) > ld_w_threshold, unique(CL_id)]
  if (is.finite(min_n_loci_flag)) {
    flagged_ids <- union(flagged_ids, pruned_stage1$clusters[n_snps >= min_n_loci_flag, CL_id])
  }
  chr_ids <- pruned_stage1$map_snp[Chr == chr, unique(CL_id)]
  needs_merge_ids <- if (direction == "high") flagged_ids else setdiff(chr_ids, flagged_ids)
  stage1_snp <- pruned_stage1$map_snp[Chr == chr & CL_id %in% needs_merge_ids]
  message(chr, " Stage 1 (", direction, "): ", data.table::uniqueN(stage1_snp$CL_id), " clusters")
  p_stage1 <- plot_clusters(
    stage1_snp,
    sprintf("Stage 1: ld_complexity_reduction() -- %d clusters (%s)", data.table::uniqueN(stage1_snp$CL_id), direction)
  )

  ## combined result restricted to this chromosome's matching-direction
  ## groups only (group_id prefix "F" for flagged/high, "U" for
  ## unflagged/low) -- matches the Stage 1 scope above; result$groups
  ## otherwise includes both buckets genome-wide, making the comparison
  ## apples-to-oranges
  groups_chr <- result$groups[Chr == chr & startsWith(group_id, group_prefix)]
  final_snp <- groups_chr[, .(marker = unlist(members)), by = group_id]
  data.table::setnames(final_snp, "group_id", "CL_id")
  final_snp <- map[Chr == chr, c("marker", "Pos", ld_w_col), with = FALSE][final_snp, on = "marker"]
  message(chr, " combined (", direction, "): ", data.table::uniqueN(final_snp$CL_id), " groups")
  p_combined <- plot_clusters(
    final_snp,
    sprintf("ld_prune_and_eMLG() -- %d groups (%s)", data.table::uniqueN(final_snp$CL_id), direction)
  )

  p_compare <- p_stage1 / p_combined
  fname <- paste0(out_folder, chr, "_stage1_vs_combined_", direction, ".png")
  ggplot2::ggsave(fname, p_compare, width = width, height = height)
  message("Saved: ", fname)

  invisible(p_compare)
}
