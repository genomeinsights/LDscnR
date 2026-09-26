#' Test LD-complexity-reduction clusters and assemble significant ones into regions
#'
#' The core detection step of the current Stage-1-unit outlier workflow.
#' It tests units and assembles significant ones into reported regions.
#' It does not construct a null: [ld_outlier_perm()] compares discovery counts
#' with user-supplied surrogate phenotypes, whereas [ld_region_rotation()]
#' asks whether reported intervals overlap external annotations more than
#' their span and chromosome assignment would predict.
#'
#' `statistic = "simes"` aggregates marker-level `p_obs` (`min(n * p_(i) / i)`
#' per cluster, so the penalty scales with cluster size); `statistic = "unit"`
#' takes `p_obs` already computed per unit -- built with [ld_unit_matrix()] and
#' scanned with your own engine -- and does no aggregation, only BH and assembly.
#'
#' `min_r2_rho` for `assembly = "stage2_discovered"` is read from
#' `stage1$params$rho`, NOT a separate argument: the whole point of this
#' assembly rule is that it introduces no parameter beyond what built `stage1`,
#' and a caller-supplied value could silently disagree with that. `stage1` must
#' therefore come from an `ld_complexity_reduction()` new enough to record
#' `params` (added alongside this function); an older `stage1` object errors
#' rather than guessing.
#'
#' @param stage1 An `ld_complexity_reduction` object (must have `$params$rho`).
#' @param map data.frame/data.table with `marker`, `Chr`, `Pos`, aligned to
#'   `stage1`'s marker universe.
#' @param p_obs Per-marker p-values (if `statistic = "simes"`), aligned to
#'   `map$marker`, or per-unit p-values in unit order (if `statistic = "unit"`,
#'   aligned to the units [ld_unit_matrix()] would return at the same
#'   `size_floor`).
#' @param statistic `"simes"` or `"unit"`.
#' @param size_floor Minimum markers per tested unit (default 8L). This
#'   software default is not a universal support threshold; select it without
#'   reference to association results and assess sensitivity to the choice.
#' @param alpha BH level (default 0.05).
#' @param assembly `"stage2_discovered"` (default) -- re-run
#'   [ld_prune_and_eMLG()] over the significant clusters only, so a region can
#'   only span discovered signal -- or `"physical"`, a plain merge of
#'   significant clusters within `gap` on a chromosome. The two are not
#'   equivalent in general; `"physical"` is retained as a near-free check
#'   (agreement between an arbitrary rule and a motivated one is evidence
#'   neither is an artefact of the other), not as an equally-preferred
#'   alternative.
#' @param GTs,LD_decay Required only when `assembly = "stage2_discovered"`.
#' @param score_threshold,distance_threshold [ld_prune_and_eMLG()] parameters
#'   for `"stage2_discovered"`. `ld_w_col`/`ld_w_threshold` are fixed
#'   internally, not exposed: used only to merge already-significant
#'   clusters it is not a gate, only a speed filter, and this call flags
#'   every significant cluster regardless (`ld_w_col = "n_loci"`,
#'   `ld_w_threshold = 0`, and `min_n_loci_flag = 1` -- `n_loci` is always
#'   present and always positive, so this needs nothing beyond what
#'   [ld_complexity_reduction()] itself always produces; it deliberately
#'   does not require an `ld_w_*` column to exist on `stage1$map_snp`).
#' @param gap Physical merge distance in bp, used only when
#'   `assembly = "physical"`.
#'
#' @return An `ld_outlier_test` object: `units` (data.table, one row per tested
#'   unit -- `unit_id`, `Chr`, `from`, `to`, `n_markers`, `p`, `q`,
#'   `significant`), `regions` (data.table, one row per assembled region --
#'   `Chr`, `from`, `to`, `n_units`, `n_markers`, `occupancy`), `params`
#'   (everything above, resolved).
#'
#' @seealso [ld_unit_matrix()], [ld_outlier_perm()], [ld_region_rotation()],
#'   [ld_complexity_reduction()]
#' @export
ld_outlier_test <- function(stage1, map, p_obs,
                            statistic = c("simes", "unit"),
                            size_floor = 8L, alpha = 0.05,
                            assembly = c("stage2_discovered", "physical"),
                            GTs = NULL, LD_decay = NULL,
                            score_threshold = 0.80, distance_threshold = 1e5,
                            gap = 3e5) {
  statistic <- match.arg(statistic)
  assembly <- match.arg(assembly)
  if (assembly == "stage2_discovered" && is.null(stage1$params$rho))
    stop("stage1$params$rho is missing -- this stage1 object predates ld_complexity_reduction() ",
         "recording its own rho. Rebuild it, or pass assembly = \"physical\".")
  if (assembly == "stage2_discovered" && (is.null(GTs) || is.null(LD_decay)))
    stop("assembly = \"stage2_discovered\" needs both GTs and LD_decay.")

  ## per-unit p/q/significant, shared with ld_outlier_perm()'s "units"-level fast path --
  ## see .ld_outlier_tested_units()'s own comment for why this is split out.
  units <- .ld_outlier_tested_units(stage1, map, p_obs, statistic, size_floor, alpha)
  sig <- units[significant == TRUE]

  ## ---- assembly -----------------------------------------------------------
  if (!nrow(sig)) {
    regions <- data.table::data.table(Chr = character(), from = numeric(), to = numeric(),
                                      n_units = integer(), n_markers = integer(),
                                      occupancy = numeric())
  } else if (assembly == "physical") {
    regions <- .physical_merge(sig[, .(Chr, from, to)], gap)
    ## n_units/n_markers per merged region: join the significant units back against the
    ## now-correctly-merged spans (foverlaps by.x, so the region's OWN from/to survive
    ## in the result, prefixed nothing -- the unit's original from/to get the "i." prefix
    ## and are not kept) -- see .physical_merge()'s comment for why this used to be
    ## computed inside the same buggy cumsum step instead.
    data.table::setkey(regions, Chr, from, to)
    ov <- data.table::foverlaps(sig[, .(Chr, from, to, n_markers)], regions,
                                by.x = c("Chr","from","to"), type = "within")
    agg <- ov[, .(n_units = .N, n_markers = sum(n_markers)), by = .(Chr, from, to)]
    regions <- merge(regions, agg, by = c("Chr","from","to"), all.x = TRUE)
    ## occupancy: markers belonging to the region's units / markers physically in its span
    mp <- data.table::as.data.table(map)
    regions[, occupancy := vapply(seq_len(.N), function(i) {
      n_span <- sum(mp$Chr == Chr[i] & mp$Pos >= from[i] & mp$Pos <= to[i])
      n_markers[i] / max(n_span, 1L) }, 0)]
  } else {
    ## stage2_discovered: re-run ld_prune_and_eMLG() over the SIGNIFICANT clusters only.
    ## min_r2_rho comes from stage1$params$rho, never from a caller-supplied value.
    cl <- data.table::as.data.table(stage1$clusters)
    nl <- if ("n_loci" %in% names(cl)) cl$n_loci else cl$n_snps
    cl_sig <- cl[nl >= size_floor][sig$unit_id]
    mk_sig <- unlist(cl_sig$members, use.names = FALSE)
    missing_gt <- setdiff(mk_sig, colnames(GTs))
    if (length(missing_gt))
      stop(sprintf("%d marker(s) in the significant clusters are missing from ",
                   length(missing_gt)), "`colnames(GTs)` (e.g. ",
           paste(utils::head(missing_gt, 5), collapse = ", "), ").")
    ms_sig <- data.table::as.data.table(stage1$map_snp)[marker %chin% mk_sig]
    sub <- structure(list(map_snp = ms_sig, clusters = cl_sig, pruned = cl_sig$core_snp),
                     class = "ld_complexity_reduction")
    ## ld_w_col/ld_w_threshold only need to flag EVERY significant cluster for
    ## the merge pass (see @details/@param ld_w_threshold above) -- with
    ## min_n_loci_flag = 1 already guaranteeing that on its own (every
    ## cluster clears n_snps >= 1), the ld_w_col > ld_w_threshold condition's
    ## own contribution to the flagged set is provably redundant, not a
    ## second, independent gate. Hardcoding ld_w_col = "ld_w_095" here
    ## therefore required a column this function has no other reason to need
    ## -- ld_complexity_reduction() only computes ld_w columns when asked to
    ## separately, so a `stage1` built without that step (a fully valid,
    ## documented ld_outlier_test() input) errored inside ld_prune_and_eMLG()
    ## on a missing column, for a value that was never going to change the
    ## outcome. Using "n_loci" instead -- always present on map_snp, always a
    ## positive integer -- keeps ld_w_threshold = 0 exactly as much of a
    ## structural no-op (n_loci > 0 is guaranteed true, not just empirically
    ## true of real ld_w values) while requiring nothing beyond what
    ## ld_complexity_reduction() itself always produces.
    pr <- ld_prune_and_eMLG(GTs = GTs[, mk_sig, drop = FALSE], stage1 = sub,
                            ld_w_col = "n_loci", ld_w_threshold = 0,
                            LD_decay = LD_decay, min_r2_rho = stage1$params$rho,
                            score_threshold = score_threshold,
                            distance_threshold = distance_threshold,
                            compute_unflagged_eMLG = FALSE, min_n_loci_eMLG = 1,
                            min_n_loci_flag = 1, cores = 1)
    g <- data.table::as.data.table(pr$groups)
    mp <- data.table::as.data.table(map)
    ## .marker_positions(), not a per-region setNames()[...] lookup -- same fix as
    ## .ld_outlier_units() and the Simes aggregation, applied here too so every
    ## marker-name lookup in this package goes through the one fast path.
    g_marker <- unlist(g$members, use.names = FALSE)
    g_idx    <- .marker_positions(g_marker, mp$marker)
    g_group  <- rep.int(seq_len(nrow(g)), lengths(g$members))
    gs <- data.table::data.table(gid = g_group, Chr = as.character(mp$Chr[g_idx]),
                                 Pos = mp$Pos[g_idx])[
      , .(Chr = Chr[1], from = min(Pos), to = max(Pos), n_markers = .N), by = gid]
    data.table::setkey(gs, gid)
    regions <- gs[.(seq_len(nrow(g)))][, .(Chr, from, to, n_units = NA_integer_, n_markers)]
    ## n_units per region: how many of OUR significant units fall inside [from, to]
    regions[, n_units := vapply(seq_len(.N), function(i)
      sum(sig$Chr == Chr[i] & sig$from >= from[i] & sig$to <= to[i]), 0L)]
    regions[, occupancy := vapply(seq_len(.N), function(i) {
      n_span <- sum(mp$Chr == Chr[i] & mp$Pos >= from[i] & mp$Pos <= to[i])
      n_markers[i] / max(n_span, 1L) }, 0)]
  }
  data.table::setorder(regions, Chr, from)

  structure(list(
    units = units[, .(unit_id, Chr, from, to, n_markers, p, q, significant)],
    regions = regions,
    params = list(statistic = statistic, size_floor = size_floor, alpha = alpha,
                  assembly = assembly, score_threshold = score_threshold,
                  distance_threshold = distance_threshold, gap = gap,
                  min_r2_rho = if (assembly == "stage2_discovered") stage1$params$rho else NA)
  ), class = "ld_outlier_test")
}

#' @export
print.ld_outlier_test <- function(x, ...) {
  p <- x$params
  cat(sprintf("<ld_outlier_test> %s / %s | floor = %d, alpha = %.2f | %d/%d units -> %d regions\n",
              p$statistic, p$assembly, p$size_floor, p$alpha,
              sum(x$units$significant), nrow(x$units), nrow(x$regions)))
  invisible(x)
}
