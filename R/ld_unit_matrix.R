#' One variable per LD-complexity-reduction cluster, for your own scan
#'
#' `ld_outlier_test()` tests clusters, not markers, and a cluster is not itself
#' something you can hand to an association test -- it has to become ONE variable
#' per individual first. This builds that variable, three ways, and stops there:
#' running the test on it (EMMAX, LFMM, a GLM, whatever) is your job, exactly as
#' with the per-marker p-values [ld_scan()] expects. The package supplies LD
#' structure and region assembly; it does not run your scan.
#'
#' The three representations answer different questions and are not
#' interchangeable in general, only observed to be close on some panels:
#'
#' - `"consensus_dosage"`: polarise every member to a common allele, then
#'   row-mean. No size penalty -- a larger cluster is a BETTER-ESTIMATED
#'   consensus, so size helps rather than hurts. Averages away signal that is
#'   genuinely SNP-specific (differing even between strongly linked markers).
#' - `"eMLG"`: [make_eMLGs()]'s own block consensus. Conceptually close to
#'   `consensus_dosage` but not verified identical; included so that question is
#'   answered empirically rather than assumed.
#' - `"representative"`: the cluster's core SNP, chosen by [ld_complexity_reduction()]
#'   for CENTRALITY (highest median r^2 to the rest of the cluster) -- not for
#'   signal. A real, single SNP, at whatever centrality picked it.
#'
#' A fourth SNP-level choice, the member most correlated with a block's eMLG
#' consensus, is available via [ld_group_matrix()] instead -- not here,
#' because it operates on [ld_prune_and_eMLG()]'s groups, a different (and
#' possibly coarser, where clusters were merged) partition from this
#' function's Stage-1 units. See [ld_group_matrix()]'s documentation for why
#' its output cannot be substituted for one of this function's columns.
#'
#' @param GTs Genotype dosage matrix (individuals x SNPs, column names = markers).
#' @param stage1 An `ld_complexity_reduction` object (or equivalent list with
#'   `clusters$members`, `clusters$core_snp`).
#' @param map data.frame/data.table with `marker`, `Chr`, `Pos`, aligned to
#'   `stage1`'s marker universe.
#' @param size_floor Minimum markers per cluster to be included (default 8L).
#'   `"eMLG"`'s own default inside [make_eMLGs()] is `l_min = 10`, LARGER than a
#'   typical `size_floor` -- this function passes `l_min = 1` to [make_eMLGs()]
#'   internally so filtering happens once, here, at `size_floor`; override via
#'   `emlg_args$l_min` only if you want a second, stricter gate inside eMLG
#'   itself.
#' @param repr One of `"consensus_dosage"`, `"eMLG"`, `"representative"`.
#' @param emlg_args Extra arguments to [make_eMLGs()] (used only when
#'   `repr = "eMLG"`).
#'
#' @return For `repr %in% c("consensus_dosage", "eMLG")`: an individuals x units
#'   numeric matrix. For `repr = "representative"`: a character vector of marker
#'   names (one per unit; index into `GTs` yourself). Every return carries
#'   `attr(, "units")`, the [ld_outlier_test()]-shaped units table (`unit_id`,
#'   `Chr`, `from`, `to`, `n_markers`), aligned to columns (or to the vector, for
#'   `"representative"`), built from `stage1`.
#'
#' @seealso [ld_outlier_test()], [ld_complexity_reduction()], [make_eMLGs()],
#'   [ld_group_matrix()]
#' @export
ld_unit_matrix <- function(GTs, stage1, map, size_floor = 8L,
                           repr = c("consensus_dosage", "eMLG", "representative"),
                           emlg_args = list()) {
  repr <- match.arg(repr)
  units <- .ld_outlier_units(stage1, map, size_floor)
  units_summary <- units[, .(unit_id, Chr, from, to, n_markers)]
  if (repr %in% c("consensus_dosage", "eMLG")) {
    ## Checked once, here, before either arm's per-unit loop: GTs[, mk]
    ## on a marker absent from colnames(GTs) throws R's own "subscript out
    ## of bounds", which does not say which marker or which unit, and (for
    ## "eMLG") would otherwise fail partway through an expensive
    ## make_eMLGs() call instead of immediately.
    missing_gt <- setdiff(unlist(units$members, use.names = FALSE), colnames(GTs))
    if (length(missing_gt))
      stop(sprintf("%d marker(s) in the tested units are missing from `colnames(GTs)` ",
                   length(missing_gt)), "(e.g. ",
           paste(utils::head(missing_gt, 5), collapse = ", "), ").")
  }

  out <- switch(repr,
    consensus_dosage = {
      m <- vapply(units$members, function(mk) consensus_dosage(GTs, mk),
                  numeric(nrow(GTs)))
      colnames(m) <- units$unit_id
      m
    },
    representative = {
      v <- units$core_snp
      names(v) <- units$unit_id
      v
    },
    eMLG = {
      ## make_eMLGs()'s own cluster_level_map() unconditionally needs Chr and Pos in
      ## map_cl, not just marker/CL_id/n_loci as ld_complexity_reduction()'s own
      ## @return implied -- found by running this against a real bundle, not by
      ## reading the docs alone. Pulled from `map` by marker name.
      mp <- data.table::as.data.table(map)
      map_cl <- data.table::rbindlist(lapply(seq_len(nrow(units)), function(i) {
        mk <- units$members[[i]]
        d <- mp[marker %chin% mk]
        data.table::data.table(marker = d$marker, Chr = d$Chr, Pos = d$Pos,
                               CL_id = units$unit_id[i], n_loci = units$n_markers[i]) }))
      args <- utils::modifyList(list(GTs = GTs, map_cl = map_cl, l_min = 1L), emlg_args)
      res <- do.call(make_eMLGs, args)
      ## make_eMLGs applies its OWN l_min gate; with l_min = 1 that should be a
      ## no-op, but confirm rather than assume -- a stricter emlg_args$l_min
      ## silently drops units size_floor already let through.
      got <- colnames(res$eMLG)
      missing <- setdiff(as.character(units$unit_id), got)
      if (length(missing)) warning(sprintf(
        "make_eMLGs() dropped %d unit(s) that cleared size_floor -- check emlg_args$l_min.",
        length(missing)))
      units_summary <- units_summary[as.character(unit_id) %in% got]
      res$eMLG[, as.character(units_summary$unit_id), drop = FALSE]
    })
  attr(out, "units") <- units_summary
  out
}
