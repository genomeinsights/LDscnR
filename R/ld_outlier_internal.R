## Internal plumbing shared by ld_unit_matrix(), ld_outlier_test() and
## ld_outlier_perm(). Not exported; written once here rather than reimplemented
## per function the way the pre-package R_3sp_blocks/00_common.R had to be.

## THE PATTERN TO USE ANYWHERE THIS PACKAGE LOOKS UP VALUES BY MARKER NAME, repeated
## across many small groups (once per cluster, once per region, ...): ONE match() over
## everything at once, never a named-vector lookup (setNames(x, names)[some_names])
## called once per group.
##
## setNames(x, names)[k] looks cheap and is deceptive: R does not cache a hash on a
## plain named vector between calls, so x[k] internally re-hashes the FULL `names`
## every single time it is called. Looping that once per cluster -- 1,356 times on the
## 3sp panel -- cost 27.7s (building unit spans) and 14s (Simes aggregation) for work
## that a single match() does in milliseconds, because match() hashes its table ONCE
## and reuses that hash across every element of the query in that one call. Found and
## fixed in three places in this file and ld_outlier_test.R before being extracted here
## as the thing to reach for so a fourth instance does not get written the slow way.
##
## @param marker_names Character vector of markers to look up, in any order,
##   possibly with repeats (e.g. unlist(units$members) across every unit at once).
## @param ref_markers Character vector defining the index space that a values
##   vector (p_obs, map$Pos, ...) is aligned to (e.g. map$marker).
## @return Integer vector, same length/order as `marker_names`: each marker's
##   position in `ref_markers`. Index an aligned values vector with this --
##   `map$Pos[.marker_positions(mk, map$marker)]`, `p_obs[.marker_positions(mk, map$marker)]`
##   -- rather than looking values up by name.
.marker_positions <- function(marker_names, ref_markers) {
  match(marker_names, ref_markers)
}

## One row per cluster clearing size_floor: unit_id, Chr, from, to, n_markers,
## members (list of marker names), core_snp. `map` supplies Pos for from/to.
.ld_outlier_units <- function(stage1, map, size_floor) {
  cl <- data.table::as.data.table(stage1$clusters)
  nl <- if ("n_loci" %in% names(cl)) cl$n_loci else cl$n_snps
  keep <- nl >= size_floor
  if (!any(keep)) stop("No clusters clear size_floor = ", size_floor, ".")
  cl <- cl[keep]
  mp <- data.table::as.data.table(map)
  flat_marker <- unlist(cl$members, use.names = FALSE)
  flat_idx    <- .marker_positions(flat_marker, mp$marker)
  flat_unit   <- rep.int(seq_len(nrow(cl)), lengths(cl$members))
  flat_pos    <- mp$Pos[flat_idx]
  span <- data.table::data.table(unit_id = flat_unit, Pos = flat_pos)[
    , .(from = min(Pos), to = max(Pos)), by = unit_id]
  data.table::setkey(span, unit_id)
  units <- data.table::data.table(
    unit_id   = seq_len(nrow(cl)),
    Chr       = cl$Chr,
    from      = span[.(seq_len(nrow(cl)))]$from,
    to        = span[.(seq_len(nrow(cl)))]$to,
    n_markers = nl[keep],
    core_snp  = cl$core_snp
  )
  units$members <- cl$members
  units$member_idx <- unname(split(flat_idx, flat_unit)[as.character(seq_len(nrow(cl)))])
  data.table::setorder(units, Chr, from)
  units
}

## Simes: min(n * p_(i) / i) over a unit's member marker p-values.
.simes <- function(p) { p <- sort(p[is.finite(p)]); if (!length(p)) return(NA_real_)
  min(length(p) * p / seq_along(p)) }

## BH-adjusted q for a units table's $p column, added in place.
.bh <- function(units) { units$q <- stats::p.adjust(units$p, method = "BH"); units }


## Shared by ld_outlier_test() and ld_outlier_perm(): units table with p/q/significant
## added, NO ASSEMBLY. Extracted because assembly (especially stage2_discovered's
## ld_prune_and_eMLG(), which computes eMLG per flagged cluster) is the expensive part,
## and ld_outlier_perm()'s "units" level never needed it -- assembly cannot change which
## units are significant, only how they are grouped into regions afterwards. Before this
## split, EVERY permutation surrogate paid the full assembly cost even when the caller
## only wanted a discovery count. Found by watching a live run stay slow on the arm
## 00_config.R had documented as "cheap".
## `units` may be a PRE-BUILT table (from .ld_outlier_units(), unchanged across
## surrogates) or NULL to build it fresh -- callers in a loop should always pass the
## cached table; ld_outlier_test()'s single observed call is the only caller that
## still builds it inline, where the cost is paid once and does not matter.
.ld_outlier_tested_units <- function(stage1, map, p_obs, statistic, size_floor, alpha,
                                     units = NULL) {
  if (is.null(units)) units <- .ld_outlier_units(stage1, map, size_floor)
  if (statistic == "simes") {
    ## integer indexing into p_obs (aligned to map$marker by contract), not a
    ## character-name lookup -- see .ld_outlier_units()'s comment for why this
    ## matters. member_idx is precomputed once and carried on `units`, so every
    ## surrogate in a permutation loop reuses it for free.
    units$p <- vapply(units$member_idx, function(ix) .simes(p_obs[ix]), 0)
  } else {
    if (length(p_obs) != nrow(units))
      stop(sprintf("statistic = \"unit\": p_obs has %d values but %d units clear size_floor = %d.",
                   length(p_obs), nrow(units), size_floor))
    units$p <- p_obs
  }
  units <- .bh(units)
  units$significant <- !is.na(units$q) & units$q <= alpha
  units
}


## Physical merge: group rows within `gap` bp on the same chromosome. USE THIS, not a
## hand-rolled cumsum -- cummax(to)[-.N] computed over the WHOLE sorted table (not per
## chromosome) never resets at a chromosome boundary. Chr[-1] != Chr[-.N] being TRUE
## still forces a break exactly AT the boundary, but every row AFTER that inherits the
## previous chromosome's stale, larger cummax(to) as its floor -- and because every
## chromosome's own positions restart near zero, that stale value is almost always
## bigger than anything on the new chromosome, making from[-1] - cummax(to)[-.N]
## negative for every row and the gap test never fire again. The result is silent,
## near-total over-merging of every chromosome that sorts after one with any sizeable
## span: measured as one 18.9 Mb "region" built from five real regions whose true
## pairwise gaps were 1.0-9.6 Mb, none under the 300 kb threshold. Found in a
## visualisation script, but the identical pattern was ALSO shipped in this package's
## own ld_outlier_test(assembly = "physical") -- fixed there in the same commit, and
## every number this package has reported so far used assembly = "stage2_discovered"
## (a different, unaffected code path), so nothing already reported is wrong -- but the
## "physical" arm had never been exercised at multi-chromosome scale before this bug
## was found, and its one earlier measurement (17 regions on this panel) should be
## treated as unverified until rerun.
##
## @param D data.table with `Chr`, `from`, `to` (one row per interval to merge).
## @param gap Merge distance in bp.
## @return data.table, one row per merged region: `Chr`, `from`, `to`.
.physical_merge <- function(D, gap) {
  if (!nrow(D)) return(data.table::data.table(Chr = character(), from = numeric(), to = numeric()))
  d <- data.table::copy(D); data.table::setorder(d, Chr, from)
  d[, grp := cumsum(c(TRUE, from[-1] - cummax(to)[-.N] > gap)), by = Chr]
  out <- d[, .(from = min(from), to = max(to)), by = .(Chr, grp)]
  out[, grp := NULL][]
}
