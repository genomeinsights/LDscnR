## Internal plumbing shared by ld_unit_matrix(), ld_outlier_test() and
## ld_outlier_perm(). Not exported; written once here rather than reimplemented
## per function the way the pre-package R_3sp_blocks/00_common.R had to be.

## One row per cluster clearing size_floor: unit_id, Chr, from, to, n_markers,
## members (list of marker names), core_snp. `map` supplies Pos for from/to.
.ld_outlier_units <- function(stage1, map, size_floor) {
  cl <- data.table::as.data.table(stage1$clusters)
  nl <- if ("n_loci" %in% names(cl)) cl$n_loci else cl$n_snps
  keep <- nl >= size_floor
  if (!any(keep)) stop("No clusters clear size_floor = ", size_floor, ".")
  cl <- cl[keep]
  mp <- data.table::as.data.table(map)
  pos_of <- stats::setNames(mp$Pos, mp$marker)
  units <- data.table::data.table(
    unit_id   = seq_len(nrow(cl)),
    Chr       = cl$Chr,
    from      = vapply(cl$members, function(m) min(pos_of[m]), 0),
    to        = vapply(cl$members, function(m) max(pos_of[m]), 0),
    n_markers = nl[keep],
    core_snp  = cl$core_snp
  )
  units$members <- cl$members
  data.table::setorder(units, Chr, from)
  units
}

## Simes: min(n * p_(i) / i) over a unit's member marker p-values.
.simes <- function(p) { p <- sort(p[is.finite(p)]); if (!length(p)) return(NA_real_)
  min(length(p) * p / seq_along(p)) }

## BH-adjusted q for a units table's $p column, added in place.
.bh <- function(units) { units$q <- stats::p.adjust(units$p, method = "BH"); units }
