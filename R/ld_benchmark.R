## Outlier-region TP/FP benchmark -- a method-agnostic gold standard for scoring
## any outlier / association scan against SIMULATED ground truth on a common
## outlier-region (OR) frame. Migrated from the LDscnR-paper scoring engine.

#' Flag true-positive QTNs in a simulated map
#'
#' Marks which simulated QTNs count as detectable "true positives": a QTN passing
#' a minor-allele-frequency floor whose per-chromosome share of additive variance
#' (`Va / sum(Va)` over the chromosome's QTNs) is at least `p_va_min`. Sets a
#' logical column `true_pos_QTN` used throughout the benchmark.
#'
#' @param map A data.table with (at least) `type` (`"QTN"` for QTNs), `Chr`, a
#'   MAF column and a per-QTN additive-variance column.
#' @param va_col,maf_col Column names for additive variance and MAF.
#' @param maf_min Minimum MAF for a QTN to be detectable (default 0.1).
#' @param p_va_min Minimum per-chromosome Va share for a QTN to count as a true
#'   positive (default 0.05).
#' @return `map` (copied) with `true_pos_QTN` (and helper `sum_Va_ok`/`p_Va_ok`).
#' @export
flag_true_qtns <- function(map, va_col = "Va", maf_col = "MAF",
                           maf_min = 0.1, p_va_min = 0.05) {
  map <- data.table::copy(data.table::as.data.table(map))
  map[, true_pos_QTN := FALSE]
  qtn_ok <- map$type == "QTN" & map[[maf_col]] > maf_min
  map[qtn_ok, sum_Va_ok := sum(get(va_col)), by = Chr]
  map[qtn_ok, p_Va_ok := get(va_col) / sum_Va_ok]
  map[qtn_ok & p_Va_ok >= p_va_min, true_pos_QTN := TRUE]
  map[]
}

#' LD (r^2) and distance between candidate markers and QTNs
#'
#' For every candidate marker, its r^2 and physical distance to each QTN on the
#' same chromosome within `max_bp`. This is the lookup that decides whether an
#' outlier region is "close enough" to a QTN to count as a detection
#' ([classify_ors()], [evaluate_ors()]).
#'
#' @param GTs Genotype dosage matrix (individuals x SNPs; column names = markers).
#' @param map data.table with `marker`, `Chr`, `Pos`, `type`.
#' @param candidate_markers Markers to tabulate (the OR-membership universe).
#' @param max_bp Maximum candidate-QTN distance to retain (default 2e6).
#' @param cores Cores for the per-chromosome loop.
#' @return A data.table with `marker`, `qtn_marker`, `r2`, `dist_bp`.
#' @export
qtn_ld_table <- function(GTs, map, candidate_markers, max_bp = 2e6, cores = 1) {
  map <- data.table::as.data.table(map)
  chr_levels <- unique(map[marker %in% candidate_markers, Chr])
  out <- parallel::mclapply(chr_levels, function(ch) {
    chr_map      <- map[Chr == ch]
    qtn_markers  <- chr_map[type == "QTN", marker]
    cand_markers <- chr_map[marker %in% candidate_markers, marker]
    if (length(qtn_markers) == 0 || length(cand_markers) == 0) return(NULL)
    gts_qtn <- as.matrix(GTs[, qtn_markers, drop = FALSE]); storage.mode(gts_qtn) <- "double"
    gts_cnd <- as.matrix(GTs[, cand_markers, drop = FALSE]); storage.mode(gts_cnd) <- "double"
    R2 <- stats::cor(gts_cnd, gts_qtn, use = "pairwise.complete.obs")^2
    pos_cnd <- stats::setNames(chr_map$Pos[match(cand_markers, chr_map$marker)], cand_markers)
    pos_qtn <- stats::setNames(chr_map$Pos[match(qtn_markers, chr_map$marker)], qtn_markers)
    dt <- data.table::data.table(
      marker     = rep(cand_markers, times = length(qtn_markers)),
      qtn_marker = rep(qtn_markers, each = length(cand_markers)),
      r2         = as.vector(R2),
      dist_bp    = abs(rep(pos_cnd, times = length(qtn_markers)) -
                         rep(pos_qtn, each = length(cand_markers))))
    dt[dist_bp <= max_bp]
  }, mc.cores = cores)
  data.table::rbindlist(out[!vapply(out, is.null, logical(1))], use.names = TRUE)
}

#' Decay-relative TP-match thresholds
#'
#' The r^2 and distance a region must be within of a QTN to count as a detection,
#' derived from the (median) LD-decay parameters: `r2 = ld_from_rho(b, c, rho_r2)`,
#' `dist = min(d_from_rho(a, rho_d), dmax_cap)`.
#'
#' @param decay_sum Per-chromosome LD-decay summary (`b`, `c`, `a`).
#' @param rho_r2,rho_d Decay-relative r^2 and distance quantiles (defaults
#'   0.75, 0.95).
#' @param dmax_cap Cap on the match distance (default `Inf`).
#' @return A list with `r2min` and `dmax`.
#' @export
score_thresholds <- function(decay_sum, rho_r2 = 0.75, rho_d = 0.95, dmax_cap = Inf) {
  ds <- data.table::as.data.table(decay_sum)
  list(r2min = ld_from_rho(stats::median(ds$b), stats::median(ds$c), rho_r2),
       dmax  = min(d_from_rho(stats::median(ds$a), rho_d), dmax_cap))
}

## internal: assign one OR its focal QTN (highest-r^2 QTN within the thresholds;
## ties broken toward the QTN best supported by the region's neutral markers).
.assign_focal_qtn <- function(cluster_markers, qtn_ld_table, qtn_marker_set, r2_match, d_match) {
  rows <- qtn_ld_table[marker %in% cluster_markers & r2 > r2_match & dist_bp < d_match]
  if (nrow(rows) == 0) return(list(qtn = NA_character_, evidence = NA_real_))
  best_per_qtn <- rows[, .(max_r2 = max(r2)), by = qtn_marker]
  if (nrow(best_per_qtn) == 1) return(list(qtn = best_per_qtn$qtn_marker[1], evidence = best_per_qtn$max_r2[1]))
  neutral_markers <- setdiff(cluster_markers, qtn_marker_set)
  if (length(neutral_markers) > 0) {
    rows_neutral <- rows[marker %in% neutral_markers]
    if (nrow(rows_neutral) > 0) {
      best_neutral <- rows_neutral[, .(max_r2 = max(r2)), by = qtn_marker]
      data.table::setorder(best_neutral, -max_r2)
      winner <- best_neutral$qtn_marker[1]
      return(list(qtn = winner, evidence = best_per_qtn[qtn_marker == winner, max_r2]))
    }
  }
  data.table::setorder(best_per_qtn, -max_r2)
  list(qtn = best_per_qtn$qtn_marker[1], evidence = best_per_qtn$max_r2[1])
}

#' Classify each outlier region as true or false positive
#'
#' Assigns every region its focal QTN (nearest in r^2 within the match
#' thresholds), resolves duplicates (each QTN can be claimed by only one region --
#' the best-supported one), and flags whether the claimed QTN is a true positive.
#' One row per region -- the table an OR-level Manhattan plot needs.
#'
#' @param regions List of member-marker vectors (e.g. `ld_outlier_regions()$regions`).
#' @param map data.table from [flag_true_qtns()] (`type`, `true_pos_QTN`, `marker`).
#' @param qtn_ld_table Output of [qtn_ld_table()].
#' @param r2_match,d_match Match thresholds ([score_thresholds()]).
#' @return A data.table: `CL_id`, `n_loci`, `qtn` (claimed QTN or `NA`),
#'   `evidence` (r^2), `is_TP`.
#' @export
classify_ors <- function(regions, map, qtn_ld_table, r2_match, d_match) {
  map <- data.table::as.data.table(map)
  qtn_marker_set   <- map[type == "QTN", marker]
  true_pos_markers <- map[true_pos_QTN == TRUE, marker]
  if (length(regions) == 0)
    return(data.table::data.table(CL_id = integer(0), n_loci = integer(0),
                                  qtn = character(0), evidence = numeric(0), is_TP = logical(0)))
  assign_list <- lapply(regions, .assign_focal_qtn, qtn_ld_table = qtn_ld_table,
                        qtn_marker_set = qtn_marker_set, r2_match = r2_match, d_match = d_match)
  assignments <- data.table::data.table(
    CL_id    = seq_along(regions),
    n_loci   = vapply(regions, length, integer(1)),
    qtn      = vapply(assign_list, function(x) x$qtn, character(1)),
    evidence = vapply(assign_list, function(x) x$evidence, numeric(1)))
  assigned <- assignments[!is.na(qtn)]
  if (nrow(assigned) > 0) {
    data.table::setorder(assigned, qtn, -evidence)
    keep_id <- assigned[, .SD[1], by = qtn]$CL_id
    drop_id <- setdiff(assigned$CL_id, keep_id)
    assignments[CL_id %in% drop_id, qtn := NA_character_]
  }
  assignments[, is_TP := !is.na(qtn) & qtn %in% true_pos_markers]
  assignments[]
}

## internal: like classify_ors() but retains, for every region, the pre-dedup
## candidate QTN, whether it is itself a true positive, and which region beat it.
## evaluate_ors() uses `dropped_by_dedup` + `candidate_qtn_is_true_positive` to
## identify dedup losers that point at a real QTN (scored as neither TP nor FP).
.diagnose_ors <- function(regions, map, qtn_ld_table, r2_match, d_match) {
  map <- data.table::as.data.table(map)
  qtn_marker_set  <- map[type == "QTN", marker]
  true_pos_lookup <- stats::setNames(map$true_pos_QTN, map$marker)
  pos_lookup      <- stats::setNames(map$Pos, map$marker)
  if (length(regions) == 0) return(data.table::data.table())
  assign_list <- lapply(regions, .assign_focal_qtn, qtn_ld_table = qtn_ld_table,
                        qtn_marker_set = qtn_marker_set, r2_match = r2_match, d_match = d_match)
  dt <- data.table::data.table(
    CL_id         = seq_along(regions),
    n_loci        = vapply(regions, length, integer(1)),
    cluster_pos   = vapply(regions, function(cl) stats::median(pos_lookup[cl], na.rm = TRUE), numeric(1)),
    candidate_qtn = vapply(assign_list, function(x) x$qtn, character(1)),
    evidence      = vapply(assign_list, function(x) x$evidence, numeric(1)))
  dt[, candidate_qtn_is_true_positive := true_pos_lookup[candidate_qtn]]
  dt[, final_qtn := candidate_qtn]
  dt[, dropped_by_dedup := FALSE]
  assigned <- dt[!is.na(candidate_qtn)]
  if (nrow(assigned) > 0) {
    data.table::setorder(assigned, candidate_qtn, -evidence)
    keep_id <- assigned[, .SD[1], by = candidate_qtn]$CL_id
    drop_id <- setdiff(assigned$CL_id, keep_id)
    dt[CL_id %in% drop_id, `:=`(final_qtn = NA_character_, dropped_by_dedup = TRUE)]
  }
  dt[, is_TP := !is.na(final_qtn) & true_pos_lookup[final_qtn]]
  data.table::setorder(dt, CL_id)
  dt[]
}

#' Evaluate outlier regions against ground truth (dedup-neutral)
#'
#' Scores a set of outlier regions into TP / FP / recall / precision on the QTN
#' truth. **Dedup-neutral**: a region that matches an already-claimed
#' true-positive QTN (a duplicate for the same QTN, from over-fragmentation) is
#' counted as *neither* TP nor FP -- it is dropped. This removes the
#' clustering-parameter sensitivity of the metric (splitting one true region into
#' N pieces still yields 1 TP and 0 extra FP), while genuine false regions still
#' count and over-merging is still penalised through recall (one region can claim
#' only one QTN). This makes it a fair, method-agnostic comparison of outlier
#' callers on a common OR frame.
#'
#' @param regions List of member-marker vectors.
#' @param map data.table from [flag_true_qtns()].
#' @param qtn_ld_table Output of [qtn_ld_table()].
#' @param r2_match,d_match Match thresholds ([score_thresholds()]).
#' @return A list: `TP`, `FP`, `FN`, `extras` (dedup losers, uncounted),
#'   `Precision`, `Recall`, `PR` (= Precision x Recall).
#' @seealso [classify_ors()], [pr_auc()], [score_thresholds()]
#' @export
evaluate_ors <- function(regions, map, qtn_ld_table, r2_match, d_match) {
  map <- data.table::as.data.table(map)
  n_true <- map[true_pos_QTN == TRUE, .N]
  if (!length(regions))
    return(list(TP = 0L, FP = 0L, FN = n_true, extras = 0L, Precision = NA_real_, Recall = 0, PR = 0))
  dg <- .diagnose_ors(regions, map, qtn_ld_table, r2_match, d_match)
  extra <- dg$dropped_by_dedup == TRUE & dg$candidate_qtn_is_true_positive %in% TRUE
  TP <- sum(dg$is_TP)
  FP <- sum(!dg$is_TP & !extra)
  FN <- n_true - TP
  Precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  Recall    <- if (n_true > 0)   TP / n_true else 0
  PR        <- if (!is.na(Precision) && !is.na(Recall)) Precision * Recall else NA_real_
  list(TP = TP, FP = FP, FN = FN, extras = sum(extra),
       Precision = Precision, Recall = Recall, PR = PR)
}

#' Trapezoidal precision-recall AUC for an ordered threshold sweep
#'
#' Integrates precision over recall (from recall 0) for the (recall, precision)
#' operating points traced by a single *ordered* knob -- the C-score `tau_C` or a
#' single-SNP `alpha`. Ties in recall collapse to the maximum precision (upper
#' envelope), which absorbs the mild non-monotonicity from re-clustering at each
#' threshold. This is the standard PR-AUC the C-score enables (one ordered knob).
#'
#' @param recall,precision Equal-length numeric vectors of operating points.
#' @return The trapezoidal PR-AUC (numeric scalar), or `NA` if empty.
#' @export
pr_auc <- function(recall, precision) {
  ok <- is.finite(recall) & is.finite(precision)
  if (sum(ok) < 1L) return(NA_real_)
  d <- data.table::data.table(r = recall[ok], p = precision[ok])[, .(p = max(p)), by = r]
  data.table::setorder(d, r)
  r <- c(0, d$r); p <- c(d$p[1], d$p)
  sum(diff(r) * (utils::head(p, -1) + utils::tail(p, -1)) / 2)
}

utils::globalVariables(c(
  "true_pos_QTN", "type", "Chr", "sum_Va_ok", "p_Va_ok", "marker", "Pos",
  "r2", "dist_bp", "qtn_marker", "max_r2", "qtn", "evidence", "CL_id", "is_TP",
  "candidate_qtn", "candidate_qtn_is_true_positive", "final_qtn", "dropped_by_dedup",
  ".SD", "r", "p", "n_loci", "cluster_pos"))
