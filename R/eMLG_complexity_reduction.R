# Shared genotype/eMLG helpers -----------------------------------------------
#
# Low-level building blocks shared by make_eMLGs() below, and by
# dynamic_cut_eMLG()/ld_prune_and_eMLG() (see those files) -- consensus
# genotype construction, polarization, and the score_eMLG/weighted_row_mean
# primitives used throughout eMLG-based complexity reduction.

default_cluster_colours <- function() {
  c(
    "#B2DF8A", "#FFD92F", "firebrick", "#33A02C", "#7FC97F", "#CAB2D6",
    "#FB8072", "grey30", "#E6AB02", "#FDC086", "steelblue", "#1F78B4",
    "#FB9A99", "#1B9E77", "#BC80BD", "#E31A1C", "#7570B3", "#A6761D",
    "#A6CEE3", "salmon", "#FFFF33", "forestgreen", "#FDCDAC", "#BF5B17",
    "#A6761D", "#FBB4AE", "#4DAF4A", "#B3E2CD", "#FDDAEC", "#BEBADA",
    "#FFF2AE", "#1F78B4", "#66C2A5", "#F0027F", "#E6AB02", "#E78AC3",
    "#FF7F00", "#8DA0CB", "#6A3D9A", "#B15928", "#E41A1C"
  )
}


detect_gt_input <- function(x, tol = 1e-8) {
  vals <- as.numeric(x)
  vals <- vals[is.finite(vals)]

  if (!length(vals)) return("hard")

  is_hard <- all(abs(vals - round(vals)) < tol) &&
    all(vals %in% c(0, 1, 2))

  if (is_hard) "hard" else "dosage"
}


## Flips (2 - x) any column negatively correlated with the best-covered
## reference column, so downstream row-means/consensus don't cancel out
## alleles that are simply coded on opposite strands/reference alleles.
## Identity for ncol <= 1 (nothing to polarize against).
polarize_genotypes <- function(x, cor_threshold = 0) {
  x <- as.matrix(x)

  if (ncol(x) <= 1L) return(x)

  ref_idx <- which.max(colSums(!is.na(x)))
  ref <- x[, ref_idx]

  flip <- vapply(seq_len(ncol(x)), function(j) {
    r <- suppressWarnings(
      stats::cor(ref, x[, j], use = "pairwise.complete.obs")
    )
    is.finite(r) && r < cor_threshold
  }, logical(1))

  if (any(flip)) {
    x[, flip] <- 2 - x[, flip]
  }

  x
}


expected_gt_hard <- function(x) {
  x <- polarize_genotypes(x)

  n <- rowSums(!is.na(x))
  c1 <- rowSums(x == 1, na.rm = TRUE)
  c2 <- rowSums(x == 2, na.rm = TRUE)

  E <- (c1 + 2 * c2) / n
  E[n == 0] <- NA_real_

  E
}


## The consensus dosage genotype for a cluster: mean dosage across its
## (polarized) member markers per individual. Identity for a single-column
## input (polarize_genotypes() short-circuits, rowMeans() over one column
## returns that column unchanged).
expected_gt_dosage <- function(x) {
  x <- polarize_genotypes(x)

  E <- rowMeans(x, na.rm = TRUE)
  E[is.nan(E)] <- NA_real_

  E
}


## Does a consensus genotype survive hard-calling? cor(round(x), x)^2 --
## the round-trip fidelity downstream LD/Ohta statistics need, since they
## consume hard-called genotypes, not the continuous consensus signal.
score_eMLG <- function(x) {
  r2 <- suppressWarnings(
    stats::cor(round(x), x, use = "pairwise.complete.obs")^2
  )

  if (!is.finite(r2)) NA_real_ else r2
}


## Per-individual weighted average across cluster/group consensus signals,
## weighted by each side's supporting evidence (n_loci) so a single-marker
## cluster doesn't get the same say as one backed by hundreds of markers.
weighted_row_mean <- function(x, w) {
  x <- as.matrix(x)
  w <- as.numeric(w)

  ok <- !is.na(x)
  num <- rowSums(sweep(x, 2, w, `*`), na.rm = TRUE)
  den <- rowSums(sweep(ok, 2, w, `*`), na.rm = TRUE)

  y <- num / den
  y[den == 0] <- NA_real_

  y
}


cluster_level_map <- function(map_snp, id_col = "CL_id") {
  data.table::setDT(map_snp)

  map_snp[
    !is.na(get(id_col)),
    {
      pos <- as.numeric(Pos)

      .(
        Chr = Chr[1],
        Pos_min = min(pos, na.rm = TRUE),
        Pos_max = max(pos, na.rm = TRUE),
        Pos_mid = median(pos, na.rm = TRUE),
        span_bp = max(pos, na.rm = TRUE) - min(pos, na.rm = TRUE),
        n_loci = .N,
        r2_eMLG = if ("r2_eMLG" %in% names(.SD)) r2_eMLG[1] else NA_real_,
        col = if ("CL_col" %in% names(.SD)) CL_col[1] else NA_character_
      )
    },
    by = id_col
  ]
}

add_snp_positions <- function(map_snp, map_ref) {
  map_ref <- data.table::as.data.table(map_ref)
  map_snp <- data.table::as.data.table(map_snp)

  pos_cols <- intersect(c("marker", "Chr", "Pos"), names(map_ref))

  merge(
    unique(map_ref[, ..pos_cols]),
    map_snp,
    by = "marker",
    all.y = TRUE,
    sort = FALSE
  )
}

# Main function ----------------------------------------------------------

#' Build eMLGs (Expected Multi-Locus Genotypes) From Pre-Defined Clusters
#'
#' Reduces a pre-existing set of LD clusters (e.g. from [ld_complexity_reduction()]
#' or any other clustering that assigns a `CL_id` per marker) directly to one
#' consensus genotype per cluster -- no merging or restructuring of the
#' cluster boundaries themselves. For that, see [dynamic_cut_eMLG()] /
#' [ld_prune_and_eMLG()], which additionally consolidate clusters that are
#' still correlated with each other via a quality-gated dynamic tree cut;
#' `make_eMLGs()` is the simpler building block that just summarizes clusters
#' as they already stand.
#'
#' @param GTs Numeric matrix of genotypes, individuals x markers (hard-called
#'   0/1/2 or continuous dosage -- see `input`).
#' @param map_cl A data.table/data.frame with at least `marker`, `CL_id`, and
#'   `n_loci` columns (rows with `NA` `CL_id` are dropped).
#' @param input `"auto"` (default) detects hard-called vs. dosage genotypes
#'   from the data itself; `"hard"`/`"dosage"` force one or the other.
#' @param cor_th Stored on the returned object for downstream reference; not
#'   used to filter clusters here (compare [dynamic_cut_eMLG()]'s
#'   `threshold`, which does gate on it).
#' @param l_min Minimum cluster size (raw marker count) to keep; smaller
#'   clusters are dropped entirely.
#' @param ncores Number of cores for the per-cluster consensus computation
#'   (`parallel::mclapply()`, fork-based -- not available on Windows).
#' @param col_vector Optional vector of colours to cycle across clusters for
#'   plotting; defaults to [default_cluster_colours()].
#'
#' @return A list: `eMLG` (matrix, individuals x clusters), `map_eMLG`
#'   (one row per cluster: position span, `n_loci`, `r2_eMLG`), `map_SNP`
#'   (one row per marker, with its cluster's `r2_eMLG`/colour attached),
#'   `clusters` (named list of marker vectors, one per cluster), `input`,
#'   `cor_th`.
#'
#' @examples
#' \dontrun{
#' eMLGs <- make_eMLGs(
#'   GTs, map_cl[!is.na(CL_id) & n_loci >= 10, ],
#'   cor_th = 0.8, l_min = 1, ncores = 4
#' )
#' }
#'
#' @export
make_eMLGs <- function(GTs,
                       map_cl,
                       input = c("auto", "hard", "dosage"),
                       cor_th = 0.8,
                       l_min = 10,
                       ncores = 1,
                       col_vector = NULL) {
  input <- match.arg(input)

  if (is.null(col_vector)) {
    col_vector <- default_cluster_colours()
  }

  map_cl <- data.table::copy(map_cl)
  data.table::setDT(map_cl)

  map_cl <- map_cl[!is.na(CL_id) & n_loci >= l_min]

  cls <- split(map_cl$marker, map_cl$CL_id)
  cls <- cls[lengths(cls) >= l_min]

  if (!length(cls)) {
    stop("No LD clusters passed l_min.")
  }

  all_markers <- unique(unlist(cls, use.names = FALSE))

  if (input == "auto") {
    input <- detect_gt_input(GTs[, all_markers, drop = FALSE])
    message("Detected genotype input type: ", input)
  }

  expected_fun <- switch(
    input,
    hard = expected_gt_hard,
    dosage = expected_gt_dosage
  )

  process_cluster <- function(markers) {
    expected_fun(GTs[, markers, drop = FALSE])
  }

  res <- parallel::mclapply(cls, process_cluster, mc.cores = ncores)

  eMLG <- do.call(cbind, res)
  colnames(eMLG) <- names(cls)
  rownames(eMLG) <- rownames(GTs)

  r2_eMLG <- apply(eMLG, 2, score_eMLG)

  cl_cols <- rep(col_vector, length.out = length(cls))

  map_snp <- data.table::data.table(
    CL_id = rep(names(cls), lengths(cls)),
    marker = unlist(cls, use.names = FALSE),
    CL_col = rep(cl_cols, lengths(cls)),
    n_loci = rep(lengths(cls), lengths(cls)),
    r2_eMLG = rep(r2_eMLG, lengths(cls))
  )

  map_snp <- add_snp_positions(map_snp, map_cl)

  map_eMLG <- cluster_level_map(
    map_snp,
    id_col = "CL_id"
  )

  map_eMLG <- map_eMLG[match(colnames(eMLG), CL_id)]


  list(
    eMLG = eMLG,
    map_eMLG = map_eMLG,
    map_SNP = map_snp,
    clusters = cls,
    input = input,
    cor_th = cor_th
  )
}
