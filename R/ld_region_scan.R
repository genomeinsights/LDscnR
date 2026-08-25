## Internal helpers shared by the gate, the region scan and the second-tier score.

## Marker -> position and marker -> chromosome, read straight off an `ld_edges`
## object so downstream functions never need the map again.
.edge_coords <- function(edges) {
  mk <- unlist(lapply(edges, `[[`, "marker"), use.names = FALSE)
  ps <- unlist(lapply(edges, `[[`, "pos"), use.names = FALSE)
  ch <- rep(names(edges), vapply(edges, function(E) length(E$marker), integer(1)))
  list(pos = stats::setNames(ps, mk), chr = stats::setNames(ch, mk))
}

.empty_regions <- function()
  data.table::data.table(chr = character(), lo = numeric(), hi = numeric(),
                         size = integer(), score = numeric())

## Cluster the markers of one C-vector at threshold `tau` into LD regions, keep
## those with at least `l_min` markers, and summarise each by its span and its
## summed C-mass (the region score s_R).
.region_table <- function(C, tau, l_min, edges, co) {
  mk <- names(C)[which(C >= tau)]                     # which() also drops NAs
  if (!length(mk)) return(.empty_regions())
  r <- ld_regions(mk, edges)
  r <- r[lengths(r) >= l_min]
  if (!length(r)) return(.empty_regions())
  data.table::data.table(
    chr   = unname(co$chr[vapply(r, `[`, character(1), 1L)]),
    lo    = vapply(r, function(x) min(co$pos[x]), numeric(1)),
    hi    = vapply(r, function(x) max(co$pos[x]), numeric(1)),
    size  = lengths(r),
    score = vapply(r, function(x) sum(C[x]), numeric(1)))
}

## Location-matched empirical p and BH q for a set of observed regions `O`
## against the surrogate regions `S` (one table, with a `b` surrogate column).
## Shared by ld_region_scan() and ld_region_c2() so the two cannot drift apart.
.region_pq <- function(O, S, B, fdr) {
  if (!nrow(O)) return(cbind(O, data.table::data.table(
    id = integer(), n_null_ge = integer(), p_R = numeric(), q_R = numeric(), sig = logical())))
  O <- data.table::copy(O)
  O[, "id" := seq_len(.N)]
  O[, "n_null_ge" := 0L]
  if (nrow(S)) {
    ## For each observed region, the strongest surrogate region overlapping its
    ## span, per surrogate; then count the surrogates reaching s_R.
    data.table::setkeyv(S, c("chr", "lo", "hi"))
    ov <- data.table::foverlaps(O[, c("id", "s_R", "chr", "lo", "hi"), with = FALSE], S,
                                by.x = c("chr", "lo", "hi"), type = "any", nomatch = NULL)
    if (nrow(ov)) {
      best <- ov[, list(bs = max(get("score"))), by = c("id", "s_R", "b")]
      nh <- best[get("bs") >= get("s_R"), list(n = data.table::uniqueN(get("b"))), by = "id"]
      if (nrow(nh)) O[nh, "n_null_ge" := get("i.n"), on = "id"]
    }
  }
  O[, "p_R" := (1 + get("n_null_ge")) / (1 + B)]
  O[, "q_R" := stats::p.adjust(get("p_R"), "BH")]
  O[, "sig" := get("q_R") < fdr]
  O[]
}

## All surrogate regions at one (tau, l_min), stacked with a `b` column.
.surr_table <- function(null, tau, l_min, edges, co) {
  data.table::rbindlist(lapply(seq_len(null$B), function(b) {
    ## `b` is set unconditionally: a surrogate with no regions must still
    ## contribute a 0-row table of the SAME shape, or rbindlist() errors.
    .region_table(null$C_surr[[b]], tau, l_min, edges, co)[, "b" := b]
  }))
}

#' The background gate -- check a null before reading any p-value from it
#'
#' For each null bundle, the **median per-surrogate** count of markers with
#' `C > 0` and of LD-regions of at least `l_min` markers, alongside the observed
#' counts. A basis whose surrogates approach the observed counts is mis-specified:
#' it reproduces the very signal it is supposed to be a null for, and its
#' p-values are not interpretable. Run this *before* [ld_region_scan()], and
#' report it -- a mis-specified null that is disclosed is a result, whereas one
#' that is quietly dropped is not.
#'
#' Report the median per surrogate, never a count pooled across surrogates: a
#' pooled total confounds "every surrogate is mildly noisy" with "one surrogate
#' exploded", and the two have opposite implications. For the same reason the
#' returned table carries `q3_regions`, `max_regions` and `frac_surr_any_region`
#' next to the median -- permutation nulls are often heavy-tailed, silent in four
#' surrogates out of five and explosive in the fifth, and a median alone reports
#' that as "silent".
#'
#' @param null An `ld_null` object ([ld_null_from_p()] or [structured_null()]),
#'   or a **named list** of them to tabulate several bases at once.
#' @param edges An [ld_edges()] object. When `null` is a list, build one edge
#'   graph over the union of their `universe`s and pass it here: edge membership
#'   is a property of a marker pair, not of the marker set, so a shared graph is
#'   what makes the rows comparable.
#' @param tau C-score threshold at which regions are formed (default 0.05).
#' @param l_min Minimum region size in markers (default 3).
#' @param warn_at Fraction of the observed region count at which a basis is
#'   flagged as failing the gate (default 0.5).
#'
#' @return A `data.table`, one row per basis, with the observed counts, the
#'   median/IQR/max per-surrogate counts, `frac_surr_any_region`, the
#'   median-to-observed ratio, and a logical `pass`. `pass` is `TRUE` only when
#'   the observed data produced at least one region *and* the median surrogate
#'   count stays below `warn_at` of it; zero observed regions is a failure, not
#'   a pass, and leaves `ratio` as `NA`.
#' @seealso [ld_region_scan()], [ld_null_from_p()]
#' @export
ld_gate <- function(null, edges, tau = 0.05, l_min = 3L, warn_at = 0.5) {
  if (inherits(null, "ld_null")) null <- list(null)
  if (is.null(names(null))) names(null) <- vapply(null, function(x) x$basis %||% "?", character(1))
  co <- .edge_coords(edges)
  n_reg <- function(C) nrow(.region_table(C, tau, l_min, edges, co))

  out <- data.table::rbindlist(lapply(seq_along(null), function(i) {
    x <- null[[i]]
    cg <- vapply(x$C_surr, function(s) sum(s > 0), numeric(1))
    rg <- vapply(x$C_surr, n_reg, numeric(1))
    data.table::data.table(
      basis = names(null)[i], engine = x$engine %||% NA_character_, B = x$B,
      obs_Cgt0 = sum(x$C_obs > 0), med_Cgt0 = stats::median(cg),
      obs_regions = n_reg(x$C_obs), med_regions = stats::median(rg),
      q3_regions = stats::quantile(rg, 0.75, names = FALSE), max_regions = max(rg),
      frac_surr_any_region = mean(rg > 0))
  }))
  out[, "ratio" := ifelse(get("obs_regions") > 0, get("med_regions") / get("obs_regions"), NA_real_)]
  ## A basis passes only if the observed data produced regions AND the median
  ## surrogate stays below `warn_at` of that count. "Observed found nothing,
  ## surrogates found plenty" is the most emphatic gate failure there is, and it
  ## must not be reported as a pass just because the ratio is undefined. `ratio`
  ## is left NA -- it genuinely is undefined at zero observed regions, and the
  ## obs_regions column already tells the reader why `pass` is FALSE.
  out[, "pass" := get("obs_regions") > 0 & !is.na(get("ratio")) & get("ratio") < warn_at]
  zero  <- out$basis[out$obs_regions == 0]
  noisy <- out$basis[out$obs_regions > 0 & !out$pass]
  if (length(zero))
    warning(sprintf("Gate FAILED for: %s -- no observed regions at this (tau, l_min), so there is nothing for the surrogates to be measured against.",
                    paste(zero, collapse = ", ")), call. = FALSE)
  if (length(noisy))
    warning(sprintf("Gate FAILED for: %s -- median surrogate regions reach >= %.0f%% of observed. Do not read p-values from these bases.",
                    paste(noisy, collapse = ", "), 100 * warn_at), call. = FALSE)
  structure(out, class = c("ld_gate", class(out)),
            params = list(tau = tau, l_min = l_min, warn_at = warn_at))
}

#' Call outlier regions and test each against the null at its own locus
#'
#' The one instrument of the LD-scan pipeline. Regions are formed by clustering
#' the markers whose C-score reaches `tau` ([ld_regions()]) and keeping those with
#' at least `l_min` markers. Each observed region `R` is scored by its summed
#' C-mass, and its p-value is the fraction of surrogates that produce a region
#' *overlapping R's own genomic span* with at least as much mass:
#'
#' \deqn{s_R = \sum_{j \in R} C_j, \qquad
#'       p_R = \frac{1 + \#\{b : s^{null}_{b,R} \ge s_R\}}{1 + B}}
#'
#' followed by Benjamini-Hochberg across the discovered regions to give `q_R`.
#'
#' The test is deliberately **location-matched**, and that is what makes it work
#' where a genome-wide count does not. A signal confined to part of the sampled
#' range need not be reproduced by a structure-preserving surrogate *at its own
#' locus*, so a genuinely range-restricted peak -- which a pooled count would
#' charge to structure -- can still separate. It is also robust to a globally
#' inflated null, because excess null peaks are scattered and rarely land on any
#' particular observed locus.
#'
#' `q_R`, not a Bonferroni correction over every marker, is the multiple-testing
#' correction this pipeline reports: the per-marker burden was already carried at
#' the C-score and region-size stages.
#'
#' **Run [ld_gate()] first.** A p-value from a null whose surrogates reproduce the
#' observed background is not interpretable, and nothing in this function can
#' detect that for you.
#'
#' @param null An `ld_null` object ([ld_null_from_p()] or [structured_null()]).
#' @param edges An [ld_edges()] object built over `null$universe`.
#' @param tau C-score threshold for forming regions (default 0.05).
#' @param l_min Minimum region size in markers (default 3).
#' @param fdr BH level at which `sig` is flagged (default 0.05).
#'
#' @return An `ld_region_scan` object: `regions` (a `data.table` with `chr`,
#'   `lo`, `hi`, `size`, `s_R`, `n_null_ge`, `p_R`, `q_R`, `sig`), the `p_floor`
#'   `1/(1+B)`, and `params`.
#' @seealso [ld_gate()], [ld_region_c2()], [ld_scan()], [ld_null_from_p()]
#' @export
ld_region_scan <- function(null, edges, tau = 0.05, l_min = 3L, fdr = 0.05) {
  stopifnot(inherits(null, "ld_null"))
  co <- .edge_coords(edges)
  B <- null$B
  O <- .region_table(null$C_obs, tau, l_min, edges, co)
  data.table::setnames(O, "score", "s_R")
  O <- .region_pq(O, .surr_table(null, tau, l_min, edges, co), B, fdr)
  if (nrow(O)) {
    data.table::setorderv(O, c("q_R", "s_R"), c(1L, -1L))
    data.table::setcolorder(O, c("id", "chr", "lo", "hi", "size", "s_R",
                                 "n_null_ge", "p_R", "q_R", "sig"))
  }
  structure(list(regions = O, p_floor = 1 / (1 + B),
                 params = list(tau = tau, l_min = l_min, fdr = fdr, B = B,
                               basis = null$basis %||% NA_character_,
                               engine = null$engine %||% NA_character_)),
            class = "ld_region_scan")
}

#' @export
print.ld_region_scan <- function(x, ...) {
  p <- x$params
  cat(sprintf("<ld_region_scan> %d region(s), %d significant at q_R < %.2f | tau = %.2f, l_min = %d, B = %d\n",
              nrow(x$regions), sum(x$regions$sig %||% FALSE), p$fdr, p$tau, p$l_min, p$B))
  cat(sprintf("  basis: %s | engine: %s | p-floor = %.4f\n",
              p$basis, p$engine, x$p_floor))
  if (nrow(x$regions)) {
    if (all(x$regions$p_R == x$p_floor))
      cat("  NOTE: every region sits at the p-floor. q_R says 'not a structure artifact';\n",
          "       it cannot rank them. Use ld_region_c2() for that.\n", sep = "")
    print(utils::head(x$regions, 10L))
  }
  invisible(x)
}

#' @export
print.ld_gate <- function(x, ...) {
  p <- attr(x, "params")
  cat(sprintf("<ld_gate> tau = %.2f, l_min = %d, median per surrogate (flag at ratio >= %.2f)\n",
              p$tau, p$l_min, p$warn_at))
  print(data.table::as.data.table(x)[, c("basis", "engine", "B", "obs_regions", "med_regions",
                                         "max_regions", "frac_surr_any_region", "ratio", "pass"),
                                     with = FALSE])
  invisible(x)
}
