#' Compute Chromosome-wise LD Structure and Decay
#'
#' Estimates linkage disequilibrium (LD) structure and LD-decay parameters
#' across chromosomes from genotype data stored in a GDS file. The function
#' combines background LD estimation, chromosome-wise decay fitting, and
#' genome-wide summarization into a unified workflow.
#'
#' The procedure consists of:
#' \enumerate{
#'   \item Estimation of background LD from inter-chromosomal SNP pairs.
#'   \item Sliding-window estimation of LD-decay parameters per chromosome.
#'   \item Robust aggregation of decay parameters across windows.
#'   \item Prediction of decay parameters based on chromosome size.
#'   \item Derivation of recommended LD window sizes for user-defined LD thresholds (\eqn{\rho}).
#' }
#'
#' LD decay is modeled as:
#' \deqn{r^2(d) = b + \frac{c - b}{1 + a d}}
#'
#' where:
#' \itemize{
#'   \item \eqn{a} controls the rate of decay,
#'   \item \eqn{b} is background LD (long-distance baseline),
#'   \item \eqn{c} is short-range LD,
#'   \item \eqn{d} is physical distance (bp).
#' }
#'
#' @param gds A GDS file handle containing genotype data.
#' @param el_data_folder If not `NULL`, a directory to which the per-chromosome LD
#'   edge lists are written (one `<chr>.el` file each). These edge lists are the
#'   input to [compute_ld_w()] for the local-LD (`ld_w`) statistic; see
#'   `max_SNPs_decay` for the subsampling caveat that governs how complete they are.
#' @param q Quantile of \eqn{r^2} used for decay fitting (default = 0.95).
#' @param seed Optional integer. Fixes the random draws this function makes --
#'   the background-LD subsample, the `max_SNPs_decay` thinning and the
#'   stratified pair sampling -- so that repeat calls on the same data return
#'   the same fit. The caller's RNG stream is saved and restored, so a seed
#'   changes this function's behaviour and nothing else in a pipeline.
#'
#'   Worth setting whenever two runs are to be COMPARED rather than merely
#'   produced, because \eqn{b} propagates: it sets every decay-relative
#'   threshold through [ld_from_rho()], including the stage-1 join in
#'   [ld_complexity_reduction()] and the `min_r2` derived in
#'   [ld_prune_and_eMLG()]. Unseeded repeat parses of one simulated dataset
#'   moved \eqn{b} from 0.03038 to 0.03045, which was enough to move stage-1
#'   cluster counts by 0.75% (12,435 against 12,529) and GRM marker counts by
#'   0.70% -- larger than the differences some settings changes produce, so an
#'   unseeded A/B comparison in that range measures the subsample, not the
#'   setting.
#' @param n_sub_bg Number of SNPs used to estimate background LD.
#' @param n_win_decay Number of sliding windows per chromosome.
#' @param overlap Proportion of overlap between consecutive windows (0–1).
#' @param max_SNPs_decay Maximum number of SNPs sampled per chromosome when
#'   estimating the LD-decay curve (default `Inf`). \strong{Caution:} this cap does
#'   \emph{not} apply only to the decay fit -- the per-chromosome edge list is built
#'   from the same subsample, so when the edge lists are retained (`keep_el = TRUE`)
#'   or written (`el_data_folder`) the cap \emph{carries over into them}. Those edge
#'   lists are the input to [compute_ld_w()] downstream, so a finite
#'   `max_SNPs_decay` thins the `ld_w` support (markers outside the subsample get no
#'   neighbours) and reduces C-score sensitivity. Cap it (for speed) \emph{only}
#'   when the run is used solely to estimate LD decay; leave it at `Inf` whenever the
#'   edge lists feed `ld_w` estimation downstream.
#' @param prob_robust Central proportion of windows retained for robust summarization.
#' @param max_pairs Maximum number of SNP pairs per window used in decay fitting.
#' @param n_strata Number of distance strata used when subsampling SNP pairs.
#' @param ld_method LD statistic passed to \code{SNPRelate::snpgdsLDMat} (its
#'   result is squared to r^2 for decay/background fitting). Default \code{"corr"}
#'   -- the non-EM Pearson correlation of genotype dosages (equivalent to
#'   \code{cor()^2}, differing only in missing-data handling); \code{"r"} is the
#'   EM-based estimate.
#' @param keep_el Logical; whether to retain the per-chromosome LD edge lists in the
#'   returned object (consumed by [compute_ld_w()]). Subject to the same
#'   `max_SNPs_decay` subsampling caveat.
#' @param slide Sliding window size in number of SNPs used for LD estimation.
#' @param rho_targets Numeric vector of target LD thresholds used to derive recommended window sizes.
#' @param cores Number of CPU cores for parallel computation.
#' @param min_maf_decay Minor-allele-frequency threshold for decay fitting
#'   (default 0.1). When non-\code{NULL}, minor allele frequencies are computed
#'   directly from \code{gds} (via \code{SNPRelate::snpgdsSNPRateFreq()}) and only
#'   SNP pairs where BOTH members have \code{MAF > min_maf_decay} are used for
#'   background-LD and decay-curve fitting (low-MAF pairs mechanically deflate
#'   \eqn{r^2} regardless of true linkage, which biases the fitted a/b/c
#'   parameters). Filtering is applied only to the curve-fitting step; the
#'   returned edge lists (\code{by_chr[[ch]]$el}) still contain ALL SNPs,
#'   unfiltered, for downstream use (e.g. \code{compute_ld_w()}). Set to
#'   \code{NULL} to disable MAF filtering.
#' @param ld_w_rho Optional numeric vector of relative-LD levels \eqn{\rho}. When
#'   supplied, the local-LD statistic \code{ld_w} is computed \emph{in place} from
#'   the per-chromosome edge lists already built here (via [compute_ld_w()]) and
#'   returned as \code{$ld_ws}; the edge lists are then dropped (unless
#'   \code{keep_el = TRUE}). This is the cheapest route to \code{ld_w} -- the edge
#'   lists are reused, not recomputed or saved -- so with \code{keep_el = FALSE} and
#'   no \code{el_data_folder} nothing large is written or retained. \code{NULL}
#'   (default) skips it.
#'
#' @return An object of class \code{"ld_decay"} containing (and, when
#'   \code{ld_w_rho} is set, an additional \code{ld_ws} matrix of local-LD support,
#'   markers x \eqn{\rho}):
#' \describe{
#'   \item{by_chr}{List of per-chromosome results including decay fits and optional LD edge lists.}
#'   \item{decay_sum}{Data table of chromosome-wise decay parameters and derived quantities.}
#'   \item{decay_model}{Robust regression model linking decay rate \eqn{a} to chromosome size.}
#'   \item{recommendation}{Suggested LD window sizes for specified \eqn{\rho} thresholds.}
#'   \item{params}{List of parameters used in the computation.}
#' }
#'
#' @details
#' The function estimates LD decay independently per chromosome and then
#' fits a genome-wide model relating decay rate (\eqn{a}) to chromosome size.
#' This allows extrapolation of LD behavior and derivation of consistent
#' window sizes across heterogeneous genomic architectures.
#'
#' Recommended window sizes are provided in SNP units to match downstream
#' functions that operate on marker indices rather than physical distance.
#'
#' @section Comparing decay estimates across runs and datasets:
#' The fitted parameters (\eqn{a}, \eqn{b}, \eqn{c}, and the \code{a_pred}
#' derived from them) are \strong{conditional on the settings used}, not
#' properties of the data alone. They are stable when those settings are held
#' fixed: on one simulated dataset, repeat runs differing only in random seed
#' gave \eqn{a} = 1.798e-05 and 1.799e-05. Note that this stability is a
#' property of \eqn{a}, not of the whole fit -- \eqn{b} is estimated from a
#' subsample and moves enough between unseeded runs to shift the thresholds
#' derived from it (see `seed`), so use `seed` when comparing runs rather than
#' relying on the fit being reproducible by default. But they move substantially when the
#' settings -- or the marker set -- change, and \eqn{a} propagates: it sets the
#' \code{d_from_rho()} window behind \code{ld_w}, and \code{decay_sum} feeds
#' \code{score_thresholds()}. Treat estimates as comparable only across runs
#' that share the whole configuration.
#'
#' Matching the settings is necessary but \strong{not sufficient}, because
#' \code{slide} is counted in SNPs: the physical span it covers is
#' \code{slide * bp_per_snp}, so the same nominal window means different things
#' at different marker densities. Halving the markers of one dataset while
#' changing nothing else (via \code{max_SNPs_decay}) moved \eqn{a} from
#' 1.80e-05 to 2.2-2.6e-05 and pushed \eqn{c} from ~0.97 to ~0.999 -- density
#' alone. To compare across datasets, check that \code{decay_sum$slide_bp} and
#' \code{rho_slide_pred} are similar, rather than only that \code{slide}
#' matches; \code{recommendation$suggested_slide_by_chr} gives the
#' density-aware slide for a target \eqn{\rho}.
#'
#' A related consequence: the estimator is sensitive to how pairs are
#' \emph{composed}, not just how many there are, since
#' \code{subsample_pairs_for_decay()} derives its distance strata from the
#' observed range of the pairs it is given. Feeding it a differently sampled
#' pair set (a marker subsample, or pairs drawn directly) shifts the strata and
#' therefore the fit, even at equal pair counts.
#'
#' @export
compute_LD_decay <- function(
    gds,
    el_data_folder = NULL,
    q = 0.95,
    n_sub_bg = 5000,
    n_win_decay = 20,
    overlap = 0.5,
    max_SNPs_decay = Inf,
    prob_robust = 0.95,
    max_pairs = 5000,
    ld_method = "corr",
    n_strata = 20,
    keep_el = FALSE,
    slide = 1000,
    rho_targets = c(0.90, 0.95, 0.99),
    cores = 1,
    min_maf_decay = 0.1,
    ld_w_rho = NULL,
    seed = NULL
) {

  ## Everything random in this function happens downstream of here: the
  ## background-LD subsample (n_sub_bg), the per-chromosome marker thinning
  ## (max_SNPs_decay) and the stratified pair sampling (max_pairs). Seeding once
  ## at the top therefore fixes the whole fit.
  ##
  ## The caller's RNG stream is restored on exit, so passing a seed changes what
  ## THIS function does and nothing else. Without that, seeding here would
  ## silently re-align every later draw in a pipeline that also samples.
  if (!is.null(seed)) {
    if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) stats::runif(1)
    .old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", .old_seed, envir = globalenv()), add = TRUE)
    set.seed(seed)
  }

  if (!is.null(el_data_folder)) {
    if (!dir.exists(el_data_folder)) dir.create(el_data_folder, recursive = TRUE)
    keep_el = FALSE # don't keep if saved
    message("Savinge edge list to folder ",el_data_folder)
  }

  ids  <- .read_gds_ids(gds)
  chrs <- unique(ids$snp_chr)

  ## ---- MAF filtering for decay-curve fitting.
  ## Active whenever `min_maf_decay` is non-NULL. Minor allele frequencies are
  ## computed directly from the GDS (via snpgdsSNPRateFreq), so callers do not
  ## supply them. Filtering affects only background-LD and decay-curve fitting;
  ## edge lists (by_chr[[ch]]$el) remain unfiltered for downstream use.
  maf <- NULL
  if (!is.null(min_maf_decay)) {
    maf <- SNPRelate::snpgdsSNPRateFreq(gds)$MinorFreq
    names(maf) <- ids$snp_id
    message("Filtering SNPs to MAF > ", min_maf_decay,
            " for background-LD and decay-curve fitting (edge lists remain unfiltered).")
  }

  b <- estimate_background_ld(gds, n_sub = n_sub_bg, q = q, ld_method = ld_method,maf = maf, min_maf = min_maf_decay)

  out_by_chr <- vector("list", length(chrs))
  names(out_by_chr) <- chrs


  message("Estimating LD-decay")
  pb <- txtProgressBar(min = 0, max = length(chrs)-1, style = 3)
  setTxtProgressBar(pb, 0)
  #ch = "Chr1"
  for (ch in chrs) {

    chr_idx  <- which(ids$snp_chr == ch)
    pos_chr  <- ids$snp_pos[chr_idx]
    snps_chr <- ids$snp_id[chr_idx]

    ## robust chromosome-specific bp per SNP scale
    chr_size_bp <- max(pos_chr) - min(pos_chr)
    n_snps_chr  <- length(chr_idx)

    bp_per_snp <- if (n_snps_chr > 1) {
      chr_size_bp / (n_snps_chr - 1)
    } else {
      NA_real_
    }

    ## current slide converted to bp
    slide_bp <- slide * bp_per_snp

    sample_idx <- sort(sample(
      chr_idx,
      size = min(max_SNPs_decay, length(chr_idx)),
      replace = FALSE
    ))

    el <- get_el(
      gds = gds,
      idx = sample_idx,
      slide_win_ld = slide,   # still in SNPs for get_el
      cores = cores,
      by_chr = TRUE,
      method = ld_method
    )

    data.table::setkey(el, SNP1)

    window_size <- chr_size_bp / n_win_decay
    step_size   <- window_size * overlap

    ## ---- MAF-filtered view of the edge list, used ONLY for decay-curve
    ## fitting below. `el` itself is left untouched (all SNPs, all MAFs)
    ## so downstream ld_w computation can apply its own MAF filter later.
    ##
    ## Both filters (MAF and the window-size cut) are applied in ONE pass, and
    ## only the four columns the fit actually reads are carried over. Previously
    ## this made a full-width copy of the edge list and then a second copy for
    ## `d < window_size`, dragging SNP1/SNP2/Chr1/Chr2 -- four character columns,
    ## roughly two thirds of the table's memory -- through both, despite
    ## estimate_decay_chr() and coef_ld_dec() only ever using pos1/pos2/d/r2.
    if (!is.null(maf)) {
      hi_maf_ids <- snps_chr[maf[snps_chr] > min_maf_decay & !is.na(maf[snps_chr])]
      el_decay <- el[SNP1 %in% hi_maf_ids & SNP2 %in% hi_maf_ids & d < window_size,
                     list(pos1, pos2, r2, d)]
    } else {
      el_decay <- el[d < window_size, list(pos1, pos2, r2, d)]
    }

    ## ---- LD-decay
    decay <- suppressWarnings(
      estimate_decay_chr(
        el = el_decay,
        b = b,
        window_size = window_size,
        step_size   = step_size,
        max_pairs   = max_pairs,
        n_strata    = n_strata,
        q = q,
        cores = cores
      )
    )


    if(!is.null(el_data_folder)){
      ## save as text file to el_data_folder
      fwrite(el, paste0(el_data_folder,ch,".el"))
    }

    decay[, contrast := c - b]

    decay[, regime := ifelse(
      contrast < 0.05,
      "weak",
      "structured"
    )]

    decay_sum_chr <- summarize_decay(decay[regime == "structured", ])

    if (is.null(decay_sum_chr)) next

    decay_sum_chr[, chr_size := max(pos_chr)]
    decay_sum_chr[, n_snp_chr := n_snps_chr]
    decay_sum_chr[, bp_per_snp := bp_per_snp]
    decay_sum_chr[, slide_snp := slide]
    decay_sum_chr[, slide_bp := slide_bp]
    decay_sum_chr[, n_w_used := sum(na.omit(decay$regime == "structured"))]
    decay_sum_chr[, Chr := ch]

    ## rho actually covered by the user-supplied slide, using raw chromosome a
    if (is.finite(decay_sum_chr$a) && decay_sum_chr$a > 0 &&
        is.finite(slide_bp) && slide_bp > 0) {
      decay_sum_chr[, rho_slide_raw := (a * slide_bp) / (1 + a * slide_bp)]
    } else {
      decay_sum_chr[, rho_slide_raw := NA_real_]
    }
    setTxtProgressBar(pb, which(chrs==ch))
    out_by_chr[[ch]] <- list(
      snp_ids   = snps_chr,
      ## retain the edge list in memory only when the caller asked for it
      ## (keep_el); otherwise store the on-disk path -- but only when there IS one:
      ## with no el_data_folder, paste0() would produce a bare "<chr>.el" that was
      ## never written, which downstream code then tries (and fails) to fread.
      ## NULL correctly says "no edges here".
      ##
      ## The in-place ld_w pass below deliberately does NOT hold these: it cannot
      ## run inside this loop, because ld_w is defined against `a_pred` -- the
      ## genome-wide regression of decay rate on chromosome size, which only
      ## exists once every chromosome has been fitted. Retaining every
      ## chromosome's edges until then made peak memory scale with the number of
      ## chromosomes; instead the ld_w pass rebuilds them one chromosome at a
      ## time from `gds` and discards each after use, so the peak is one edge
      ## list regardless of how many chromosomes there are.
      el        = if (keep_el) el
                  else if (!is.null(el_data_folder)) paste0(el_data_folder, ch, ".el")
                  else NULL,
      decay     = decay,
      decay_sum = decay_sum_chr
    )
  }
  close(pb)

  ## A chromosome whose decay fit failed summarize_decay()'s gate (fewer than 5
  ## windows with a usable fit) was skipped by `next` above -- but out_by_chr is
  ## PRE-ALLOCATED with names, so that leaves a NULL element still present in the
  ## list. Consumers then dereference it: compute_ld_w() evaluating
  ## `decay_sum[Chr == chr_obj$decay_sum$Chr]` on a NULL entry yields a length-0
  ## RHS and data.table's opaque "RHS of == is length 0" error, far from the
  ## actual cause. Drop the holes here and say plainly which chromosomes are gone.
  failed_chr <- names(out_by_chr)[vapply(out_by_chr, is.null, logical(1))]
  if (length(failed_chr)) {
    out_by_chr <- out_by_chr[!vapply(out_by_chr, is.null, logical(1))]
    if (!length(out_by_chr)) {
      stop("LD-decay estimation failed for every chromosome (", paste(failed_chr, collapse = ", "),
           "): no chromosome had >= 5 windows with a valid fit. Nothing can be estimated from ",
           "this dataset at these settings -- see summarize_decay() and n_win_decay.", call. = FALSE)
    }
    warning("LD-decay estimation failed for chromosome(s) ", paste(failed_chr, collapse = ", "),
            " (fewer than 5 windows with a valid fit). They are ABSENT from by_chr, decay_sum ",
            "and ld_ws, so downstream results cover only: ",
            paste(names(out_by_chr), collapse = ", "), ".", call. = FALSE)
  }

  message("Predicting a from chromosome size")
  decay_sum <- data.table::rbindlist(
    lapply(out_by_chr, function(x) x$decay_sum),
    fill = TRUE
  )


    ## robust background decay model
    if(length(unique(decay_sum$chr_size))<3){
      message("Not enough unique chromosome sizes to predict, using median")
      d_mod <- NA
      decay_sum[, a_pred := median(a, na.rm = TRUE)]
      decay_sum[, c_pred := median(c, na.rm = TRUE)]
    }else{
      d_mod <- MASS::rlm(log(a) ~ log(chr_size), data = decay_sum)

      decay_sum[, a_pred := exp(predict(d_mod, decay_sum))]
      decay_sum[, c_pred := median(c, na.rm = TRUE)]

    }


  ## rho covered by the current slide using predicted background a
  decay_sum[, rho_slide_pred := (a_pred * slide_bp) / (1 + a_pred * slide_bp)]

  ## recommended slide for each requested rho target
  rho_targets <- sort(unique(rho_targets))
  rho_targets <- rho_targets[is.finite(rho_targets) & rho_targets > 0 & rho_targets < 1]

  for (rho_t in rho_targets) {
    nm_bp  <- paste0("slide_bp_rho_",  formatC(rho_t, format = "f", digits = 2))
    nm_snp <- paste0("slide_snp_rho_", formatC(rho_t, format = "f", digits = 2))

    ## required physical window
    decay_sum[, (nm_bp) := rho_t / (a_pred * (1 - rho_t))]

    ## convert required bp window back to SNP window
    decay_sum[, (nm_snp) := get(nm_bp) / bp_per_snp]
  }

  ## optional compact summary for user-facing reporting
  rec_cols_snp <- grep("^slide_snp_rho_", names(decay_sum), value = TRUE)

  recommendation <- list(
    slide_input_snp = slide,
    rho_targets = rho_targets,

    rho_covered = decay_sum[, .(
      Chr,
      chr_size,
      n_snp_chr,
      bp_per_snp,
      a_pred,
      slide_snp,
      slide_bp,
      rho_slide_pred
    )],

    suggested_slide_by_chr = decay_sum[, c(
      "Chr", "chr_size", "n_snp_chr", "bp_per_snp", "a_pred", rec_cols_snp
    ), with = FALSE],

    suggested_slide_summary = if (length(rec_cols_snp) > 0) {
      tmp <- lapply(rec_cols_snp, function(cc) {
        vals <- decay_sum[[cc]]
        data.table::data.table(
          target = sub("^slide_snp_rho_", "", cc),
          median = stats::median(vals, na.rm = TRUE),
          p90    = stats::quantile(vals, probs = 0.90, na.rm = TRUE),
          max    = max(vals, na.rm = TRUE)
        )
      })
      data.table::rbindlist(tmp)
    } else {
      NULL
    }
  )

  out <- list(
    by_chr         = out_by_chr,
    decay_sum      = decay_sum,
    decay_model    = d_mod,
    recommendation = recommendation,
    params         = list(
      q = q,
      n_sub_bg = n_sub_bg,
      n_win_decay = n_win_decay,
      overlap = overlap,
      el_data_folder = el_data_folder,
      n_strata = n_strata,
      max_pairs = max_pairs,
      keep_el = keep_el,
      prob_robust = prob_robust,
      slide = slide,             # in SNPs
      ld_method = ld_method,     # so edges rebuilt later match the ones fitted here
      min_maf_decay = min_maf_decay,
      max_SNPs_decay = max_SNPs_decay,
      rho_targets = rho_targets,
      cores = cores
    )
  )

  class(out) <- "ld_decay"

  ## optional in-place ld_w. Edge lists are reused when they are actually held
  ## (keep_el, or an el_data_folder path); otherwise each chromosome's is rebuilt
  ## from `gds` inside compute_ld_w() and dropped again, keeping peak memory at
  ## one chromosome's edges rather than all of them.
  if (!is.null(ld_w_rho)) {
    reuse <- keep_el || !is.null(el_data_folder)
    out$ld_ws <- compute_ld_w(out, rho = ld_w_rho, cores = cores,
                              gds = if (reuse) NULL else gds,
                              slide_win_ld = slide, ld_method = ld_method)
    if (!isTRUE(keep_el)) for (ch in names(out$by_chr)) out$by_chr[[ch]]$el <- NULL
  }

  if (!is.null(recommendation$suggested_slide_summary)) {
    message("Current slide covers the following rho range (using predicted background a):")
    print(decay_sum[, .(
      min_rho = min(rho_slide_pred, na.rm = TRUE),
      median_rho = median(rho_slide_pred, na.rm = TRUE),
      max_rho = max(rho_slide_pred, na.rm = TRUE)
    )])

    message("Suggested slide windows in SNPs for target rho:")
    print(recommendation$suggested_slide_summary)
  }

  out
}


#' Summarize LD Decay Parameters Across Windows
#'
#' Computes robust summary statistics of LD-decay parameters estimated
#' across sliding windows within a chromosome.
#'
#' Extreme windows are removed using symmetric quantile trimming of the
#' decay parameter \eqn{a}, and median values are returned.
#'
#' @param decay_dt Data table containing window-wise decay parameters
#'   (columns \code{a}, \code{b}, \code{c}).
#' @param prob_robust Central proportion of windows retained (default = 0.95).
#'
#' @return A \code{data.table} with median estimates of:
#' \describe{
#'   \item{a}{Decay rate.}
#'   \item{b}{Background LD.}
#'   \item{c}{Short-range LD.}
#' }
#'
#' @details
#' Trimming is based on the distribution of \eqn{a}, which is typically
#' the most variable parameter across genomic windows.
#'
#' Returns \code{NULL} if insufficient valid windows are available.
summarize_decay <- function(decay_dt, prob_robust = 0.95) {

  decay_valid <- decay_dt[!is.na(a) & a > 0 & !is.na(c)]

  if (nrow(decay_valid) < 5)
    return(NULL)

  # compute symmetric trimming bounds
  alpha <- (1 - prob_robust) / 2

  q_lo <- quantile(decay_valid$a, alpha, na.rm = TRUE)
  q_hi <- quantile(decay_valid$a, 1 - alpha, na.rm = TRUE)

  decay_trim <- decay_valid[a >= q_lo & a <= q_hi]

  if (nrow(decay_trim) < 3)
    return(NULL)

  data.table(
    c = median(decay_trim$c, na.rm = TRUE),
    a = median(decay_trim$a, na.rm = TRUE),
    b = decay_trim$b[1]
  )
}

subsample_pairs_for_decay <- function(sub,
                                      max_pairs = 5000,
                                      n_strata = 20) {


  if (nrow(sub) <= max_pairs) max_pairs <- nrow(sub)
  # log-distance strata
  log_d <- log(sub$d)
  breaks <- seq(min(log_d), max(log_d), length.out = n_strata + 1)

  sub[, strata := cut(log_d, breaks = breaks, include.lowest = TRUE)]

  target_per_stratum <- floor(max_pairs / n_strata)

  sub_sampled <- sub[, {
    if (.N > target_per_stratum) {
      .SD[sample(.N, target_per_stratum)]
    } else {
      .SD
    }
  }, by = strata]


  return(sub_sampled)
}

#' Estimate LD Decay Within a Chromosome
#'
#' Fits LD-decay models in sliding genomic windows for a single chromosome
#' using pairwise LD data (edge list format).
#'
#' For each window:
#' \enumerate{
#'   \item SNP pairs are restricted to the window.
#'   \item Pairs are subsampled across distance strata.
#'   \item LD decay is fitted using nonlinear regression.
#' }
#'
#' @param el Data table of SNP pairs with LD values and positions.
#' @param b Background LD.
#' @param window_size Window size in base pairs.
#' @param step_size Step size between windows in base pairs.
#' @param q Quantile of \eqn{r^2} used for fitting.
#' @param max_pairs Maximum number of SNP pairs per window.
#' @param n_strata Number of distance strata for subsampling.
#' @param cores Number of CPU cores.
#'
#' @return A \code{data.table} with window coordinates and decay parameters.
#'
#' @details
#' Windows with insufficient data or failed fits return partial rows with
#' missing parameter estimates.
estimate_decay_chr <- function(el,
                               b,
                               window_size,
                               step_size,
                               q = 0.95,
                               max_pairs=5000,
                               n_strata = 20,
                               cores = 1) {



  min_pos <- min(c(el$pos1, el$pos2))
  max_pos <- max(c(el$pos1, el$pos2))

  starts <- seq(min_pos, max_pos - window_size, by = step_size)
  ends   <- starts + window_size
  #i <- 1
  decay <- suppressWarnings(rbindlist(lapply(seq_along(starts), function(i) {

    sub <- el[
      pos1 >= starts[i] & pos1 < ends[i] &
        pos2 >= starts[i] & pos2 < ends[i]
    ]

    if(!nrow(sub)>100) return(NULL)

    sub <- subsample_pairs_for_decay(
      sub,max_pairs = max_pairs, n_strata = n_strata
    )

    coefs <- tryCatch(
      coef_ld_dec(dt_strata=sub, q = q,  b = b),
      error = function(e) NULL
    )

    if (is.null(coefs)) {
      data.table(start = starts[i], end = ends[i])
    } else {
      data.table(start = starts[i], end = ends[i],coefs)
    }
  }),fill=TRUE,use.names = TRUE))

  return(decay)
}

#' Estimate Background LD from Inter-chromosomal SNP Pairs
#'
#' Computes background LD as a high quantile of \eqn{r^2} between SNPs
#' located on different chromosomes.
#'
#' @param gds GDS file handle.
#' @param idx Optional vector of SNP indices.
#' @param n_sub Number of SNPs sampled for estimation.
#' @param q Quantile used to define background LD.
#' @param maf Optional named numeric vector of MAF, named by SNP id. If supplied
#'   together with \code{min_maf}, only SNPs with \code{MAF > min_maf} are
#'   eligible for sampling (low-MAF pairs mechanically deflate r^2 and bias the
#'   background estimate).
#' @param min_maf MAF threshold used with \code{maf}.
#' @param ld_method LD statistic passed to \code{SNPRelate::snpgdsLDMat} (squared
#'   to r^2). Default \code{"corr"} (non-EM Pearson correlation of dosages).
#'
#' @return Numeric scalar representing background LD (\eqn{b}).
#'
#' @details
#' Inter-chromosomal LD provides an empirical estimate of baseline
#' correlation unrelated to physical linkage. This value anchors the
#' LD-decay model and stabilizes parameter estimation.
#'
#' Sampling is performed proportionally across chromosomes.
estimate_background_ld <- function(gds,
                                   idx=NULL,
                                   n_sub = 5000,
                                   q = 0.95,
                                   ld_method="corr",
                                   maf = NULL,
                                   min_maf = NULL) {


  ids <- .read_gds_ids(gds)

  if (is.null(idx)) idx <- seq_along(ids$snp_id)

  if (!is.null(maf) && !is.null(min_maf)) {
    maf_ok <- maf[ids$snp_id[idx]] > min_maf
    maf_ok[is.na(maf_ok)] <- FALSE
    idx <- idx[maf_ok]
  }

  n_snps <- length(idx)
  if (n_snps < 2L) {
    stop("Not enough SNPs to estimate background LD (after MAF filtering, if applied).")
  }

  chr_vec <- ids$snp_chr[idx]
  chr_levels <- unique(chr_vec)

  if (length(chr_levels) < 2L) {
    stop("Need at least two chromosomes to estimate background LD.")
  }


  ## proportional sampling across chromosomes
  snp_pool <- unlist(lapply(chr_levels, function(ch) {

    ix <- which(chr_vec == ch)
    n_ch <- length(ix)

    sample(ix,size = min(n_ch, ceiling(n_sub * n_ch / n_snps)),replace = FALSE)

  }))

  if (length(unique(chr_vec[snp_pool])) < 2L) {
    stop("Sampled SNPs fall on a single chromosome; increase n_sub.")
  }

  ld <- SNPRelate::snpgdsLDMat(
    gds,
    snp.id = ids$snp_id[idx][snp_pool],
    method = ld_method,
    slide = -1,
    verbose = FALSE
  )

  r2 <- ld$LD^2
  chr_sub <- chr_vec[snp_pool]

  ## upper triangle only
  inter_idx <- which(outer(chr_sub, chr_sub, FUN = "!=") & upper.tri(r2))
  r2_inter <- r2[inter_idx]

  if (length(r2_inter) == 0L) {
    stop("No inter-chromosomal SNP pairs found.")
  }

  b <- quantile(r2_inter, probs = q, na.rm = TRUE)

  message("Background LD: ", round(b, 3))
  return(b)
}



#' Fit LD Decay Model to Stratified LD Data
#'
#' Fits a nonlinear LD-decay model to binned LD values derived from
#' distance-stratified SNP pairs.
#'
#' The model is:
#' \deqn{r^2(d) = b + \frac{c - b}{1 + a d}}
#'
#' @param dt_strata Data table of SNP pairs with distance strata.
#' @param q Quantile used to summarize LD within strata.
#' @param b Background LD.
#'
#' @return A \code{data.table} containing:
#' \describe{
#'   \item{a}{Decay rate.}
#'   \item{c}{Short-range LD.}
#'   \item{b}{Background LD.}
#'   \item{agg}{Binned data used for fitting.}
#'   \item{raw}{Original stratified data.}
#' }
#'
#' @details
#' The model is fitted using weighted nonlinear least squares,
#' where weights correspond to the number of SNP pairs per stratum.
#'
#' Returns \code{NULL} if fitting fails or insufficient data are available.
coef_ld_dec <- function(dt_strata,
                        q = 0.95,
                        b) {

  if (nrow(dt_strata) < 100) return(NULL)
  #dt_strata <- sub
  agg <- dt_strata[, .(
    d_mid = exp(mean(log(d))),   # geometric midpoint
    r2_q  = quantile(r2, q, na.rm = TRUE),
    n     = .N
  ), by = strata]


  if (nrow(agg) < 10) return(NULL)

  # Initial parameter guesses
  c_start <- max(agg$r2_q, na.rm = TRUE)
  a_start <- 1 / median(agg$d_mid)

  ld_from_rho
  fit <- tryCatch(
    nls(
      r2_q ~ b + (c - b)/(1 + a * d_mid),
      data = agg,
      start = list(c = min(max(agg$r2_q), 1), a = a_start),
      algorithm = "port",
      lower = c(c = b, a = 0),
      upper = c(c = 1, a = Inf),
      weights = sqrt(n),
      control = nls.control(warnOnly = TRUE)
    ),
    error = function(e) NULL
  )


  if (is.null(fit)) return(NULL)

  coefs <- coef(fit)

  return(data.table(c=coefs["c"],a=coefs["a"],b,agg=list(agg),raw=list(dt_strata)))

}



#' Convert Relative LD Threshold to Physical Distance
#'
#' Computes the physical distance corresponding to a relative LD threshold \eqn{\rho}.
#'
#' @param a Decay rate.
#' @param rho Relative LD threshold (0 < rho < 1).
#'
#' @return Distance in base pairs.
#'
#' @export
d_from_rho <- function(a, rho){
  rho / (a * (1 - rho))
}

#' Convert Relative LD Threshold to Expected LD Value
#'
#' Computes the expected \eqn{r^2} corresponding to a relative LD threshold \eqn{\rho}.
#'
#' @param b Background LD.
#' @param c Short-range LD (default = 1).
#' @param rho Relative LD threshold.
#'
#' @return Expected \eqn{r^2}.
#'
#' @export
ld_from_rho <- function(b, c = 1, rho){
  b + (c - b) * (1 - rho)
}

#' Interpolate Genetic Map Position (cM)
#'
#' Linearly interpolates centimorgan (cM) position from a genetic map, for
#' arbitrary physical positions. Converts physical (bp) distance between
#' markers into genetic (cM) distance -- a much more direct measure of
#' recombination probability than physical distance, which can be
#' misleading wherever recombination rate itself varies (e.g. across a
#' low-recombination block like a centromere or inversion versus the
#' surrounding chromosome). Checked directly on real data: some marker
#' pairs over 1 Mb apart physically sat at cM distance 0 (fully linked, no
#' measurable recombination between them), while others only a few hundred
#' kb apart already spanned several cM.
#'
#' @param genetic_map A data.table/data.frame with `Chr`, `Pos`, and `cM`
#'   columns (at least 2 rows per chromosome to interpolate). Any
#'   chromosome-naming differences from `Chr` below (e.g. a genetic map
#'   using `"chromosome_1"` instead of `"Chr1"`) must be reconciled by the
#'   caller before calling this function.
#' @param Chr Character vector of chromosome names to query.
#' @param Pos Numeric vector of physical positions to query, same length
#'   and order as `Chr`.
#'
#' @return Numeric vector of interpolated cM positions, same length/order
#'   as `Chr`/`Pos`. Positions outside a chromosome's mapped range are
#'   clamped to the nearest endpoint's cM value (`rule = 2` in
#'   [stats::approx()]) -- fine for markers just past the map's first/last
#'   point, but a chromosome entirely absent from `genetic_map`, or with
#'   fewer than 2 points, returns `NA` for all its positions.
#'
#' @examples
#' \dontrun{
#' genetic_map <- data.table::data.table(
#'   Chr = "Chr1", Pos = c(0, 1e6, 2e6), cM = c(0, 5, 20)
#' )
#' interpolate_cM(genetic_map, Chr = "Chr1", Pos = c(5e5, 1.5e6))
#' }
#'
#' @export
interpolate_cM <- function(genetic_map, Chr, Pos) {
  genetic_map <- data.table::as.data.table(genetic_map)
  stopifnot(all(c("Chr", "Pos", "cM") %in% names(genetic_map)))

  by_chr <- split(genetic_map, by = "Chr")
  out <- rep(NA_real_, length(Pos))

  for (ch in unique(Chr)) {
    idx <- which(Chr == ch)
    m <- by_chr[[ch]]
    if (is.null(m) || nrow(m) < 2L) next
    out[idx] <- stats::approx(m$Pos, m$cM, xout = Pos[idx], rule = 2)$y
  }

  out
}

parallel_apply <- function(X, FUN, cores = 1) {

  if (cores > 1 && .Platform$OS.type != "windows") {
    parallel::mclapply(X, FUN, mc.cores = cores)
  } else {
    lapply(X, FUN)
  }
}

#' Print an `ld_decay` object
#'
#' @param x An object of class `ld_decay`.
#' @param digits Number of significant digits to display.
#' @param ... Ignored (for S3 `print` consistency).
#' @return Invisibly, `x`.
#' @method print ld_decay
#' @export
print.ld_decay <- function(x, digits = 3, ...) {

  cat("<ld_decay>\n")

  ## ---- parameters
  if (!is.null(x$params)) {

    cat("\nRun parameters:\n")

    cat(
      "  slide window:",
      format(x$params$slide, big.mark = ","),
      "SNPs\n"
    )

    cat(
      "  background LD quantile q:",
      signif(x$params$q, digits),
      "\n"
    )

    cat(
      "  chromosomes analysed:",
      length(x$by_chr),
      "\n"
    )

    cat(
      "  keep_el:",
      x$params$keep_el,
      "\n"
    )
  }

  ds <- x$decay_sum

  if (!is.null(ds) && nrow(ds) > 0) {

    cat("\nDecay parameter summary:\n")

    if ("a_pred" %in% names(ds)) {

      cat(
        "  predicted a:",
        "median =", signif(stats::median(ds$a_pred, na.rm = TRUE), digits),
        " range = [",
        signif(min(ds$a_pred, na.rm = TRUE), digits), ", ",
        signif(max(ds$a_pred, na.rm = TRUE), digits), "]\n",
        sep = ""
      )
    }

    if ("c_pred" %in% names(ds)) {

      cat(
        "  predicted c:",
        signif(stats::median(ds$c_pred, na.rm = TRUE), digits),
        "\n"
      )
    }

    if ("n_w_used" %in% names(ds)) {

      cat(
        "  structured decay windows:",
        "median =", stats::median(ds$n_w_used, na.rm = TRUE),
        " range = [",
        min(ds$n_w_used, na.rm = TRUE), ", ",
        max(ds$n_w_used, na.rm = TRUE), "]\n",
        sep = ""
      )
    }
  }

  ## ---- rho coverage
  rec <- x$recommendation

  if (!is.null(rec)) {

    rho_cov <- rec$rho_covered

    if (!is.null(rho_cov) && "rho_slide_pred" %in% names(rho_cov)) {

      cat("\nCurrent slide window coverage:\n")

      cat(
        "  rho covered:",
        "median =", signif(stats::median(rho_cov$rho_slide_pred, na.rm = TRUE), digits),
        " range = [",
        signif(min(rho_cov$rho_slide_pred, na.rm = TRUE), digits), ", ",
        signif(max(rho_cov$rho_slide_pred, na.rm = TRUE), digits), "]\n",
        sep = ""
      )

      cat(
        "  (slide =", format(x$params$slide, big.mark = ","), "SNPs)\n"
      )
    }

    ## ---- suggested slide windows
    ss <- rec$suggested_slide_summary

    if (!is.null(ss) && nrow(ss) > 0) {

      cat("\nSuggested slide windows for target rho:\n")

      for (i in seq_len(nrow(ss))) {

        cat(
          "  rho =", ss$target[i],
          ": median =", format(round(ss$median[i]), big.mark = ","), "SNPs",
          ", p90 =", format(round(ss$p90[i]), big.mark = ","), "SNPs",
          ", max =", format(round(ss$max[i]), big.mark = ","), "SNPs",
          "\n"
        )
      }

      cat(
        "  (median = typical chromosome, p90 = covers most chromosomes)\n"
      )
    }
  }

  ## ---- stored components
  cat("\nStored components:\n")

  cat(
    "  by_chr, decay_sum, decay_model, recommendation, params\n"
  )

  invisible(x)
}


#' Compute Local LD Support (ld_w)
#'
#' Computes per-SNP local LD support within a distance defined by one or
#' more relative LD thresholds \eqn{\rho}.
#'
#' For each SNP, LD support is defined as the median \eqn{r^2} with
#' neighboring SNPs within the corresponding distance window.
#'
#' @param ld_decay Object of class \code{"ld_decay"}.
#' @param rho Relative LD threshold used to define the window. May be a
#'   vector of multiple thresholds -- each chromosome's edge list is read
#'   (or pulled from memory) and symmetrized once and reused for every
#'   \code{rho}, rather than repeating that work once per threshold. This
#'   matters most when \code{keep_el = FALSE} was used in
#'   \code{compute_LD_decay()}, so each chromosome's edge list has to be
#'   re-read from disk.
#' @param cores Number of CPU cores.
#' @param gds Optional GDS handle. If supplied, each chromosome's edge list is
#'   recomputed \emph{on the fly} from the genotypes (via [get_el()]) and discarded
#'   after use, instead of reading a stored edge list -- so \code{compute_LD_decay()}
#'   can be run with \code{keep_el = FALSE} and no \code{el_data_folder}, avoiding
#'   the (large) edge-list files entirely. Any stored edge list is ignored when
#'   \code{gds} is given.
#' @param slide_win_ld,ld_method Passed to [get_el()] for the on-the-fly mode: the
#'   sliding-window size (in SNPs) for \code{SNPRelate::snpgdsLDMat} and the LD
#'   statistic (default \code{"corr"}). Use the same values as the
#'   \code{compute_LD_decay()} run that produced \code{ld_decay}.
#'
#' @return If \code{rho} has length 1, a named numeric vector of LD support
#'   values (names = SNP/marker ids, in the same SNP order as
#'   \code{ld_decay$by_chr}). If \code{rho} has length > 1, a numeric matrix with
#'   one row per SNP (rownames = marker ids, same order) and one column per
#'   \code{rho}, named \code{"rho_<value>"}. The row names let callers index/align
#'   by marker directly (e.g. \code{ld_ws[map$marker, ]}) without setting them.
#'
#' @details
#' The physical window is derived using:
#' \deqn{d = \frac{\rho}{a(1 - \rho)}}
#'
#' where \eqn{a} is the chromosome-specific decay rate.
#'
#' The edge lists come from \code{compute_LD_decay()} (with \code{keep_el = TRUE} or
#' an \code{el_data_folder}); alternatively, supply \code{gds} to recompute them on
#' the fly and avoid saving them at all (see the \code{gds} argument).
#'
#' @export
compute_ld_w <- function(
    ld_decay,
    rho = 0.95,
    cores = 1,
    gds = NULL,
    slide_win_ld = 1000,
    ld_method = "corr"
) {
  rho <- unique(rho)

  #chr_obj <- ld_decay$by_chr$Chr1
  by_chr_w <- parallel_apply(ld_decay$by_chr, function(chr_obj) {

    a <- ld_decay$decay_sum[Chr==chr_obj$decay_sum$Chr,a_pred]

    ## Edge list: either the one stored by compute_LD_decay (in memory or on disk),
    ## or -- when `gds` is supplied -- recomputed on the fly for this chromosome's
    ## markers and discarded after use, so the edge lists never have to be saved.
    el <- if (!is.null(gds)) {
      get_el(gds, SNP_id = chr_obj$snp_ids, slide_win_ld = slide_win_ld,
             method = ld_method, cores = 1)
    } else if (!is.null(chr_obj$el)) {
      if (is.character(chr_obj$el)) fread(chr_obj$el, showProgress = FALSE) else chr_obj$el
    } else {
      stop("No edge list present; supply `gds` to estimate ld_w on the fly, or ",
           "run compute_LD_decay() with keep_el = TRUE or an el_data_folder.")
    }

    ## Make symmetric -- once per chromosome, reused across all rho below. This
    ## table has TWICE the edge list's rows, so it carries only the three columns
    ## the median below actually reads; pos/pos_other used to be included and were
    ## never used, costing ~40% of the largest allocation in this function.
    el <- data.table::rbindlist(list(
      el[, .(SNP = SNP1, r2, d)],
      el[, .(SNP = SNP2, r2, d)]
    ))

    w <- vapply(rho, function(r) {
      d_window <- d_from_rho(a, r)
      ld_w <- el[d<d_window,.(r2_median=median(r2)),by=SNP]
      ld_w[match(chr_obj$snp_ids,ld_w$SNP),r2_median]
    }, FUN.VALUE = numeric(length(chr_obj$snp_ids)))

    matrix(w, nrow = length(chr_obj$snp_ids), ncol = length(rho),
           dimnames = list(chr_obj$snp_ids, paste0("rho_", rho)))

  }, cores = cores)

  ld_w <- do.call(rbind, by_chr_w)

  if (length(rho) == 1L) return(stats::setNames(as.vector(ld_w), rownames(ld_w)))

  ld_w
}

#' Plot LD-decay results
#'
#' S3 plotting method for objects of class \code{"ld_decay"} produced by
#' \code{compute_LD_decay()}.
#'
#' @param x Object of class \code{"ld_decay"}.
#' @param type Type of plot. One of \code{"summary"}, \code{"chr"},
#'   or \code{"recommendation"}.
#' @param chr Optional chromosome name for \code{type = "chr"}.
#'   Defaults to the first chromosome in \code{x$by_chr}.
#' @param rho Optional target rho value for highlighting in
#'   \code{type = "recommendation"}.
#' @param ... Further graphical arguments passed to low-level plotting functions.
#'
#' @details
#' \describe{
#'   \item{\code{type = "summary"}}{
#'     Shows a multi-panel overview of chromosome-wise LD decay, including
#'     observed and predicted decay rates, slide-window rho coverage, and the
#'     number of informative windows retained per chromosome.
#'   }
#'   \item{\code{type = "chr"}}{
#'     Shows chromosome-specific decay estimates across windows, including
#'     window-wise parameter estimates and fitted LD-decay curves from the
#'     stratified data stored in \code{agg}.
#'   }
#'   \item{\code{type = "recommendation"}}{
#'     Visualizes suggested slide-window sizes in SNP units for the requested
#'     target rho values.freco
#'   }
#' }
#'
#' @return Invisibly returns \code{x}.
#'
#' @method plot ld_decay
#' @export
plot.ld_decay <- function(
    x,
    type = c("summary", "chr", "recommendation"),
    chr = NULL,
    rho = NULL,
    ...
){

  type <- match.arg(type)

  if (!inherits(x, "ld_decay")) {
    stop("`x` must be an object of class 'ld_decay'.")
  }

  ds <- x$decay_sum

  if (is.null(ds) || !nrow(ds)) {
    stop("No `decay_sum` found in `x`.")
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)



  if (type == "summary") {

    par(
      mfrow = c(2, 2),
      mar = c(3, 3, 2, 1),   # tighter margins
      mgp = c(2, 0.7, 0),    # axis label spacing
      tcl = -0.3             # shorter ticks
    )

    ## ------------------------------------------------------------
    ## 1. Observed vs predicted decay rate by chromosome size
    ## ------------------------------------------------------------
    ok <- is.finite(ds$chr_size) & ds$chr_size > 0 &
      is.finite(ds$a) & ds$a > 0

    plot(
      log(ds$chr_size[ok]),
      log(ds$a[ok]),
      xlab = "log chromosome size (bp)",
      ylab = "log observed decay rate (a)",
      main = "Decay rate vs chromosome size",
      pch = 19,
      ...
    )

    if ("a_pred" %in% names(ds)) {
      ord <- order(ds$chr_size)
      lines(
        log(ds$chr_size[ord]),
        log(ds$a_pred[ord]),
        lwd = 2
      )
    }

    if ("Chr" %in% names(ds)) {
      text(
        log(ds$chr_size[ok]),
        log(ds$a[ok]),
        labels = ds$Chr[ok],
        pos = 3,
        cex = 0.8
      )
    }

    ## ------------------------------------------------------------
    ## 2. Rho covered by current slide
    ## ------------------------------------------------------------
    if ("rho_slide_pred" %in% names(ds)) {
      plot(
        seq_len(nrow(ds)),
        ds$rho_slide_pred,
        xaxt = "n",
        xlab = "Chromosome",
        ylab = expression(rho),
        main = "Rho covered by current slide",
        pch = 19,
        ...
      )
      #axis(1, at = seq_len(nrow(ds)), labels = ds$Chr, las = 2)
      axis(1, at = seq_len(nrow(ds)), labels = ds$Chr, las = 2, cex.axis = 0.7)
      abline(h = median(ds$rho_slide_pred, na.rm = TRUE), lty = 2)
    } else {
      plot.new()
      title("Rho covered by current slide")
      text(0.5, 0.5, "No rho_slide_pred available")
    }

    ## ------------------------------------------------------------
    ## 3. Structured windows used per chromosome
    ## ------------------------------------------------------------
    if ("n_w_used" %in% names(ds)) {
      barplot(
        height = ds$n_w_used,
        names.arg = ds$Chr,
        las = 2,
        ylab = "Number of structured windows",
        main = "Informative windows per chromosome",
        ...
      )
    } else {
      plot.new()
      title("Informative windows per chromosome")
      text(0.5, 0.5, "No n_w_used available")
    }

    ## ------------------------------------------------------------
    ## 4. Suggested slide windows across rho targets
    ## ------------------------------------------------------------
    rec_cols <- grep("^slide_snp_rho_", names(ds), value = TRUE)

    if (length(rec_cols) > 0) {
      mat <- as.matrix(ds[, ..rec_cols])
      colnames(mat) <- sub("^slide_snp_rho_", "", rec_cols)

      matplot(
        x = seq_len(nrow(mat)),
        y = mat,
        type = "b",
        pch = 19,
        lty = 1,
        xaxt = "n",
        xlab = "Chromosome",
        ylab = "Suggested slide (SNPs)",
        main = "Recommended slide sizes"
      )
      axis(1, at = seq_len(nrow(ds)), labels = ds$Chr, las = 2)
      legend(
        "topleft",
        legend = paste0("rho=", colnames(mat)),
        lty = 1,
        pch = 19,
        bty = "n"
      )
    } else {
      plot.new()
      title("Recommended slide sizes")
      text(0.5, 0.5, "No slide recommendations available")
    }
  }

  if (type == "recommendation") {

    rec_cols <- grep("^slide_snp_rho_", names(ds), value = TRUE)

    if (!length(rec_cols)) {
      stop("No recommendation columns found in `x$decay_sum`.")
    }

    if (!is.null(rho)) {
      rho_lab <- formatC(rho, format = "f", digits = 2)
      rec_cols <- rec_cols[sub("^slide_snp_rho_", "", rec_cols) == rho_lab]

      if (!length(rec_cols)) {
        stop("Requested `rho` not found among recommendation columns.")
      }
    }

    mat <- as.matrix(ds[, ..rec_cols])
    colnames(mat) <- sub("^slide_snp_rho_", "", rec_cols)

    if (ncol(mat) == 1) {
      plot(
        seq_len(nrow(ds)),
        mat[, 1],
        type = "b",
        pch = 19,
        xaxt = "n",
        xlab = "Chromosome",
        ylab = "Suggested slide (SNPs)",
        main = paste0("Recommended slide size (rho=", colnames(mat)[1], ")"),
        ...
      )
      axis(1, at = seq_len(nrow(ds)), labels = ds$Chr, las = 2)
      abline(h = median(mat[, 1], na.rm = TRUE), lty = 2)
    } else {

      cols <- c("salmon","steelblue","black","firebrick","darkorange")[seq_len(ncol(mat))]

      matplot(
        x = seq_len(nrow(ds)),
        y = mat,
        type = "b",
        pch = 19,
        col = cols,
        lty = 1,
        xaxt = "n",
        xlab = "Chromosome",
        ylab = "Suggested slide (SNPs)",
        main = "Recommended slide sizes",
        ...
      )
      axis(1, at = seq_len(nrow(ds)), labels = ds$Chr, las = 2)
      legend(
        "topleft",
        legend = paste0("rho=", colnames(mat)),
        lty = 1,
        pch = 19,
        col = cols,
        bty = "n"
      )
    }
  }

  if (type == "chr") {

    if (is.null(chr)) {
      chr <- names(x$by_chr)[1]
    }

    if (!chr %in% names(x$by_chr)) {
      stop("`chr` not found in `x$by_chr`.")
    }

    chr_obj <- x$by_chr[[chr]]
    decay   <- chr_obj$decay

    if (is.null(decay) || !nrow(decay)) {
      stop("No chromosome-specific decay data found for `chr`.")
    }


    par(
      mfrow = c(2, 2),
      mar = c(3, 3, 2, 1),   # tighter margins
      mgp = c(2, 0.7, 0),    # axis label spacing
      tcl = -0.3             # shorter ticks
    )

    ## ------------------------------------------------------------
    ## 1. Window-wise a
    ## ------------------------------------------------------------
    if ("a" %in% names(decay)) {
      mid <- rowMeans(decay[, .(start, end)], na.rm = TRUE)

      plot(
        mid,
        decay$a,
        xlab = "Genomic position (bp)",
        ylab = "Decay rate (a)",
        main = paste0(chr, ": window-wise decay rate"),
        pch = 19,
        ...
      )

      chr_sum <- chr_obj$decay_sum
      if (!is.null(chr_sum) && "a" %in% names(chr_sum)) {
        abline(h = chr_sum$a[1], lty = 2, lwd = 2)
      }
    } else {
      plot.new()
      title(paste0(chr, ": window-wise decay rate"))
      text(0.5, 0.5, "No window-wise a estimates")
    }

    ## ------------------------------------------------------------
    ## 2. Window-wise c
    ## ------------------------------------------------------------
    if ("c" %in% names(decay)) {
      mid <- rowMeans(decay[, .(start, end)], na.rm = TRUE)

      plot(
        mid,
        decay$c,
        xlab = "Genomic position (bp)",
        ylab = "Short-range LD (c)",
        main = paste0(chr, ": window-wise short-range LD"),
        pch = 19,
        ...
      )

      chr_sum <- chr_obj$decay_sum
      if (!is.null(chr_sum) && "c" %in% names(chr_sum)) {
        abline(h = chr_sum$c[1], lty = 2, lwd = 2)
      }
    } else {
      plot.new()
      title(paste0(chr, ": window-wise short-range LD"))
      text(0.5, 0.5, "No window-wise c estimates")
    }

    ## ------------------------------------------------------------
    ## 3. Contrast by window
    ## ------------------------------------------------------------
    if (all(c("c", "b") %in% names(decay))) {
      mid <- rowMeans(decay[, .(start, end)], na.rm = TRUE)
      contrast <- decay$c - decay$b

      plot(
        mid,
        contrast,
        xlab = "Genomic position (bp)",
        ylab = expression(c - b),
        main = paste0(chr, ": LD contrast across windows"),
        pch = 19,
        ...
      )
      abline(h = 0.05, lty = 2)
    } else {
      plot.new()
      title(paste0(chr, ": LD contrast across windows"))
      text(0.5, 0.5, "No c/b estimates available")
    }

    ## ------------------------------------------------------------
    ## 4. Stratified decay curves from stored agg data
    ## ------------------------------------------------------------
    agg_list <- decay$agg[!vapply(decay$agg, is.null, logical(1))]
    agg_list <- agg_list[lengths(agg_list) > 0]

    if (length(agg_list) > 0) {

      xlim <- range(unlist(lapply(agg_list, function(z) z$d_mid)), na.rm = TRUE)
      ylim <- range(unlist(lapply(agg_list, function(z) z$r2_q)), na.rm = TRUE)

      plot(
        NA,
        xlim = xlim,
        ylim = ylim,
        xlab = "Distance (bp)",
        ylab = expression(r^2),
        main = paste0(chr, ": fitted decay curves"),
        ...
      )

      for (i in seq_along(agg_list)) {
        agg <- agg_list[[i]]
        points(agg$d_mid, agg$r2_q, pch = 16, cex = 0.6)

        if (all(c("a", "c", "b") %in% names(decay))) {
          ai <- decay$a[match(i, which(!vapply(decay$agg, is.null, logical(1))))]
          ci <- decay$c[match(i, which(!vapply(decay$agg, is.null, logical(1))))]
          bi <- decay$b[match(i, which(!vapply(decay$agg, is.null, logical(1))))]

          if (all(is.finite(c(ai, ci, bi)))) {
            dseq <- seq(min(agg$d_mid), max(agg$d_mid), length.out = 200)
            yhat <- bi + (ci - bi) / (1 + ai * dseq)
            lines(dseq, yhat)
          }
        }
      }

      chr_sum <- chr_obj$decay_sum
      if (!is.null(chr_sum) && all(c("a", "c", "b") %in% names(chr_sum))) {
        dseq <- seq(xlim[1], xlim[2], length.out = 300)
        yhat <- chr_sum$b[1] + (chr_sum$c[1] - chr_sum$b[1]) / (1 + chr_sum$a[1] * dseq)
        lines(dseq, yhat, lwd = 3,col="salmon")
      }

    } else {
      plot.new()
      title(paste0(chr, ": fitted decay curves"))
      text(0.5, 0.5, "No stored stratified fits available")
    }
  }

  invisible(x)
}

