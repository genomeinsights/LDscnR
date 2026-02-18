#' Compute chromosome-wise LD structure
#'
#' Computes and stores pairwise LD information for each chromosome in a GDS
#' object. For each chromosome, all SNP pairs within a sliding window are
#' extracted together with their physical distance and LD measure.
#'
#' The resulting object is of class \code{"ld_structure"} and is designed to
#' be reused across downstream analyses such as LD-decay estimation,
#' LD-weight computation, and outlier region detection.
#'
#' @param gds An open GDS object created with \code{SNPRelate}.
#' @param slide_win_ld Integer. Sliding window size (number of SNPs)
#'   used when computing LD via \code{snpgdsLDMat}.
#' @param cores Integer. Number of cores used for LD computation.
#'
#' @return An object of class \code{"ld_structure"} containing:
#' \describe{
#'   \item{by_chr}{Named list with one element per chromosome.}
#'   \item{by_chr[[ch]]$snp_ids}{Vector of SNP identifiers for chromosome \code{ch}.}
#'   \item{by_chr[[ch]]$edges}{A \code{data.table} with pairwise LD information,
#'   including at least \code{SNP1}, \code{SNP2}, \code{r2}, and \code{d}.}
#' }
#'
#' @details
#' The LD structure is computed once and can be reused across multiple
#' LD-scaling and outlier detection analyses, avoiding repeated LD
#' recalculation.
#'
#' @examples
#' \dontrun{
#' ld_struct <- compute_ld_structure(gds, slide_win_ld = 1000)
#' ld_struct
#' }
#'
#' @export
compute_ld_structure <- function(gds,
                                 slide_win_ld = 1000,
                                 cores = 1) {

  ids <- .read_gds_ids(gds)
  chrs <- unique(ids$snp_chr)

  out <- vector("list", length(chrs))
  names(out) <- chrs

  for (ch in chrs) {
    message(ch,"..")
    chr_idx <- which(ids$snp_chr == ch)

    el <- get_el(gds,
                 idx = chr_idx,
                 slide_win_ld = slide_win_ld,
                 cores = cores)

    data.table::setorder(el, d)  # sort by distance once

    out[[ch]] <- list(
      snp_ids = ids$snp_id[chr_idx],
      edges   = el
    )
  }

  structure(
    list(
      by_chr = out
    ),
    class = "ld_structure"
  )
}

.unbiased_ld_w <- function(sub) {

  sub <- data.table::rbindlist(list(
    sub[, .(SNP = SNP1, OTHER = SNP2, pos = pos1, pos_other = pos2, r2)],
    sub[, .(SNP = SNP2, OTHER = SNP1, pos = pos2, pos_other = pos1, r2)]
  ))

  sub[, d := abs(pos_other - pos)]
  sub[, side := ifelse(pos_other > pos, "R", "L")]

  sub[, `:=`(
    maxL = if (any(side == "L")) max(d[side == "L"]) else 0,
    maxR = if (any(side == "R")) max(d[side == "R"]) else 0
  ), by = SNP]

  sub[, S := pmin(maxL, maxR)]

  sub[, long_side := data.table::fifelse(
    maxL > maxR, "L",
    data.table::fifelse(maxR > maxL, "R", NA_character_)
  )]

  sub[, tail := (side == long_side & d > S)]
  sub[is.na(tail), tail := FALSE]

  use <- data.table::rbindlist(list(
    sub[, .(SNP, r2)],
    sub[tail == TRUE, .(SNP, r2)]
  ))

  ld_w2 <- use[, .(
    ld_w = median(r2, na.rm = TRUE),
    ld_w_mean = mean(r2, na.rm = TRUE),
    N = .N
  ), by = SNP]

  return(ld_w2)
}

#' @export
compute_ld_w <- function(ld_struct,
                         decay_obj,
                         rho_w) {

  result <- vector("list", length(ld_struct$by_chr))
  names(result) <- names(ld_struct$by_chr)

  for (ch in names(ld_struct$by_chr)) {

    a_chr <- decay_obj$summary[Chr == ch, a]
    d0_chr <- decay_obj$summary[Chr == ch, d0]
    d_th  <- d_from_rho(a_chr, rho_w,d0 = d0_chr)

    el <- ld_struct$by_chr[[ch]]$edges

    sub <- el[d < d_th]

    med_dt <- data.table::rbindlist(list(
      sub[, .(SNP = SNP1, r2)],
      sub[, .(SNP = SNP2, r2)]
    ))[, .(median_r2 = median(r2, na.rm=TRUE)), by=SNP]

    res <- data.table::data.table(
      SNP = ld_struct$by_chr[[ch]]$snp_ids
    )

    res <- med_dt[res, on="SNP"]

    result[[ch]] <- res$median_r2
  }

  unlist(result)
}

#' @export
print.ld_structure <- function(x, ...) {

  cat("\nLD structure object\n")
  cat("-------------------\n")

  n_chr <- length(x$by_chr)
  cat("Number of chromosomes:", n_chr, "\n")

  total_snps <- sum(vapply(x$by_chr, function(z) length(z$snp_ids), numeric(1)))
  total_edges <- sum(vapply(x$by_chr, function(z) nrow(z$edges), numeric(1)))

  cat("Total SNPs:", total_snps, "\n")
  cat("Total LD pairs:", format(total_edges, big.mark=","), "\n")

  cat("\nPer chromosome summary:\n")

  for (ch in names(x$by_chr)) {
    cat("  ", ch,
        "- SNPs:", length(x$by_chr[[ch]]$snp_ids),
        "| LD pairs:", nrow(x$by_chr[[ch]]$edges), "\n")
  }

  invisible(x)
}

#' @export
plot.ld_structure <- function(x,
                              chr = NULL,
                              n_points = 5000,
                              alpha = 0.5,
                              ...) {

  if (is.null(chr)) {
    chr <- names(x$by_chr)[1]
  }

  if (!chr %in% names(x$by_chr))
    stop("Chromosome not found.")

  el <- x$by_chr[[chr]]$edges

  if (nrow(el) == 0)
    stop("No LD edges available for this chromosome.")

  if (nrow(el) > n_points) {
    set.seed(1)
    el <- el[sample(.N, n_points)]
  }

  plot(el$d,
       el$r2,
       pch = 16,
       cex = 1,
       col = grDevices::adjustcolor("steelblue", alpha.f = alpha),
       xlab = "Distance (bp)",
       ylab = expression(r^2),
       main = paste("LD structure -", chr))

  invisible(x)
}
