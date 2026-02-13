#' @export
compute_ld_structure <- function(gds,
                                 slide_win_ld = 1000,
                                 n_cores = 1) {

  ids <- .read_gds_ids(gds)
  chrs <- unique(ids$snp_chr)

  out <- vector("list", length(chrs))
  names(out) <- chrs

  for (ch in chrs) {
    chr_idx <- which(ids$snp_chr == ch)

    el <- get_el(gds,
                 idx = chr_idx,
                 slide_win_ld = slide_win_ld,
                 n_cores = n_cores)

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
#' @export
compute_ld_w <- function(ld_struct,
                         decay_obj,
                         rho_w) {

  result <- vector("list", length(ld_struct$by_chr))
  names(result) <- names(ld_struct$by_chr)

  for (ch in names(ld_struct$by_chr)) {

    a_chr <- decay_obj$summary[Chr == ch, a]
    d_th  <- d_from_rho(a_chr, rho_w)

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
d_from_rho <- function(a, rho) {
  (1 / a) * (1 / (1 - rho) - 1)
}
