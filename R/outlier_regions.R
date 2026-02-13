#' Detect outlier regions from LD-scaled scan
#'
#' @export
detect_or <- function(scan_obj,
                      ld_struct,
                      decay_obj,
                      alpha,
                      rho_d,
                      rho_ld,
                      l_min = 2) {

  if (!inherits(scan_obj, "ld_scan"))
    stop("scan_obj must be of class 'ld_scan'.")

  ids <- scan_obj$SNP_ids
  #res <- scan_obj$result[[1]]
  out <- lapply(scan_obj$result,function(res){
    q_vals <- res$q_prime

    q_vec <- if (is.matrix(q_vals)) q_vals else q_vals
    outliers <- ids[q_vec < alpha]

    if (length(outliers) == 0) {
      out[[i]] <- list()
      next
    }

    or_list <- list()

    # loop chromosomes
    # ch = "Chr1"
    for (ch in names(ld_struct$by_chr)) {

      snp_ids_chr <- ld_struct$by_chr[[ch]]$snp_ids
      el_chr      <- ld_struct$by_chr[[ch]]$edges

      # restrict to outliers on this chr
      out_chr <- intersect(outliers, snp_ids_chr)
      if (length(out_chr) < l_min) next

      # compute thresholds
      a_chr <- decay_obj$summary[Chr == ch, a]
      b_chr <- decay_obj$summary[Chr == ch, b]

      d_th  <- d_from_rho(a_chr, rho = rho_d)
      ld_th <- b_chr + (1 - b_chr) * (1 - rho_ld)

      # filter edges
      ed <- el_chr[
        d < d_th &
          r2 > ld_th &
          SNP1 %in% out_chr &
          SNP2 %in% out_chr,
        .(SNP1, SNP2)
      ]

      if (nrow(ed) == 0) next

      g <- igraph::graph_from_data_frame(ed, directed = FALSE)
      comps <- igraph::components(g)

      ors_chr <- split(names(comps$membership), comps$membership)
      ors_chr <- ors_chr[vapply(ors_chr, length, integer(1)) >= l_min]

      if (length(ors_chr) > 0)
        or_list <- c(or_list, ors_chr)
    }

    or_list
  })




  structure(
    list(
      ORs = out,
      alpha = alpha,
      rho_d = rho_d,
      rho_ld = rho_ld,
      l_min = l_min
    ),
    class = "ld_or"
  )
}

or_draws <- function(scan_obj,
                     ld_struct,
                     decay_obj,
                     n_draws = 25,
                     rho_d_lim=list(min=0.5,max=0.999),
                     rho_ld_lim=list(min=0.9,max=0.999),
                     alpha_lim=list(min=1.31,max=4),
                     lmin_lim=list(min=1,max=10)) {

  out <- vector("list", n_draws)

  for (i in seq_len(n_draws)) {

    rho_d  <- runif(1, rho_OR_lim$min, rho_OR_lim$max)
    rho_ld  <- runif(1, rho_OR_lim$min, rho_OR_lim$max)


    alpha  <- 1 / 10^(runif(1, alpha_lim$min, alpha_lim$max))
    l_min  <- sample(seq(lmin_lim$min, lmin_lim$max), 1)

    out[[i]] <- detect_or(
      scan_obj,
      ld_struct,
      decay_obj,
      alpha = alpha,
      rho_d = rho_d,
      rho_ld = rho_ld,
      l_min = l_min
    )
  }

  structure(
    list(draws = out),
    class = "ld_or_draws"
  )
}

