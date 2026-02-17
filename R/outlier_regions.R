#' Detect outlier regions from LD-scaled scan
#'
#' @export
detect_or <- function(q_vals,
                      ld_struct,
                      decay_obj,
                      sign_th,
                      sign_if = c("less", "greater"),
                      rho_d,
                      rho_ld,
                      l_min = 2,
                      ret_table=FALSE) {

  #sign_if <- match.arg(sign_if,c("less", "greater"))

  if (is.null(colnames(q_vals)))
    stop("q_vals must have column names.")

  if (is.null(names(ld_struct$by_chr)))
    stop("ld_struct must contain chromosome-wise LD structure.")

  ids <- unlist(sapply(ld_struct$by_chr,function(x)x$snp_ids))

  # ------------------------------------------------------------
  # Ensure either all q-values or all C-values
  # ------------------------------------------------------------
  is_C <- grepl("^C", colnames(q_vals))
  value_types <- table(is_C) / ncol(q_vals)

  if (length(value_types) > 1)
    stop("Must be either only q-values or only C-values.")

  # ------------------------------------------------------------
  # Define outlier rule
  # ------------------------------------------------------------
  # "less" is default, assumes q-values
  outlier_fun <- switch(
    sign_if[1],
    less    = function(x) x < sign_th,
    greater = function(x) x > sign_th
  )

  # ------------------------------------------------------------
  # Detect ORs per column
  # ------------------------------------------------------------


  out <- apply(q_vals, 2, function(qs) {

    outliers <- as.vector(na.omit(ids[!is.na(qs)][outlier_fun(qs[!is.na(qs)])]))

    if (length(outliers) == 0)
      return(list())

    or_list <- list()

    # Loop chromosomes
    for (ch in names(ld_struct$by_chr)) {

      snp_ids_chr <- ld_struct$by_chr[[ch]]$snp_ids
      el_chr      <- ld_struct$by_chr[[ch]]$edges

      # restrict to outliers on this chr
      out_chr <- intersect(outliers, snp_ids_chr)
      if (length(out_chr) < l_min)
        next

      # compute thresholds
      a_chr <- decay_obj$summary[Chr == ch, a]
      b_chr <- decay_obj$summary[Chr == ch, b]
      c_chr <- decay_obj$summary[Chr == ch, c]
      d0_chr <- decay_obj$summary[Chr == ch, d0]

      d_th  <- d_from_rho(a_chr, rho = rho_d,d0 = d0_chr)
      ld_th <- ld_from_rho(b_chr, c_chr, rho = rho_d)

      # filter edges
      ed <- el_chr[
        d < d_th &
          r2 > ld_th &
          SNP1 %in% out_chr &
          SNP2 %in% out_chr,
        .(SNP1, SNP2)
      ]

      if (nrow(ed) == 0)
        next

      g <- igraph::graph_from_data_frame(ed, directed = FALSE)
      comps <- igraph::components(g)

      ors_chr <- split(names(comps$membership), comps$membership)
      ors_chr <- ors_chr[vapply(ors_chr, length, integer(1)) >= l_min]

      if (length(ors_chr) > 0)
        or_list <- c(or_list, ors_chr)
    }

    or_list
  })

  out <- structure(
    list(
      ORs     = out,
      sign_th = sign_th,
      sign_if = sign_if,
      rho_d   = rho_d,
      rho_ld  = rho_ld,
      l_min   = l_min
    ),
    class = "ld_or"
  )

  ifelse(ret_table,return(or_to_table(out)),return(out))


}


or_draws <- function(q_vals,
                     ld_struct,
                     decay_obj,
                     n_draws = 25,
                     rho_d_lim=list(min=0.5,max=0.999),
                     rho_ld_lim=list(min=0.9,max=0.999),
                     alpha_lim=list(min=1.31,max=4),
                     lmin_lim=list(min=1,max=10)) {

  out <- vector("list", n_draws)
  #i <- 13
  for (i in seq_len(n_draws)) {

    rho_d  <- runif(1, rho_d_lim$min, rho_d_lim$max)
    rho_ld  <- runif(1, rho_ld_lim$min, rho_ld_lim$max)


    alpha  <- 1 / 10^(runif(1, alpha_lim$min, alpha_lim$max))
    l_min  <- sample(seq(lmin_lim$min, lmin_lim$max), 1)

    out[[i]] <- detect_or(
      q_vals   = q_vals,
      ld_struct = ld_struct,
      decay_obj = decay_obj,
      sign_th   = alpha,
      rho_d     = rho_d,
      rho_ld    = rho_ld,
      l_min     = l_min
    )
  }

  structure(
    list(draws = out),
    class = "ld_or_draws"
  )
}


#' Convert ld_or object to long table
#'
#' @export
or_to_table <- function(or_obj) {

  if (!inherits(or_obj, "ld_or"))
    stop("Input must be class 'ld_or'.")

  res <- list()

  for (method in names(or_obj$ORs)) {

    or_list <- or_obj$ORs[[method]]
    if (length(or_list) == 0) next

    for (i in seq_along(or_list)) {

      snps <- or_list[[i]]

      # Extract chromosome from first SNP
      chr <- sub(":.*", "", snps[1])

      or_name <- paste0(chr, "_OR", i)

      res[[length(res) + 1]] <- data.table::data.table(
        SNP    = snps,
        method = method,
        OR_id  = or_name
      )
    }
  }

  if (length(res) == 0)
    return(data.table::data.table(SNP=character(), method=character(), OR_id=character()))

  data.table::rbindlist(res)
}
