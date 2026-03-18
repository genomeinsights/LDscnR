#' Detect outlier regions from LD-scaled scan
#'
#' @export
detect_or <- function(el,
                      vals,
                      ld_decay,
                      SNP_ids,
                      SNP_chr,
                      sign_th = 0.05,
                      sign_if = c("less", "greater"),
                      rho_d = 0.99,
                      rho_ld = 0.99,
                      l_min = 1,
                      mode = c("per_method","joint"),
                      ret_table = FALSE) {


  if (is.null(colnames(vals)))
    stop("q_vals must have column names.")

  methods <- colnames(vals)

  mode <- match.arg(mode)

  stat_type <- infer_stat_type(methods)

  outlier_fun <- switch(
    sign_if,
    less    = function(x) x < sign_th,
    greater = function(x) x > sign_th
  )

  chrs <- unique(SNP_chr)

  # Outlier rule
  outlier_fun <- switch(
    sign_if[1],
    less    = function(x) x < sign_th,
    greater = function(x) x > sign_th
  )
  #sign_if = "greater"

  # ------------------------------------------------------------
  # Helper: build ORs from SNP set
  # ------------------------------------------------------------
  build_or_from_snps <- function(outliers) {

    if (length(outliers) == 0)
      return(list())

    idx <- which(SNP_ids %in% outliers)

    or_list <- list()
    #ch = "Chr1"
    for (ch in chrs) {
    #cat(ch," -- ")

      out_chr <- SNP_ids[idx][SNP_chr[idx] == ch]

      if (length(out_chr) < l_min)
        next

      # thresholds
      a_chr  <- ld_decay$decay_sum[Chr == ch, a_pred]
      b_chr  <- ld_decay$decay_sum[Chr == ch, b]
      c_chr  <- ld_decay$decay_sum[Chr == ch, c_pred]


      d_th  <- d_from_rho(a_chr, rho = rho_d)
      ld_th <- ld_from_rho(b_chr, c_chr, rho = rho_ld)

      ed <- el[
        d < d_th &
          r2 > ld_th &
          SNP1 %in% out_chr &
          SNP2 %in% out_chr,
        .(SNP1, SNP2)
      ]

      if (nrow(ed) == 0)
        next

      g     <- igraph::graph_from_data_frame(ed, directed = FALSE)
      comps <- igraph::components(g)

      ors_chr <- split(names(comps$membership), comps$membership)
      ors_chr <- ors_chr[vapply(ors_chr, length, integer(1)) >= l_min]

      if (length(ors_chr) > 0)
        or_list <- c(or_list, ors_chr)
    }

    or_list
  }

  # ------------------------------------------------------------
  # PER-METHOD MODE (existing behavior)
  # ------------------------------------------------------------
  if (mode == "per_method") {
    #qs <- unlist(q_vals[,1])

    out <- apply(as.matrix(vals), 2, function(qs) {

      outliers <- SNP_ids[!is.na(qs) & outlier_fun(qs)]
      build_or_from_snps(outliers)
    })

  }

  # ------------------------------------------------------------
  # JOINT MODE (UNION BEFORE CLUSTERING)
  # ------------------------------------------------------------
  if (mode == "joint") {

    # union of outliers across methods
    qs <- apply(vals, 1, min)

    outliers <- SNP_ids[!is.na(qs) & outlier_fun(qs)]

    #union_outliers <- unique(unlist(outlier_matrix))

    joint_or <- build_or_from_snps(outliers)

    # return as single named element
    out <- list(Joint = joint_or)
  }

  out <- structure(
    list(
      ORs     = out,
      params= list (sign_th = sign_th,
                    sign_if = sign_if,
                    rho_d   = rho_d,
                    rho_ld  = rho_ld,
                    l_min   = l_min,
                    mode    = mode)

    ),
    class = "ld_or"
  )

  if (ret_table)
    return(or_to_table(out))

  out
}



#' Draw ORs across parameter space (flat output)
or_draws <- function(el,
                     vals,
                     SNP_ids,
                     SNP_chr,
                     ld_decay,
                     n_draws    = 25,
                     stat_type  = c("q","C"),
                     mode       = c("per_method","joint"),
                     rho_d_lim  = list(min=0.5,max=0.999),
                     rho_ld_lim = list(min=0.9,max=0.999),
                     alpha_lim  = list(min=1.31,max=4),
                     C_lim      = list(min=0,max=0.5),
                     lmin_lim   = list(min=1,max=10),
                     cores      = 1
                     ) {

  stat_type <- match.arg(stat_type)

  if (is.null(colnames(vals)))
    stop("Values must have column names.")
  #mod <- "joint"

  out <- rbindlist(lapply(mode,function(mod){

    methods <- colnames(vals)
    cat("replicate:")
    out <- rbindlist(parallel_apply(seq_len(n_draws),function(i) {
      cat(i, "- ")

      rho_d  <- runif(1, rho_d_lim$min, rho_d_lim$max)
      rho_ld <- runif(1, rho_ld_lim$min, rho_ld_lim$max)

      l_min    <- if (stat_type[1] == "C") 1 else sample(seq(lmin_lim$min, lmin_lim$max), 1)
      alpha    <- if (stat_type[1] == "C") runif(1, C_lim$min, C_lim$max) else 1 / 10^(runif(1, alpha_lim$min, alpha_lim$max))
      sign_if  <- if (stat_type[1] == "C") "greater" else "less"


      if(length(el)==0)
        return(data.table(
          draw_id = i,
          method  = methods,
          rho_d   = rho_d,
          rho_ld  = rho_ld,
          alpha   = alpha,
          l_min   = l_min,
          OR      = replicate(length(methods), list(list()), simplify = FALSE),
          OR_size = 0L
        ))
      #mod = "joint"
      or_obj <- detect_or(
        el        = el,
        vals      = vals,
        SNP_ids   = SNP_ids,
        SNP_chr   = SNP_chr,
        ld_decay  = ld_decay,
        sign_th   = alpha,
        sign_if   = sign_if,
        rho_d     = rho_d,
        rho_ld    = rho_ld,
        l_min     = l_min,
        ret_table = FALSE,
        mode      = mod
      )

      OR_methods <- names(or_obj$ORs)
      if (all(lengths(or_obj$ORs) == 0)) {
        data.table(
          draw_id = i,
          method  = OR_methods,
          rho_d   = rho_d,
          rho_ld  = rho_ld,
          alpha   = alpha,
          l_min   = l_min,
          OR      = replicate(length(OR_methods), list(list()), simplify = FALSE),
          OR_size = 0L
        )
        #meth  <- OR_methods[1]
      }else{
        rbindlist(lapply(OR_methods, function(meth) {

          or_list <- or_obj$ORs[[meth]]
          if (length(or_list) == 0) return(
            data.table(
              draw_id = i,
              method  = meth,
              rho_d   = rho_d,
              rho_ld  = rho_ld,
              alpha   = alpha,
              l_min   = l_min,
              OR      = list(list()),
              OR_size = 0L
            ))

          data.table(
            draw_id = i,
            method  = meth,
            rho_d   = rho_d,
            rho_ld  = rho_ld,
            alpha   = alpha,
            l_min   = l_min,
            OR      = list(or_list),
            OR_size = length(or_list)
          )

        }), fill = TRUE)
      }

    },cores=cores))
    cat("\n")
    return(out)
  }))


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


infer_stat_type <- function(method_names) {

  has_C  <- grepl("^C", method_names) | grepl("_C", method_names)
  has_q  <- !has_C

  if (any(has_C) && any(has_q)) {
    stop("Cannot mix C-based and q-value-based methods in the same call.")
  }

  if (all(has_C)) return("C")
  if (all(has_q)) return("q")

  stop("Unable to infer statistic type.")
}
