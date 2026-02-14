#' Compute consistency score across rho_w and OR draws
#'
#' Calculates SNP-level consistency across all rho_w and OR draws,
#' separately for each method.
#'
#' @param rho_draws_obj Object of class "ld_rho_draws".
#'
#' @return Object of class "ld_consistency".
#' @export
consistency_score <- function(rho_draws_obj) {

  if (!inherits(rho_draws_obj, "ld_rho_draws"))
    stop("Input must be of class 'ld_rho_draws'.")

  rho_draws <- rho_draws_obj$draws
  n_rho     <- rho_draws_obj$n_rho
  n_or      <- rho_draws_obj$n_or_draws

  total_draws <- n_rho * n_or

  if (total_draws == 0)
    stop("No draws found.")

  # Determine method names from first element
  first_or_draw <- rho_draws[[1]]$draws[[1]]
  methods <- names(first_or_draw$ORs)
  if (is.null(methods))
    methods <- paste0("Method_", seq_along(first_or_draw$ORs))

  # Collect all SNPs that ever appear
  all_snps <- unique(unlist(
    lapply(rho_draws, function(rho_obj)
      unlist(
        lapply(rho_obj$draws, function(or_obj)
          unlist(or_obj$ORs)
        )
      )
    )
  ))

  if (length(all_snps) == 0) {
    return(structure(
      list(consistency = list(), total_draws = total_draws),
      class = "ld_consistency"
    ))
  }

  out_list <- vector("list", length(methods))
  names(out_list) <- methods
  #m <- "lfmm_q"
  #rho_obj <- rho_draws[[1]]
  #rho_obj$draws[[1]]$ORs$emx_q
  #rho_obj$draws[[1]]$ORs$emx_F
  for (m in methods) {

    counts <- integer(length(all_snps))
    names(counts) <- all_snps

    for (rho_obj in rho_draws) {
      for (or_obj in rho_obj$draws) {

        snps_in_draw <- unique(unlist(or_obj$ORs[[m]]))

        if (length(snps_in_draw) > 0)
          counts[snps_in_draw] <- counts[snps_in_draw] + 1
      }
    }

    out_list[[m]] <- data.table::data.table(
      SNP = names(counts),
      consistency = counts / total_draws
    )
  }

  structure(
    list(
      consistency = out_list,
      total_draws = total_draws,
      n_rho       = n_rho,
      n_or        = n_or
    ),
    class = "ld_consistency"
  )
}

#' @export
print.ld_consistency <- function(x, ...) {

  cat("\nLD Consistency Score\n")
  cat("--------------------\n")
  cat("Total draws:", x$total_draws, "\n")
  cat("rho_w draws:", x$n_rho, "\n")
  cat("OR draws per rho_w:", x$n_or, "\n\n")
  if (!is.null(x$combined)) {
    cat("Combined C-score available\n")
  }
  for (m in names(x$consistency)) {

    dt <- x$consistency[[m]]

    cat(m, "\n")
    cat("  SNPs with C > 0:", sum(dt$consistency > 0), "\n")
    cat("  Max C:", round(max(dt$consistency), 3), "\n\n")
  }

  invisible(x)
}

#' Add consistency scores to SNP map
#'
#' @param map Data.table with column `marker`
#' @param consistency_obj Object of class "ld_consistency"
#'
#' @return Updated data.table
#' @export
#' @export
combine_consistency <- function(consistency_obj,
                                method = c("mean","min","max","median","geometric"),
                                add_all = TRUE) {

  if (!inherits(consistency_obj, "ld_consistency"))
    stop("Input must be of class 'ld_consistency'.")

  method <- match.arg(method)

  per_method <- consistency_obj$consistency

  renamed <- Map(function(dt, nm) {
    dt2 <- data.table::copy(dt)
    data.table::setnames(dt2, "consistency", nm)
    dt2
  }, per_method, names(per_method))

  merged <- Reduce(function(x, y)
    merge(x, y, by = "SNP", all = TRUE),
    renamed
  )

  method_names <- names(per_method)
  colnames(merged)[-1] <- method_names

  X <- as.matrix(merged[, -1, with = FALSE])

  agg_fun <- switch(
    method,
    mean     = rowMeans,
    min      = function(x) do.call(pmin, c(as.data.frame(x), na.rm = TRUE)),
    max      = function(x) do.call(pmax, c(as.data.frame(x), na.rm = TRUE)),
    median   = function(x) apply(x, 1, median, na.rm = TRUE),
    geometric = function(x) exp(rowMeans(log(pmax(x, 1e-12)), na.rm = TRUE))
  )

  C_joint <- agg_fun(X)

  merged[, C_joint := C_joint]

  if (add_all) {

    merged[, C_mean := rowMeans(X, na.rm = TRUE)]
    merged[, C_min  := do.call(pmin, c(as.data.frame(X), na.rm = TRUE))]
    merged[, C_max  := do.call(pmax, c(as.data.frame(X), na.rm = TRUE))]
  }

  consistency_obj$combined <- merged

  consistency_obj
}

