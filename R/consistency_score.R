#' Compute consistency score across rho_w and OR draws
#'
#' Calculates SNP-level consistency across all rho_w and OR draws,
#' separately for each method.
#'
#' @param rho_draws_obj Object of class "ld_rho_draws".
#'
#' @return Object of class "ld_consistency".
#' @export
consistency_score <- function(or_dt) {

if (!inherits(or_dt, "ld_rho_draws"))
  stop("Input must be class 'ld_rho_draws'.")

dt <- or_dt$draws

if (!all(c("draw_id","method","OR") %in% names(dt)))
  stop("draws table must contain draw_id, method, OR columns.")

# Total draws per method
total_draws_dt <- dt[, .(total_draws = uniqueN(draw_id)), by = method]

# ------------------------------------------------------------
# Expand OR list column safely
# ------------------------------------------------------------
long_dt <- dt[
  ,
  {
    # OR is stored as list of ORs
    # Each OR is a character vector of SNPs
    snps <- unique(unlist(OR))

    if (length(snps) == 0)
      NULL
    else
      .(SNP = snps)
  },
  by = .(draw_id, method)
]

if (nrow(long_dt) == 0) {

  C_dt <- data.table::data.table(
    SNP = character(),
    method = character(),
    C = numeric()
  )

} else {

  # Count in how many draws each SNP appears
  C_dt <- long_dt[
    ,
    .(N = uniqueN(draw_id)),
    by = .(SNP, method)
  ]

  # Attach denominator
  C_dt <- merge(C_dt, total_draws_dt, by = "method")

  C_dt[, C := N / total_draws]

  C_dt[, c("N","total_draws") := NULL]
}

# Store consistency object
C_obj <- structure(
  list(
    consistency = C_dt,
    total_draws = unique(total_draws_dt$total_draws),
    mode        = unique(dt$method)
  ),
  class = "ld_consistency"
)

C_obj
}
# consistency_score <- function(or_dt,combine=TRUE) {
#
#   long_dt <- or_dt$draws[
#     , .(SNP = unlist(OR)),
#     by = .(method)
#   ]
#
#   C_dt <- long_dt[, .N, by = .(SNP,method)]
#
#   total_draws <- or_dt$n_rho*or_dt$n_or_draws
#
#   C_dt[, C := N / total_draws]
#
#
#   C_obj <- structure(
#     list(
#       consistency = C_dt,
#       total_draws = total_draws,
#       n_rho       = or_dt$n_rho,
#       n_or        = or_dt$n_or_draws
#     ),
#     class = "ld_consistency"
#   )
#   if(combine) C_obj <- combine_consistency(C_obj)
#   return(C_obj)
# }

#' @export
add_consistency_to_map <- function(map, consistency_obj) {

  if (!inherits(consistency_obj, "ld_consistency"))
    stop("consistency_obj must be of class 'ld_consistency'.")

  if (!"marker" %in% names(map))
    stop("map must contain column 'marker'.")

  cons_dt <- data.table::copy(consistency_obj$combined)

  data.table::setnames(cons_dt, "SNP", "marker")

  setDT(map)
  setDT(cons_dt)

  out <- cons_dt[map, on = "marker"]
  new_col_order <- unique(c(colnames(map),colnames(cons_dt)))
  return(out[,..new_col_order])

}
#consistency_obj <- C_obj

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
  for (m in x$consistency[,unique(method)]) {

    dt <- x$consistency[method==m]

    cat(m, "\n")
    cat("  SNPs with C > 0:", sum(dt$C > 0), "\n")
    cat("  Max C:", round(max(dt$C), 3), "\n\n")

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
combine_consistency <- function(consistency_obj) {

  if (!inherits(consistency_obj, "ld_consistency"))
    stop("Input must be of class 'ld_consistency'.")

  library(data.table)

  dt <- as.data.table(consistency_obj$consistency)

  if (!all(c("SNP", "method", "C") %in% names(dt)))
    stop("consistency table must contain columns: SNP, method, C")

  # Wide format: one column per method
  wide <- dcast(
    dt,
    SNP ~ method,
    value.var = "C",
    fill = NA_real_
  )


  # Extract method columns
  method_cols <- setdiff(names(wide), "SNP")

  if (length(method_cols) == 0)
    stop("No method columns found.")

  X <- as.matrix(wide[, ..method_cols])

  # Combined summaries
  wide[, C_mean := rowMeans(X, na.rm = TRUE)]
  wide[, C_min  := do.call(pmin, c(as.data.frame(X), na.rm = TRUE))]
  wide[, C_max  := do.call(pmax, c(as.data.frame(X), na.rm = TRUE))]

  setnames(wide, method_cols,paste(method_cols,"C",sep="_"))
  consistency_obj$combined <- wide

  return(consistency_obj)


}

# combine_consistency <- function(consistency_obj) {
#
#   if (!inherits(consistency_obj, "ld_consistency"))
#     stop("Input must be of class 'ld_consistency'.")
#
#   dt <- consistency_obj$consistency
#
#   renamed <- Map(function(dt, nm) {
#     dt2 <- data.table::copy(dt)
#     data.table::setnames(dt2, "consistency", nm)
#     dt2
#   }, per_method, names(per_method))
#
#   merged <- Reduce(function(x, y)
#     merge(x, y, by = "SNP", all = TRUE),
#     renamed
#   )
#
#   X <- as.matrix(merged[, -1, with = FALSE])
#
#   merged[, C_mean := rowMeans(X, na.rm = TRUE)]
#   merged[, C_min  := do.call(pmin, c(as.data.frame(X), na.rm = TRUE))]
#   merged[, C_max  := do.call(pmax, c(as.data.frame(X), na.rm = TRUE))]
#
#   consistency_obj$combined <- merged
#
#   consistency_obj
# }
