#' Compute consistency score across rho_w and OR draws
#'
#' Calculates SNP-level consistency across all rho_w and OR draws,
#' separately for each method.
#'
#' @param rho_draws_obj Object of class "ld_rho_draws".
#'
#' @return Object of class "ld_consistency".
#' @export
consistency_score <- function(draws,combine=TRUE) {


  size_th <- draws[OR_size>0,quantile(OR_size,0.75)+1.5*IQR(OR_size)]

  long_dt <- draws[OR_size<size_th
    , .(SNP = unlist(OR)),
    by = .(method)
  ]

  n_draws <- draws[method==method[1],.N]
  C_dt <- long_dt[, .N, by = .(SNP,method)]

  #total_draws <- or_dt$n_rho*or_dt$n_or_draws

  C_dt[, C := N / n_draws]

  C_obj <- structure(
    list(
      consistency = C_dt,
      total_draws = n_draws
    ),
    class = "ld_consistency"
  )

  return(C_obj)
}

# consistency_obj <- consistency
#' @export
add_consistency_to_map <- function(map, consistency_obj) {

  if (!inherits(consistency_obj, "ld_consistency"))
    stop("consistency_obj must be of class 'ld_consistency'.")

  if (!"marker" %in% names(map))
    stop("map must contain column 'marker'.")

  setDT(map)


  # Combine on the fly
  cons_raw <- data.table::copy(consistency_obj$consistency)

  cons_raw[,method:=paste(method,"C",sep="_")]
  # Wide format
  cons_dt <- data.table::dcast(
    cons_raw,
    SNP ~ method,
    value.var = "C"
  )

  # Rename SNP → marker
  setnames(cons_dt, "SNP", "marker")

  # Merge (left join preserving map order)
  out <- merge(map, cons_dt, by = "marker", all.x = TRUE, sort = FALSE)

  # Preserve original map column order first
  new_cols <- setdiff(names(out), names(map))
  out <- out[, c(names(map), new_cols), with = FALSE]

  out
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
