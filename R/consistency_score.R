#' Compute SNP-level consistency across OR draws
#'
#' Calculates SNP-level consistency scores from repeated outlier-region (OR)
#' detection draws.
#'
#' For each method, the consistency score of a SNP is defined as the proportion
#' of retained draws in which that SNP appears in at least one detected outlier
#' region.
#'
#' Extremely large OR solutions are excluded before scoring using an
#' upper-tail outlier filter based on:
#' \deqn{Q3 + 1.5 \times IQR}
#' applied to \code{OR_size} across draws with at least one detected region.
#'
#' @param draws A \code{data.table} of OR draws, typically returned by
#'   \code{ld_rho_draws()} or \code{or_draws()}, containing at least the columns
#'   \code{method}, \code{OR}, and \code{OR_size}.
#' @param combine Logical; currently unused.
#'
#' @return An object of class \code{"ld_consistency"} containing:
#' \describe{
#'   \item{consistency}{A long-format \code{data.table} with columns
#'   \code{SNP}, \code{method}, \code{N}, and \code{C}, where \code{C} is the
#'   SNP-level consistency score.}
#'   \item{total_draws}{Number of draws per method used to compute consistency.}
#' }
#'
#' @details
#' The score is computed separately for each method as:
#' \deqn{C = N / n_{draws}}
#' where \eqn{N} is the number of retained draws in which the SNP occurs in an
#' outlier region and \eqn{n_{draws}} is the total number of draws for that
#' method.
#'
#' This measure is intended to highlight SNPs that recur consistently across
#' parameter settings, rather than SNPs appearing only in isolated OR solutions.
#'
#' @export
consistency_score <- function(draws,combine=TRUE) {

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

#' Add per-method consistency scores to a SNP map
#'
#' Merges SNP-level consistency scores into a SNP annotation table.
#'
#' Each method-specific consistency score is added as a new column with suffix
#' \code{"_C"}.
#'
#' @param map A data.frame or \code{data.table} containing a column named
#'   \code{marker}.
#' @param consistency_obj Object of class \code{"ld_consistency"}.
#'
#' @return A \code{data.table} containing the original map columns plus one
#' additional column per method-specific consistency score.
#'
#' @details
#' Consistency scores are reshaped to wide format before merging. The merge is
#' performed as a left join on \code{marker}, preserving all rows from
#' \code{map}.
#'
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


#' Combine consistency scores across methods
#'
#' Computes cross-method summary statistics from an object of class
#' \code{"ld_consistency"}.
#'
#' Method-specific consistency scores are reshaped to wide format and combined
#' into summary measures across methods.
#'
#' @param consistency_obj Object of class \code{"ld_consistency"}.
#'
#' @return The input \code{ld_consistency} object with an added
#' \code{combined} component containing:
#' \describe{
#'   \item{SNP}{SNP identifier.}
#'   \item{<method>_C}{Method-specific consistency scores in wide format.}
#'   \item{C_mean}{Mean consistency across methods.}
#'   \item{C_min}{Minimum consistency across methods.}
#'   \item{C_max}{Maximum consistency across methods.}
#' }
#'
#' @details
#' Missing method-specific values are retained as \code{NA}. Summary statistics
#' are computed across the available method-specific consistency scores for each
#' SNP.
#'
#' This function is useful when consistency should be interpreted jointly across
#' multiple methods rather than method by method.
#'
#' @export
combine_consistency <- function(consistency_obj) {

  if (!inherits(consistency_obj, "ld_consistency"))
    stop("Input must be of class 'ld_consistency'.")

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
