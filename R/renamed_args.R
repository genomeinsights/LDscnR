## Shared helpers for arguments renamed in 0.9.2 so that one parameter has one
## name across the package. Each function keeps the old name as a trailing
## argument and, when it is supplied, calls this helper, which warns and
## returns the old value for the new name to take.
.renamed_arg <- function(old_value, old_name, new_name) {
  warning(sprintf("`%s` is deprecated. Use `%s` instead.", old_name, new_name),
          call. = FALSE)
  old_value
}

## The decay-relative level a Stage-1 result was built with. Objects from
## before 0.9.2 stored it as `params$rho`.
.stage1_r2_rho <- function(stage1) {
  p <- stage1$params
  if (!is.null(p$r2_rho)) p$r2_rho else p$rho
}
