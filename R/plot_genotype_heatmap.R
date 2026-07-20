#' Plot a genotype dosage heatmap with optional column/row annotations
#'
#' A `ComplexHeatmap::Heatmap()` of a genotype dosage matrix (markers x
#' individuals), with optional colour-bar annotations for both dimensions --
#' e.g. Stage-1 LD-cluster id per marker, population per individual -- and an
#' optional fixed row order. Built to answer the question `plot_pruning_comparison()`
#' can't: whether an apparent LD block is one coherent signal across all
#' individuals, or is instead driven by a subset (e.g. one population being
#' different from everything else).
#'
#' @param GTs Individuals (rows) x markers (cols) dosage matrix (0/1/2),
#'   already restricted to the markers/individuals of interest (e.g. one
#'   [ld_prune_and_eMLG()] group's members).
#' @param polarize Logical; if `TRUE` (default), flips (2 - x) any marker
#'   negatively correlated with the most-complete reference marker via
#'   `polarize_genotypes()` before display, so a real LD block renders as a
#'   visually coherent stripe rather than a checkerboard caused purely by
#'   strand/reference-allele choice. Purely cosmetic -- r^2 (and so any
#'   distance/clustering computed on it) is sign-invariant, so this never
#'   changes `cluster_columns`, only what's drawn.
#' @param cor_threshold Passed to `polarize_genotypes()`'s `cor_threshold`.
#' @param cluster_columns Logical; hierarchically cluster markers (columns)
#'   on 1 - r^2 pairwise LD, computed on the unpolarized genotypes (identical
#'   either way, since r^2 is sign-invariant). Default `TRUE`.
#' @param col_cluster_method Linkage method for column clustering, passed to
#'   `stats::hclust()`. Default `"single"`, consistent with
#'   [ld_complexity_reduction()]'s own Stage 1 approach.
#' @param row_order Optional explicit character vector of `rownames(GTs)`
#'   (e.g. individuals ordered by population), used verbatim and taking
#'   precedence over `cluster_rows`. Useful for keeping one consistent
#'   individual order across many figures.
#' @param cluster_rows Logical; hierarchically cluster individuals (rows) on
#'   Euclidean genotype distance. Defaults to `TRUE` only when `row_order` is
#'   `NULL`.
#' @param col_annotation,row_annotation Optional named vector or factor giving
#'   a categorical grouping to show as a colour bar above (`col_annotation`,
#'   names = markers) or to the left (`row_annotation`, names = individuals)
#'   of the heatmap. Every id in `colnames(GTs)` / `rownames(GTs)` (after any
#'   `row_order` reordering) must have an entry.
#' @param col_annotation_colors,row_annotation_colors Optional named colour
#'   vector (names = levels of the corresponding annotation). Auto-generated
#'   from `grDevices::hcl.colors()` if not supplied.
#' @param col_annotation_name,row_annotation_name Legend/axis label for the
#'   respective annotation. Default `"annotation"`.
#' @param geno_colors Length-3 colour vector for dosage 0/1/2, passed to
#'   `circlize::colorRamp2()`. Default low-mid-high diverging blue/white/red.
#' @param legend_name Heatmap legend title. Defaults to
#'   `"genotype (polarized)"` / `"genotype"` depending on `polarize`.
#' @param title Optional `column_title` for the heatmap.
#' @param out_file Optional path to save the heatmap as a PNG (directory
#'   created if missing). If `NULL` (default), nothing is saved -- draw the
#'   returned object yourself, e.g. `ComplexHeatmap::draw(ht)`.
#' @param width,height PNG dimensions in inches, used only when `out_file` is
#'   given.
#' @param ... Additional arguments passed through to `ComplexHeatmap::Heatmap()`.
#'
#' @return The `Heatmap` object, invisibly. Printing/drawing it (automatic at
#'   the console, or via `ComplexHeatmap::draw()`) renders the figure;
#'   `out_file` (if given) has already been saved as a side effect.
#'
#' @examples
#' \dontrun{
#' grp <- result$groups[group_id == "F3968"]
#' mk <- grp$members[[1]]
#'
#' pop_order <- sample_data[order(Population), Sample_ID]
#' cl_by_marker <- setNames(map_snp$CL_id, map_snp$marker)
#' pop_by_ind <- setNames(sample_data$Population, sample_data$Sample_ID)
#'
#' ht <- plot_genotype_heatmap(
#'   GTs[, mk],
#'   row_order = pop_order,
#'   col_annotation = cl_by_marker[mk], col_annotation_name = "Stage 1 cluster",
#'   row_annotation = pop_by_ind[pop_order], row_annotation_name = "Population",
#'   title = "Chr2: F3968", out_file = "./Figures/Chr2_F3968_heatmap.png"
#' )
#' }
#'
#' @export
plot_genotype_heatmap <- function(GTs, polarize = TRUE, cor_threshold = 0,
                                   cluster_columns = TRUE, col_cluster_method = "single",
                                   row_order = NULL, cluster_rows = is.null(row_order),
                                   col_annotation = NULL, col_annotation_colors = NULL,
                                   col_annotation_name = "annotation",
                                   row_annotation = NULL, row_annotation_colors = NULL,
                                   row_annotation_name = "annotation",
                                   geno_colors = c("#2166AC", "#F7F7F7", "#B2182B"),
                                   legend_name = if (polarize) "genotype (polarized)" else "genotype",
                                   title = NULL, out_file = NULL, width = 10, height = 7, ...) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)) {
    stop(
      "plot_genotype_heatmap() requires the Bioconductor packages ",
      "'ComplexHeatmap' and 'circlize'. Install with:\n",
      '  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")\n',
      '  BiocManager::install(c("ComplexHeatmap", "circlize"))',
      call. = FALSE
    )
  }

  GTs <- as.matrix(GTs)

  if (!is.null(row_order)) {
    stopifnot(all(row_order %in% rownames(GTs)))
    GTs <- GTs[row_order, , drop = FALSE]
  }

  ## column (marker) clustering distance is computed on the *unpolarized*
  ## genotypes -- r^2 is sign-invariant, so this is identical to computing it
  ## after polarization, but doing it first keeps that invariance obvious
  hc_col <- NULL
  if (isTRUE(cluster_columns)) {
    r2 <- suppressWarnings(stats::cor(GTs, use = "pairwise.complete.obs")^2)
    r2[is.na(r2)] <- 0
    hc_col <- stats::hclust(stats::as.dist(1 - r2), method = col_cluster_method)
  }

  GT_plot <- if (isTRUE(polarize)) polarize_genotypes(GTs, cor_threshold = cor_threshold) else GTs

  ## auto-assigns one hcl.colors() palette per annotation level if the caller
  ## didn't supply explicit colours -- keeps this generic (no assumption
  ## about how many levels a given grouping has, unlike a fixed-length
  ## palette constant)
  resolve_annotation <- function(values, ids, colors, arg_name) {
    if (is.null(values)) return(NULL)
    missing <- setdiff(ids, names(values))
    if (length(missing) > 0) {
      stop(
        arg_name, " is missing ", length(missing), " of ", length(ids),
        " required id(s), e.g.: ", paste(utils::head(missing, 3), collapse = ", "),
        call. = FALSE
      )
    }
    v <- values[ids]
    if (is.null(colors)) {
      lv <- if (is.factor(v)) levels(v) else sort(unique(as.character(v)))
      colors <- stats::setNames(grDevices::hcl.colors(length(lv), "Dark 3"), lv)
    }
    list(values = as.character(v), colors = colors)
  }

  col_anno_r <- resolve_annotation(col_annotation, colnames(GT_plot), col_annotation_colors, "col_annotation")
  row_anno_r <- resolve_annotation(row_annotation, rownames(GT_plot), row_annotation_colors, "row_annotation")

  top_annotation <- NULL
  if (!is.null(col_anno_r)) {
    anno_args <- stats::setNames(list(col_anno_r$values), col_annotation_name)
    col_args <- stats::setNames(list(col_anno_r$colors), col_annotation_name)
    top_annotation <- do.call(
      ComplexHeatmap::HeatmapAnnotation,
      c(anno_args, list(col = col_args, annotation_name_side = "left"))
    )
  }

  left_annotation <- NULL
  if (!is.null(row_anno_r)) {
    anno_args <- stats::setNames(list(row_anno_r$values), row_annotation_name)
    col_args <- stats::setNames(list(row_anno_r$colors), row_annotation_name)
    left_annotation <- do.call(
      ComplexHeatmap::rowAnnotation,
      c(anno_args, list(col = col_args, annotation_name_side = "top"))
    )
  }

  geno_col <- circlize::colorRamp2(c(0, 1, 2), geno_colors)

  ## when cluster_rows is FALSE, ComplexHeatmap keeps the matrix's current
  ## row order by default -- which is already the caller's `row_order` (or
  ## GTs's original order), applied above -- so no explicit row_order= is
  ## needed here
  ht <- ComplexHeatmap::Heatmap(
    GT_plot, name = legend_name, col = geno_col,
    cluster_columns = if (!is.null(hc_col)) hc_col else FALSE,
    cluster_rows = cluster_rows,
    show_row_names = FALSE, show_column_names = FALSE,
    top_annotation = top_annotation, left_annotation = left_annotation,
    column_title = title,
    heatmap_legend_param = list(at = c(0, 1, 2), labels = c("0", "1", "2")),
    use_raster = TRUE, raster_quality = 4,
    ...
  )

  if (!is.null(out_file)) {
    dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
    grDevices::png(out_file, width = width, height = height, units = "in", res = 150)
    ComplexHeatmap::draw(ht, merge_legend = TRUE)
    grDevices::dev.off()
    message("Saved: ", out_file)
  }

  invisible(ht)
}
