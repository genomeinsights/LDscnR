#' Plot Manhattan-style scan summaries
#'
#' Produces one or more Manhattan-style plots from a SNP map, optionally
#' augmenting the map with consistency scores and detected outlier regions (ORs).
#'
#' The function:
#' \enumerate{
#'   \item computes SNP-level consistency scores from repeated OR draws,
#'   \item adds the selected consistency statistic to the SNP map,
#'   \item detects outlier regions for the selected statistic,
#'   \item prepares chromosome-wise cumulative genomic coordinates,
#'   \item generates one or more Manhattan-style panels.
#' }
#'
#' @param map SNP annotation table containing at least the columns
#'   \code{marker}, \code{Chr}, and \code{Pos}.
#' @param gds Open GDS object containing genotype data.
#' @param ld_decay Object of class \code{"ld_decay"}.
#' @param draws Object containing OR draws; typically an \code{ld_rho_draws}
#'   result with a \code{$draws} component.
#' @param rho_d Relative distance threshold used for OR detection.
#' @param rho_ld Relative LD threshold used for OR detection.
#' @param sign_th Significance threshold used for OR detection.
#' @param mode OR detection mode passed to \code{detect_or()}.
#' @param sign_if Direction of significance; typically \code{"greater"} for
#'   consistency scores and \code{"less"} for q-like statistics.
#' @param l_min Minimum number of SNPs required for an outlier region.
#' @param y_vars Character vector of column names in the augmented map to plot
#'   on the y-axis.
#' @param y_labels Character vector of y-axis labels corresponding to
#'   \code{y_vars}.
#' @param titles Optional panel titles.
#' @param thresholds Optional horizontal threshold values, one per panel.
#' @param col_var Optional column used for point colour, typically an OR label
#'   such as \code{"OR_id"}.
#' @param shape_var Optional column used to highlight special SNP classes
#'   (for example QTNs).
#' @param col_vector Optional vector of colours used for OR categories.
#' @param use_identity Logical; if \code{TRUE}, values in \code{col_var} are
#'   treated as actual colour values and plotted with
#'   \code{ggplot2::scale_color_identity()}. If \code{FALSE}, colours are
#'   assigned as categories, typically via \code{col_vector}.
#'
#' @return A ggplot/patchwork object.
#'
#' @details
#' This function is primarily intended for plotting consistency-based summaries
#' such as \code{Joint_C}, but any numeric column in the augmented map can be
#' used via \code{y_vars}.
#'
#' SNPs not assigned to any detected outlier region are labelled \code{"ns"}
#' and plotted in grey when \code{col_var} is used.
#' @export
plot_manhattan <- function(map,
                           gds,
                           ld_decay,
                           draws,
                           rho_d = 0.99,
                           rho_ld = 0.99,
                           sign_th = 0.05,
                           mode = "joint",
                           sign_if = "greater",
                           l_min = 1,
                           y_vars = c("Joint_C"),
                           y_labels = c("C (Joint analyses)"),
                           titles = NULL,
                           thresholds = c(0.05),
                           col_var = "OR_id",
                           shape_var = NULL,
                           col_vector = NULL,
                           use_identity = FALSE,
                           point_size = 1) {

  map_manh <- add_consistency_to_map(
    map,
    consistency_obj = consistency_score(draws$draws)
  )

  map_manh <- add_ORs(
    gds,
    ld_decay,
    map_manh,
    stat = y_vars[1],
    sign_th = sign_th,
    sign_if = sign_if,
    mode = mode,
    rho_d = rho_d,
    rho_ld = rho_ld,
    l_min = l_min
  )

  vars <- c(y_vars, col_var, shape_var)
  vars <- vars[!is.null(vars)]

  map_manh[, ..vars]

  if (map_manh[, length(which(OR_id != "ns"))] == 0) {
    stop("No outlier regins detected at given significance threshold, aborting")
  }

  layout <- prep_manhattan(
    cbind(
      map_manh[, .(Pos, Chr, marker)],
      map_manh[, ..vars]
    )
  )

  plot_manhattan_gg(
    layout,
    y_vars = y_vars,
    y_labels = y_labels,
    thresholds = thresholds,
    col_var = col_var,
    shape_var = shape_var,
    point_size = point_size,
    titles = titles,
    ncol = 1,
    col_vector = col_vector,
    use_identity = use_identity
  )
}

#' Add outlier-region labels to a SNP map
#'
#' Detects outlier regions from a selected SNP-level statistic and merges the
#' resulting region labels into a SNP map.
#'
#' @param gds Open GDS object.
#' @param ld_decay Object of class \code{"ld_decay"}.
#' @param map SNP annotation table containing at least \code{marker},
#'   \code{Chr}, and the selected statistic column.
#' @param stat Name of the statistic column in \code{map} used for OR detection.
#' @param sign_th Significance threshold used to classify outlier SNPs.
#' @param sign_if Direction of significance; either \code{"greater"} or
#'   \code{"less"}.
#' @param mode OR detection mode passed to \code{detect_or()}.
#' @param rho_d Relative distance threshold used for OR detection.
#' @param rho_ld Relative LD threshold used for OR detection.
#' @param l_min Minimum number of SNPs required for an outlier region.
#'
#' @return An updated \code{data.table} with an added \code{OR_id} column.
#'
#' @details
#' SNPs not assigned to any outlier region are labelled \code{"ns"}.
#'
#' The function first restricts LD estimation to SNPs exceeding the supplied
#' threshold, computes the corresponding LD edge list, detects outlier regions,
#' and then merges the resulting region assignments back into the full map.
#'
#' @export
add_ORs <- function(gds, ld_decay, map, stat="Joint_C", sign_th,sign_if="greater" ,mode,  rho_d, rho_ld,l_min=1){
  idx <- which(map[,..stat]>sign_th)

  el  <- get_el(gds, idx, slide_win_ld = -1,by_chr = TRUE)

  ORs_tbl <- detect_or(el,
                       vals = map[,..stat],
                       ld_decay = ld_decay,
                       SNP_ids  = map[,marker],
                       SNP_chr  = map[,Chr],
                       sign_th  = sign_th,
                       sign_if  = sign_if,
                       rho_d    = rho_d,
                       rho_ld   = rho_ld,
                       mode     = mode,
                       l_min    = l_min,
                       ret_table = TRUE)

  if(!is.null(map$OR_id)) map[,OR_id:=NULL]

  map <- merge(map,
               ORs_tbl,
               by.x  = "marker",
               by.y  = "SNP",
               all.x = TRUE)

  map <- map[order(as.numeric(gsub("Chr", "", Chr)), Pos)]
  map[is.na(OR_id),OR_id:="ns"]
  return(map)
}

#' Prepare cumulative coordinates for Manhattan plotting
#'
#' Constructs cumulative genomic coordinates and background annotation needed
#' for Manhattan-style plotting across chromosomes.
#'
#' @param data_manh Data frame or \code{data.table} containing at least the
#'   columns \code{Chr} and \code{Pos}.
#' @param spacer Non-negative spacing added between chromosomes on the cumulative
#'   x-axis.
#' @param chr_cols Alternating background colours used for chromosome bands.
#'
#' @return An object of class \code{"manhattan_layout"} containing:
#' \describe{
#'   \item{data}{Input data with cumulative genomic coordinates \code{BPcum}.}
#'   \item{axis}{Chromosome axis label positions.}
#'   \item{rect}{Chromosome background rectangles.}
#' }
#'
#' @details
#' Chromosomes are ordered by the numeric suffix extracted from labels such as
#' \code{"Chr1"}, \code{"Chr2"}, etc.
#'
#' @export
prep_manhattan <- function(data_manh,
                           spacer = 0,
                           chr_cols = c("white", "grey50")) {

  stopifnot(all(c("Chr", "Pos") %in% names(data_manh)))

  data_manh <- as.data.table(data_manh)

  ## extract numeric chromosome index once
  data_manh[, chr_idx := as.integer(sub("^Chr", "", Chr))]
  setorder(data_manh, chr_idx, Pos)

  ## fix chromosome factor order
  chr_levels <- unique(data_manh$Chr)
  data_manh[, Chr := factor(Chr, levels = chr_levels)]

  ## chromosome lengths + cumulative offsets
  chr_dt <- data_manh[, .(chr_len = max(Pos) + spacer), by = Chr]
  chr_dt[, offset := cumsum(chr_len) - chr_len]

  ## merge offsets back
  don <- merge(data_manh, chr_dt[, .(Chr, offset)], by = "Chr")
  don[, BPcum := Pos + offset]
  setorder(don, Chr, Pos)

  ## axis centers
  axisdf <- don[, .(center = (min(BPcum) + max(BPcum)) / 2), by = Chr]

  ## background rectangles
  rect_data <- don[, .(
    x1 = min(BPcum),
    x2 = max(BPcum),
    y1 = -Inf,
    y2 = Inf
  ), by = Chr]

  rect_data[, col := chr_cols[(seq_len(.N) - 1) %% length(chr_cols) + 1]]

  structure(
    list(
      data  = don,
      axis  = axisdf,
      rect  = rect_data
    ),
    class = "manhattan_layout"
  )
}

#' Render Manhattan-style panels
#'
#' Generates one or more Manhattan-style panels from a
#' \code{manhattan_layout} object.
#'
#' @param layout Object returned by \code{prep_manhattan()}.
#' @param y_vars Character vector of numeric variables to plot.
#' @param y_labels Character vector of y-axis labels corresponding to
#'   \code{y_vars}.
#' @param thresholds Optional numeric vector of horizontal threshold values,
#'   one per panel.
#' @param col_var Optional column name used for point colour.
#' @param shape_var Optional column name used to highlight special SNP classes.
#' @param titles Optional panel titles.
#' @param point_size Point size.
#' @param ncol Optional number of columns in the multi-panel layout.
#' @param col_vector Optional vector of colours for categorical point colouring.
#'
#' @return A ggplot/patchwork object.
#'
#' @details
#' Chromosome backgrounds are shown as alternating shaded rectangles. SNPs with
#' \code{col_var == "ns"} are plotted in grey when a colour variable is used.
#'
#' If \code{shape_var} is supplied, entries equal to \code{"QTN"} are highlighted
#' with larger cross-shaped markers.
#'
#' @export
plot_manhattan_gg <- function(layout,
                              y_vars,
                              y_labels,
                              thresholds = NULL,
                              col_var = NULL,
                              shape_var = NULL,
                              titles = NULL,
                              point_size = 1,
                              ncol = NULL,
                              col_vector=NULL,
                              use_identity=FALSE) {

  if(is.null(col_vector)) col_vector <- c("#B2DF8A", "#FFD92F", "firebrick", "#33A02C", "#7FC97F", "#CAB2D6",
                                          "#FB8072", "grey30", "#E6AB02", "#FDC086", "steelblue", "#1F78B4",
                                          "#FB9A99", "#1B9E77", "#BC80BD", "#E31A1C", "#7570B3", "#A6761D",
                                          "#A6CEE3", "salmon", "#FFFF33", "forestgreen", "#FDCDAC", "#BF5B17",
                                          "#A6761D", "#FBB4AE", "#4DAF4A", "#B3E2CD", "#FDDAEC", "#BEBADA",
                                          "#FFF2AE", "#1F78B4", "#66C2A5", "#F0027F", "#E6AB02", "#E78AC3",
                                          "#FF7F00", "#8DA0CB", "#6A3D9A", "#B15928", "#E41A1C")

  if (!inherits(layout, "manhattan_layout"))
    stop("layout must be from prep_manhattan().")


  don     <- data.table::copy(layout$data)
  axisdf  <- layout$axis
  rect_dt <- layout$rect

  n_panels <- length(y_vars)

  if (!is.null(thresholds) && length(thresholds) != n_panels)
    stop("thresholds must match length of y_vars")

  ## make sure we have enough colors
  tmp <- unique(unlist(don[,..col_var]))
  tmp <- tmp[tmp!="ns"]
  col_vector <- rep(col_vector,ceiling(length(tmp)/length(col_vector)))

  #i  <- 1
  plots <- lapply(seq_len(n_panels), function(i) {

    don[, yval := get(y_vars[i])]
    don[is.na(yval), yval := 0]

    p <- ggplot2::ggplot() +

      # chromosome background
      ggplot2::geom_rect(
        data = rect_dt,
        ggplot2::aes(xmin = x1, xmax = x2, ymin = -Inf, ymax = Inf),
        fill = rect_dt$col,
        alpha = 0.3,
        inherit.aes = FALSE
      )

    # --- Non-OR points ---
    if (!is.null(col_var)) {
      p <- p + ggplot2::aes(color = .data[[col_var]])

      if (isTRUE(use_identity)) {
        p <- p + ggplot2::scale_color_identity()
      } else if (!is.null(col_vector)) {
        p <- p + ggplot2::scale_color_manual(values = col_vector)
      }
    }

    # --- OR coloured points ---
    if (!is.null(col_var)) {
      p <- p +
        ggplot2::geom_point(
          data = don[get(col_var) != "ns"],
          ggplot2::aes(BPcum, yval, colour = get(col_var)),
          size = point_size
        )
    } else {
      p <- p +
        ggplot2::geom_point(
          data = don,
          ggplot2::aes(BPcum, yval),
          size = point_size,
          colour = "grey50"
        )
    }

    # --- QTN markers ---
    if (!is.null(shape_var)) {
      p <- p +
        ggplot2::geom_point(
          data = don[get(shape_var) == "QTN" & get(col_var) != "ns"],
          ggplot2::aes(BPcum, yval,colour = get(col_var)),
          shape = 3,
          size = point_size * 3
        ) +
        ggplot2::geom_point(
          data = don[get(shape_var) == "QTN" & get(col_var) == "ns"],
          ggplot2::aes(BPcum, yval),
          shape = 3,
          size = point_size * 3,
          col="black"
        )
    }

    # --- Threshold ---
    if (!is.null(thresholds) && !is.na(thresholds[i])) {
      p <- p +
        ggplot2::geom_hline(
          yintercept = thresholds[i],
          linetype = 2,
          colour = "grey30"
        )
    }

    p +
      ggplot2::scale_x_continuous(
        breaks = axisdf$center,
        labels = axisdf$Chr,
        expand = ggplot2::expansion(mult = c(0.01, 0.01))
      ) +
      ggplot2::labs(
        x = NULL,
        y = y_labels[i]
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 90),
        legend.position = "none",
        aspect.ratio = 0.25
      ) +
      ggplot2::scale_color_manual(values=col_vector)+
      ggplot2::ggtitle(titles[i])
  })

  if(is.null(ncol)){
    patchwork::wrap_plots(plots)
  }else{
    patchwork::wrap_plots(plots,ncol=ncol)
  }
}
