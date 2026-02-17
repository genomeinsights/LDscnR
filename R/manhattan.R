col_vector <- c("#B2DF8A", "#FFD92F", "firebrick", "#33A02C", "#7FC97F", "#CAB2D6",
                "#FB8072", "grey30", "#E6AB02", "#FDC086", "steelblue", "#1F78B4",
                "#FB9A99", "#1B9E77", "#BC80BD", "#E31A1C", "#7570B3", "#A6761D",
                "#A6CEE3", "salmon", "#FFFF33", "forestgreen", "#FDCDAC", "#BF5B17",
                "#A6761D", "#FBB4AE", "#4DAF4A", "#B3E2CD", "#FDDAEC", "#BEBADA",
                "#FFF2AE", "#1F78B4", "#66C2A5", "#F0027F", "#E6AB02", "#E78AC3",
                "#FF7F00", "#8DA0CB", "#6A3D9A", "#B15928", "#E41A1C")


prep_manhattan <- function(data_manh,
                           spacer = 0,
                           chr_cols = c("white", "grey50")) {

  stopifnot(all(c("Chr", "bp") %in% names(data_manh)))

  data_manh <- as.data.table(data_manh)

  ## extract numeric chromosome index once
  data_manh[, chr_idx := as.integer(sub("^Chr", "", Chr))]
  setorder(data_manh, chr_idx, bp)

  ## fix chromosome factor order
  chr_levels <- unique(data_manh$Chr)
  data_manh[, Chr := factor(Chr, levels = chr_levels)]

  ## chromosome lengths + cumulative offsets
  chr_dt <- data_manh[, .(chr_len = max(bp) + spacer), by = Chr]
  chr_dt[, offset := cumsum(chr_len) - chr_len]

  ## merge offsets back
  don <- merge(data_manh, chr_dt[, .(Chr, offset)], by = "Chr")
  don[, BPcum := bp + offset]
  setorder(don, Chr, bp)

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

plot_manhattan_gg <- function(layout,
                            y_vars,
                            y_labels,
                            thresholds = NULL,
                            or_var = NULL,
                            type_var = NULL,
                            point_size = 1,
                            ncol = NULL) {

  if (!inherits(layout, "manhattan_layout"))
    stop("layout must be from prep_manhattan().")

  don     <- data.table::copy(layout$data)
  axisdf  <- layout$axis
  rect_dt <- layout$rect

  n_panels <- length(y_vars)

  if (!is.null(thresholds) && length(thresholds) != n_panels)
    stop("thresholds must match length of y_vars")

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
    if (!is.null(or_var)) {
      p <- p +
        ggplot2::geom_point(
          data = don[get(or_var) == "ns"],
          ggplot2::aes(BPcum, yval),
          size = point_size,
          colour = "grey50"
        )
    }

    # --- OR coloured points ---
    if (!is.null(or_var)) {
      p <- p +
        ggplot2::geom_point(
          data = don[get(or_var) != "ns"],
          ggplot2::aes(BPcum, yval, colour = get(or_var)),
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
    if (!is.null(type_var)) {
      p <- p +
        ggplot2::geom_point(
          data = don[get(type_var) == "QTN"],
          ggplot2::aes(BPcum, yval),
          shape = 3,
          size = point_size * 3,
          colour = "black"
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
      )
  })

  if(is.null(ncol)){
    patchwork::wrap_plots(plots)
  }else{
    patchwork::wrap_plots(plots,ncol=ncol)
  }
}
