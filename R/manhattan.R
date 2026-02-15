col_vector <- c("#B2DF8A", "#FFD92F", "firebrick", "#33A02C", "#7FC97F", "#CAB2D6",
                "#FB8072", "grey30", "#E6AB02", "#FDC086", "steelblue", "#1F78B4",
                "#FB9A99", "#1B9E77", "#BC80BD", "#E31A1C", "#7570B3", "#A6761D",
                "#A6CEE3", "salmon", "#FFFF33", "forestgreen", "#FDCDAC", "#BF5B17",
                "#A6761D", "#FBB4AE", "#4DAF4A", "#B3E2CD", "#FDDAEC", "#BEBADA",
                "#FFF2AE", "#1F78B4", "#66C2A5", "#F0027F", "#E6AB02", "#E78AC3",
                "#FF7F00", "#8DA0CB", "#6A3D9A", "#B15928", "#E41A1C")


prep_manhattan <- function(data_manh,
                           spacer = 5e6,
                           chr_cols = c("bisque", "white")) {

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
                              y,
                              point_size = 0.4,
                              alpha = 0.8,
                              color_vect=NULL,
                              shape_vect=NULL) {

  stopifnot(inherits(layout, "manhattan_layout"))

  don      <- copy(layout$data)
  axisdf  <- layout$axis
  rect_dt <- layout$rect

  ## y aesthetic
  don[, yval := get(y)]

  ## significance coloring
  if (is.null(color_vect)) {
    don[, col_custom := "none"]
  }

  if (is.null(shape_vect)) {
    don[, shape_custom := "none"]
  }

  ggplot() +
    ## chromosome background
    geom_rect(
      data = rect_dt,
      aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
      fill = rect_dt$col,
      alpha = 0.3,
      inherit.aes = FALSE
    ) +
    ## points
    geom_point(
      data = don,
      aes(BPcum, yval, color = col_custom,shape=shape_custom),
      size = point_size,
      alpha = alpha
    ) +
    scale_color_manual(values =c("black",rep(col_vector,3)),guide="none")+
    scale_shape_manual(values =c(20,3),guide="none")+
    scale_x_continuous(
      breaks = axisdf$center,
      labels = axisdf$Chr,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(x = NULL, y = y) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 0),
      aspect.ratio = 0.25
    )
}

