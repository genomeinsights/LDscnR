## Faceted Manhattan plotting for LD-aware outlier-region analyses. A single,
## self-contained helper (migrated/streamlined from the LDscnR-paper
## plot_OR_manhattan) that covers every panel the workflow needs: a value on the
## y-axis (C-score, -log10(q), ld_w, an association statistic) coloured either by
## discrete outlier region / group, or by a continuous per-marker column.

## internal: resolve a per-marker numeric vector from either a column name in
## `map`, a named vector (names = markers), or a full-length aligned vector.
.resolve_marker_value <- function(x, markers, map) {
  if (is.character(x) && length(x) == 1L && x %in% names(map))
    return(as.numeric(map[[x]])[match(markers, map$marker)])
  if (!is.null(names(x))) return(as.numeric(x[markers]))
  if (length(x) == length(markers)) return(as.numeric(x))
  stop("Cannot resolve a per-marker vector: pass a column name in `map`, ",
       "a named vector (names = markers), or a vector aligned to `map`.")
}

#' LD-aware Manhattan plot
#'
#' A per-chromosome (faceted) Manhattan plot for outlier-region analyses. The
#' y-axis is any per-marker value -- a C-score, `-log10(q)`, an `ld_w` window, an
#' association statistic. Points are coloured in one of two modes: **discrete** --
#' by outlier region (or any per-marker group label, e.g. TP/FP), each group its
#' own colour and ungrouped markers grey; or **continuous** -- by a numeric
#' per-marker column via a viridis gradient. Highlight markers (e.g. true QTNs) can
#' be marked with crosses, and a dashed threshold line drawn.
#'
#' Coloring precedence: `group` > `regions` > `colour`. Supply exactly one.
#'
#' @param map data.frame/data.table with `marker`, `Chr`, `Pos`.
#' @param value y-axis value: a named numeric vector (names = markers) or the name
#'   of a numeric column in `map`.
#' @param value_label y-axis label.
#' @param regions Optional list of member-marker vectors (e.g.
#'   `ld_outlier_regions()$regions`); each region gets its own palette colour.
#' @param group Optional named vector (names = markers) of discrete group labels
#'   (markers absent from it are drawn grey). Takes precedence over `regions`; use
#'   it to colour by detection status (`"TP"`/`"FP"`) etc.
#' @param group_colours Optional named colour vector for the `group`/`regions`
#'   labels; defaults to [default_cluster_colours()] recycled.
#' @param colour Optional continuous colouring: a named numeric vector (names =
#'   markers) or a numeric column name in `map`, drawn with a viridis gradient.
#'   Used only when neither `regions` nor `group` is given.
#' @param colour_label Legend title for the continuous `colour`.
#' @param qtn Optional character vector of markers to mark with black crosses.
#' @param hline Optional y value for a dashed horizontal reference line
#'   (e.g. `tau_C`, or `-log10(0.05)`).
#' @param title Plot title.
#' @param pos_unit Divide `Pos` by this for the x-axis (default `1e6` -> Mb).
#' @param x_label x-axis label (default `"Position (Mbp)"`).
#' @param point_size Size of the coloured (grouped / continuous) points.
#'
#' @return A `ggplot` object.
#' @seealso [ld_outlier_regions()], [ld_cscore()], [default_cluster_colours()]
#' @export
ld_manhattan <- function(map, value, value_label = "value",
                         regions = NULL, group = NULL, group_colours = NULL,
                         colour = NULL, colour_label = NULL,
                         qtn = NULL, hline = NULL, title = NULL,
                         pos_unit = 1e6, x_label = "Position (Mbp)",
                         point_size = 1.2) {
  map <- data.table::as.data.table(map)
  markers <- map$marker
  d <- data.frame(marker = markers, Chr = map$Chr, Pos = map$Pos,
                  x = map$Pos / pos_unit,
                  y = .resolve_marker_value(value, markers, map),
                  stringsAsFactors = FALSE)

  mode <- if (!is.null(group)) "group" else if (!is.null(regions)) "regions" else
          if (!is.null(colour)) "continuous" else "none"

  base <- ggplot2::ggplot() +
    ggplot2::facet_wrap(~ Chr, nrow = 1, scales = "free_x") +
    ggplot2::labs(x = x_label, y = value_label, title = title) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(strip.background = ggplot2::element_blank())

  if (mode %in% c("group", "regions")) {
    if (mode == "regions") {
      grp <- rep(NA_character_, nrow(d))
      for (i in seq_along(regions)) grp[markers %in% regions[[i]]] <- as.character(i)
    } else {
      grp <- unname(group[markers])
    }
    d$grp <- ifelse(is.na(grp), "ns", grp)
    labs_ <- setdiff(unique(d$grp), "ns")
    pal <- if (!is.null(group_colours)) group_colours else
           stats::setNames(rep(default_cluster_colours(), length.out = length(labs_)), labs_)
    pal <- c(pal, ns = "grey78")
    show_legend <- !is.null(group)          # regions -> too many to legend; groups -> show
    p <- base +
      ggplot2::geom_point(data = d[d$grp == "ns" & is.finite(d$y), ],
                          ggplot2::aes(x = .data$x, y = .data$y), colour = "grey78", size = 0.6) +
      ggplot2::geom_point(data = d[d$grp != "ns", ],
                          ggplot2::aes(x = .data$x, y = .data$y, colour = .data$grp), size = point_size) +
      ggplot2::scale_colour_manual(values = pal, name = NULL,
                                   guide = if (show_legend) "legend" else "none")
  } else if (mode == "continuous") {
    d$col <- .resolve_marker_value(colour, markers, map)
    p <- base +
      ggplot2::geom_point(data = d[is.finite(d$y), ],
                          ggplot2::aes(x = .data$x, y = .data$y, colour = .data$col), size = point_size) +
      ggplot2::scale_colour_viridis_c(name = colour_label %||% "", na.value = "grey85")
  } else {
    p <- base + ggplot2::geom_point(data = d[is.finite(d$y), ],
                                    ggplot2::aes(x = .data$x, y = .data$y), size = 0.6, colour = "grey40")
  }

  if (!is.null(hline)) p <- p + ggplot2::geom_hline(yintercept = hline, linetype = 2, colour = "grey40")
  if (!is.null(qtn)) {
    dq <- d[d$marker %in% qtn & is.finite(d$y), ]
    if (nrow(dq)) p <- p + ggplot2::geom_point(data = dq, ggplot2::aes(x = .data$x, y = .data$y),
                                               shape = 3, size = 3, stroke = 1, colour = "black")
  }
  p
}

## small null-coalescing helper (local, not exported)
`%||%` <- function(a, b) if (is.null(a)) b else a

utils::globalVariables(c("Chr", "marker"))
