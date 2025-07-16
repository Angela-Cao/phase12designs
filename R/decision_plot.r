#' Decision Map Plot
#'
#' This function creates a decision plot containing customizable decision zones.
#'
#' @param filename File path.
#' @param filetype File type.
#' @param xlab x-axis label. (Default is "Toxicity Probability")
#' @param ylab y-axis label. (Default is "Efficacy Probability")
#' @param x_breaks Numeric vector for x-axis major ticks. (Default is 'c(0, 1')
#' @param y_breaks Numeric vector for y-axis major ticks. (Default is 'c(0, 1')
#' @param x_labels Labels corresponding to `x_breaks`. (Default is 'c(0, 1')
#' @param y_labels Labels corresponding to `y_breaks`. (Default is 'c(0, 1')
#' @param zones A list of rectangular zones to draw, where each rectangle is a list with elements `xmin`, `xmax`, `ymin`, `ymax`, and `color`.
#' @param legend_info A list with two elements: `labels` (character vector) and `colors` (character vector) for the legend.
#' @param title Title of plot. (Default is 'NULL')
#' @param title_pos A numeric vector (x, y) indicating the position of the title text.
#' @param legend_pos A numeric vector (x, y) indicating the position of the legend.
#' @param grid_lines Whether to include background grid lines. (Default is TRUE.)
#' @param plot_size A numeric vector indicating width and height. (Default is c(7, 7)).

#' @export
#' @examples
#' zones <- list(list(xmin = 0.0, xmax = 0.2, ymin = 0, ymax = 1.0, color = "#a8eea8"),
#'               list(xmin = .2, xmax = .3, ymin = 0, ymax = 0.6, color = "#a8eea8"),
#'               list(xmin = .2, xmax = .3, ymin = .6, ymax = 1, color = "#a8d5ee")))
#' plot_decision_zones("test_plot1.png",
#'                     filetype = "png", zones = zones,
#'                     title = "Decision Zones")

decision_plot <- function(filename,
                                 filetype = c("png", "pdf", "svg"), # new parameter
                                 xlab = "Toxicity Probability",
                                 ylab = "Efficacy Probability",
                                 x_breaks = c(0, 1),
                                 y_breaks = c(0, 1),
                                 x_labels = c(0, 1),
                                 y_labels = c(0, 1),
                                 zones = list(), # each entry is entered as list(xmin, xmax, ymin, ymax, color)
                                 legend_info = list(labels = NULL, colors = NULL),
                                 title = NULL,
                                 title_pos = c(0.05, 1.1),
                                 legend_pos = c(0.3, 1.2),
                                 grid_lines = TRUE,
                                 plot_size = c(7, 7)) {

  filetype <- match.arg(filetype)

  switch(filetype,
         png = png(filename, width = plot_size[1], height = plot_size[2], units = "in", res = 300),
         pdf = pdf(filename, width = plot_size[1], height = plot_size[2]),
         svg = svg(filename, width = plot_size[1], height = plot_size[2])
  )

  par(xpd = FALSE)
  plot(c(0, 1), c(0, 1), type = "n",
       xlab = xlab, ylab = ylab,
       xaxt = "n", yaxt = "n",
       cex.lab = 1.4)

  axis(1, at = x_breaks, labels = x_labels)
  axis(2, at = y_breaks, labels = y_labels)


  for (z in zones) {
    rect(z$xmin, z$ymin, z$xmax, z$ymax, col = z$color, border = "transparent")
  }

  if (grid_lines) {
    abline(h = seq(0, 1, 0.2), lty = 2)
    abline(v = seq(0, 1, 0.1), lty = 2)
  }

  par(xpd = TRUE)
  if (!is.null(legend_info$labels) && !is.null(legend_info$colors)) {
    legend(legend_pos[1], legend_pos[2],
           legend = legend_info$labels,
           pch = rep(15, length(legend_info$labels)),
           col = legend_info$colors,
           ncol = length(legend_info$labels))
  }

  if (!is.null(title)) {
    text(title_pos[1], title_pos[2], title, cex = 1.6)
  }

  dev.off()
}
