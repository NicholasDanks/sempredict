#' Observed against predicted, by method
#'
#' Pools the outcome indicators after centring each at its mean and draws
#' observed against predicted for one repetition, with the 45-degree line
#' (one observed point per predicted point) and the fitted slope. Under a
#' construction that over-disperses, such as the PLSc construct-score chain,
#' the fitted line is flatter than the 45-degree line.
#'
#' @param cv An object from [cv_predict()].
#' @param columns Result columns to draw; defaults to all non-benchmark columns.
#' @param rep Repetition to draw.
#' @return A `ggplot` object. Requires the ggplot2 package.
#' @export
plot_dispersion <- function(cv, columns = NULL, rep = 1) {
  stopifnot(inherits(cv, "sem_cv"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("plot_dispersion() needs the ggplot2 package")
  if (is.null(columns)) columns <- cv$columns[cv$column_method != ""]
  Yobs <- as.matrix(cv$data[, cv$ynames, drop = FALSE])
  dd <- do.call(rbind, lapply(columns, function(m) {
    yh <- matrix(cv$yhat[, rep, m, , drop = FALSE], nrow = nrow(Yobs))
    keep <- stats::complete.cases(yh)
    yc <- scale(Yobs[keep, , drop = FALSE], scale = FALSE); hc <- scale(yh[keep, , drop = FALSE], scale = FALSE)
    data.frame(method = m, observed = as.vector(yc), predicted = as.vector(hc),
               slope = sum(yc * hc) / sum(hc^2), stringsAsFactors = FALSE)
  }))
  dd$method <- factor(dd$method, levels = columns)
  lab <- unique(dd[, c("method", "slope")])
  lab$text <- sprintf("slope = %.2f", lab$slope)
  ggplot2::ggplot(dd, ggplot2::aes(predicted, observed)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
    ggplot2::geom_point(alpha = .35, size = 1.2, colour = "grey30") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#C0392B", linewidth = .9) +
    ggplot2::geom_text(data = lab, ggplot2::aes(x = -Inf, y = Inf, label = text),
                       hjust = -0.1, vjust = 1.5, colour = "#C0392B", size = 3.2) +
    ggplot2::facet_wrap(~ method) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "Predicted (difference from the indicator mean)",
                  y = "Observed (difference from the indicator mean)") +
    ggplot2::theme_minimal()
}
