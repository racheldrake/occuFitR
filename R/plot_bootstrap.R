#' Plot bootstrap distributions of goodness-of-fit statistics
#'
#' Creates a histogram of the parametric bootstrap distribution for a
#' goodness-of-fit statistic, with a vertical dashed line indicating
#' the observed statistic.
#'
#' @param x An object of class \code{occu_gof}, returned by
#'   \code{\link{occu_gof}}.
#' @param test A character string specifying which goodness-of-fit
#'   statistic to plot. Either \code{"mb"}, \code{"psi_p"} or
#'   \code{"p"}.
#'
#' @return A \code{ggplot} object showing the bootstrap distribution
#'   of the selected goodness-of-fit statistic.
#'
#' @details
#' The histogram displays the distribution of statistics obtained
#' from parametric bootstrap simulations. The observed statistic is
#' indicated by a vertical dashed line. The plot title includes the
#' corresponding bootstrap p-value.
#'
#' @seealso \code{\link{occu_gof}}, \code{\link{plot_residuals}},
#'   \code{\link{plot_spatial_residuals}}
#'
#' @export
plot_bootstrap <- function(x, test) {

  test <- match.arg(
    test,
    choices = c("mb", "psi_p", "p")
  )

  ggplot2::ggplot(
    x$bootstrap,
    ggplot2::aes(x = .data[[test]])
  ) +
    ggplot2::geom_histogram(
      bins = 50,
      fill = "grey70",
      colour = "white"
    ) +
    ggplot2::geom_vline(
      xintercept = x$statistic[test],
      linetype = 2,
      linewidth = 1
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = paste(test, "GOF statistic"),
      y = "Frequency",
      title = paste0(
        "Bootstrap GOF test: p = ",
        round(x$p_value[test], 4)
      )
    )
}
