#' Plot goodness-of-fit diagnostics
#'
#' Generic plotting method for goodness-of-fit results from
#' occupancy models.
#'
#' @param x An object of class \code{occu_gof}.
#' @param test A character string specifying which goodness-of-fit
#'   test to plot. Either \code{"psi_p"} or \code{"p"}.
#' @param type A character string specifying the type of diagnostic
#'   plot. Either \code{"bootstrap"}, \code{"residuals"} or \code{"spatial"}.
#' @param ... Additional arguments passed to the corresponding
#'   plotting function.
#'
#' @return The plot or plots returned by the corresponding plotting
#'   function.
#'
#' @seealso \code{\link{occu_gof}},
#'   \code{\link{plot_residuals}},
#'   \code{\link{plot_spatial_residuals}}
#'
#' @export
plot.occu_gof <- function(
    x,
    test,
    type = "bootstrap",
    ...
) {

  type <- match.arg(type)

  if (type == "bootstrap") {
    plot_bootstrap(
      x,
      test = test,
      ...
    )
  }


  if (type == "residuals") {
    plot_residuals(
      x,
      test = test,
      ...
    )
  }

  if (type == "spatial") {
    plot_spatial_residuals(
      x,
      test = test,
      ...
    )
  }
}
