#' Summarise occupancy model GOF results
#'
#' @param object An object of class \code{occu_gof}.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{summary.occu_gof}.
#'
#' @seealso
#' \code{\link[unmarked:occu]{unmarked::occu}},
#' \code{\link[unmarked:getY]{unmarked::getY}},
#' \code{\link{occu_gof}},
#' \code{\link{plot_bootstrap}},
#' \code{\link{plot_residuals}},
#' \code{\link{plot_spatial_residuals}}
#'
#' @export
summary.occu_gof <- function(object, ...) {

  fit <- object$fit

  y <- unmarked::getY(fit)

  results <- data.frame(
    test = c(
      "MacKenzie-Bailey",
      "Psi-p",
      "Detection probability"
    ),
    statistic = unname(object$statistic),
    p_value = unname(object$p_value)
  )

  out <- list(
    results = results,
    n_sites = nrow(y),
    n_visits = ncol(y),
    n_boot = object$n_boot,
    fit = fit
  )

  class(out) <- "summary.occu_gof"

  out
}
