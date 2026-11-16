#' Calculate the occupancy and detection probability goodness-of-fit statistic
#'
#' Calculates a site-level goodness-of-fit (GOF) statistic for the occupancy
#' and detection components of a fitted single-season occupancy model. The
#' statistic compares the observed number of detections at each site with
#' the expected number, accounting for both the probability of occupancy
#' and the probability of detection.
#'
#' For each site, the expected number of detections and its variance are
#' calculated using the fitted occupancy and detection probabilities.
#' Standardised site-level residuals are then calculated, and the overall
#' GOF statistic is obtained by summing their squared values.
#'
#' @param fit A fitted occupancy model of class \code{unmarkedFitOccu}.
#' @param predictions Optional list of fitted occupancy and detection
#'   probabilities, as returned by \code{.get_occu_predictions()}. Supplying
#'   precomputed predictions avoids recalculating them.
#'
#' @return A list containing:
#' \describe{
#'   \item{statistic}{The overall occupancy and detection probability GOF
#'   statistic, calculated as the sum of squared standardised site-level
#'   residuals.}
#'   \item{residuals}{A numeric vector of standardised site-level residuals.}
#'   \item{observed}{The observed number of detections at each site.}
#'   \item{expected}{The expected number of detections at each site,
#'   calculated using the fitted occupancy and detection probabilities.}
#'   \item{variance}{The variance of the number of detections at each site,
#'   calculated using the fitted occupancy and detection probabilities.}
#' }
#'
#' @details
#' The expected number of detections at site \(i\) is calculated as the
#' product of the fitted occupancy probability and the sum of the fitted
#' detection probabilities across visits.
#'
#' The variance accounts for variation in both occupancy state and
#' detection outcomes, including the dependence between detections at
#' different visits induced by the shared occupancy state.
#'
#' The overall statistic is calculated by summing the squared standardised
#' differences between observed and expected detections across sites.
#' Larger values indicate greater discrepancies between observed and
#' expected detection frequencies.
#'
#' @keywords internal
psi_p_stat <- function(fit, predictions = NULL){

  y <- unmarked::getY(fit)

  n_sites <- nrow(y)
  n_visits <- ncol(y)

  if (is.null(predictions)) {
    predictions <- .get_occu_predictions(fit)
  }

  psi_hat <- predictions$psi
  p_mat <- predictions$p_mat

  # Expected number of detections
  E_i <- psi_hat * rowSums(p_mat, na.rm = TRUE)

  # First component of variance
  term1 <- rowSums(
    psi_hat * p_mat * (1 - psi_hat * p_mat),
    na.rm = TRUE
  )

  # Second component of variance
  p_sum <- rowSums(p_mat, na.rm = TRUE)
  p_sq_sum <- rowSums(p_mat^2, na.rm = TRUE)

  pair_sum <- (p_sum^2 - p_sq_sum) / 2

  term2 <- 2 * psi_hat * (1 - psi_hat) * pair_sum

  Var_i <- term1 + term2

  # Observed detections
  O_i <- rowSums(y, na.rm = TRUE)

  # Residuals and overall statistic
  T_i <- (O_i - E_i) / sqrt(Var_i)

  T_total <- sum((O_i - E_i)^2 / Var_i, na.rm = TRUE)

  list(
    statistic = T_total,
    residuals = T_i,
    observed = O_i,
    expected = E_i,
    variance = Var_i
  )
}
