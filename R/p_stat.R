#' Calculate the detection probability goodness-of-fit statistic
#'
#' Calculates a site-level goodness-of-fit (GOF) statistic for the detection
#' component of a fitted single-season occupancy model. The statistic compares
#' the observed number of detections at each site with the expected number,
#' conditional on the species having been detected at that site.
#'
#' Only sites with at least one observed detection are included in the
#' calculation. For each included site, the expected number of detections
#' and its variance are calculated from the fitted detection probabilities.
#' Standardised site-level residuals are then calculated, and the overall
#' GOF statistic is obtained by summing their squared values.
#'
#' @param fit A fitted occupancy model of class \code{unmarkedFitOccu}.
#'
#' @return A list containing:
#' \describe{
#'   \item{statistic}{The overall detection probability GOF statistic,
#'   calculated as the sum of squared standardised site-level residuals.}
#'   \item{residuals}{A numeric vector of standardised site-level residuals
#'   for sites where the species was detected.}
#'   \item{observed}{The observed number of detections at each included site.}
#'   \item{expected}{The expected number of detections at each included site,
#'   conditional on at least one detection.}
#'   \item{variance}{The variance of the number of detections at each included
#'   site, conditional on at least one detection.}
#' }
#'
#' @details
#' The expected number of detections and corresponding variance are
#' calculated conditional on at least one detection occurring at a site.
#' This accounts for the fact that sites with no detections are excluded
#' from the detection component of the test.
#'
#' The overall statistic is calculated by summing the squared standardised
#' differences between observed and expected detections across included
#' sites. Larger values indicate greater discrepancies between observed
#' and expected detection frequencies.
#'
#' @keywords internal
p_stat <- function(fit){
  y <- unmarked::getY(fit)

  n_sites <- nrow(y)
  n_visits <- ncol(y)

  psi_hat <- predict(fit, type="state")$Predicted
  det_hat <- predict(fit, type="det")$Predicted

  p_mat <- matrix(det_hat,
                  nrow=n_sites,
                  ncol=n_visits,
                  byrow=TRUE)

  p_dot <- apply(p_mat,1,function(p){
    p <- p[!is.na(p)]
    1-prod(1-p)
  })

  detected_sites <- rowSums(y, na.rm = TRUE) > 0

  # --- keep only detected sites ---
  y <- y[detected_sites, , drop = FALSE]
  p_mat <- p_mat[detected_sites, , drop = FALSE]

  # update number of sites
  n_sites <- nrow(y)

  # compute expectations after filtering to detected only
  p_sum <- rowSums(p_mat, na.rm = TRUE)
  p_sq <- rowSums(p_mat^2, na.rm = TRUE)

  p_dot <- 1 - apply(1 - p_mat, 1, prod, na.rm = TRUE)

  E_i <- p_sum / p_dot

  Var_i <- (p_sum + p_sum^2 - p_sq)/(p_dot) -
    (p_sum/p_dot)^2

  O_i <- rowSums(y, na.rm = TRUE)

  # --- test statistic ---
  T_i <- (O_i-E_i)/sqrt(Var_i)

  T_total <- sum((O_i-E_i)^2/Var_i,na.rm=TRUE)

  list(
    statistic=T_total,
    residuals = T_i,
    observed=O_i,
    expected=E_i,
    variance=Var_i
  )

}

