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
psi_p_stat <- function(fit){

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

  E_i <- psi_hat * rowSums(p_mat,na.rm=TRUE)


  term1 <- rowSums(
    psi_hat*p_mat*(1-psi_hat*p_mat),
    na.rm=TRUE
  )

  term2 <- sapply(seq_len(n_sites), function(i){

    p <- p_mat[i,!is.na(p_mat[i,])]

    if(length(p)<2) return(0)

    pairs <- combn(p,2)

    2*psi_hat[i]*(1-psi_hat[i])*
      sum(pairs[1,]*pairs[2,])

  })

  Var_i <- term1 + term2

  O_i <- rowSums(y,na.rm=TRUE)

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
