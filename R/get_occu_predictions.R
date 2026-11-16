.get_occu_predictions <- function(fit) {

  y <- unmarked::getY(fit)

  n_sites <- nrow(y)
  n_visits <- ncol(y)

  psi_hat <- predict(fit, type = "state")$Predicted
  det_hat <- predict(fit, type = "det")$Predicted

  p_mat <- matrix(
    det_hat,
    nrow = n_sites,
    ncol = n_visits,
    byrow = TRUE
  )

  list(
    psi = psi_hat,
    det = det_hat,
    p_mat = p_mat
  )
}
