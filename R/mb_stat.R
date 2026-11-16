#' Calculate the MacKenzie-Bailey goodness-of-fit statistic
#'
#' Calculates the MacKenzie-Bailey (MB) goodness-of-fit statistic for a
#' fitted single-season occupancy model. The statistic compares the observed
#' frequencies of detection histories with their expected frequencies under
#' the fitted model.
#'
#' Expected frequencies are calculated separately for each observed
#' missing-data pattern. All possible detection histories consistent with
#' each pattern are enumerated, and their probabilities are calculated
#' using the fitted occupancy and detection probabilities. This allows
#' sites with different patterns of missing visits to contribute to the
#' statistic.
#'
#' @param fit A fitted occupancy model of class \code{unmarkedFitOccu}.
#' @param predictions Optional list of fitted occupancy and detection
#'   probabilities, as returned by \code{.get_occu_predictions()}. Supplying
#'   precomputed predictions avoids recalculating them.
#'
#' @return A list containing:
#' \describe{
#'   \item{statistic}{The overall MacKenzie-Bailey GOF statistic, calculated
#'   as the sum of squared differences between observed and expected
#'   detection-history frequencies, divided by their expected frequencies.}
#'   \item{observed}{A named vector of observed detection-history frequencies.}
#'   \item{expected}{A named vector of expected detection-history frequencies
#'   under the fitted occupancy model.}
#' }
#'
#' @details
#' Detection histories are represented as strings of zeros and ones, with
#' missing visits represented by periods. For each distinct missing-data
#' pattern, all possible combinations of detections and non-detections at
#' observed visits are enumerated.
#'
#' The probability of each detection history is calculated by combining
#' the probability of occupancy and the conditional probability of the
#' observed detection history given occupancy. For unoccupied sites, only
#' histories consisting entirely of non-detections have non-zero probability.
#'
#' Expected frequencies are obtained by summing the probabilities of each
#' history across sites sharing the corresponding missing-data pattern.
#' The overall statistic is then calculated across all detection histories.
#'
#' @keywords internal
mb_stat <- function(fit, predictions = NULL) {

  y <- unmarked::getY(fit)

  n_sites <- nrow(y)
  n_visits <- ncol(y)

  # Model predictions
  if (is.null(predictions)) {
    predictions <- .get_occu_predictions(fit)
  }

  psi_hat <- predictions$psi
  p_mat <- predictions$p_mat

  # Observed histories
  obs_hist <- apply(
    y,
    1,
    function(x) paste(ifelse(is.na(x), ".", x), collapse = "")
  )

  obs_counts <- table(obs_hist)

  # Missing-data patterns
  patterns <- unique(gsub("[01]", "x", obs_hist))

  expected_counts <- numeric()

  for (pattern in patterns) {

    fixed <- strsplit(pattern, "")[[1]]
    x_pos <- which(fixed == "x")

    # Enumerate histories
    if (length(x_pos) == 0) {
      histories <- pattern
    } else {
      combos <- expand.grid(
        rep(list(c("0", "1")), length(x_pos))
      )

      histories <- apply(combos, 1, function(z) {
        h <- fixed
        h[x_pos] <- z
        paste(h, collapse = "")
      })
    }

    site_ids <- which(gsub("[01]", "x", obs_hist) == pattern)

    # Extract probabilities for sites in this pattern
    psi <- psi_hat[site_ids]
    p <- p_mat[site_ids, , drop = FALSE]

    # Only use visits that are observed for this pattern
    observed_pos <- which(fixed != ".")

    p <- p[, observed_pos, drop = FALSE]

    for (hist in histories) {

      h <- as.integer(strsplit(hist, "")[[1]][observed_pos])

      # Probability of this history conditional on occupancy
      prob_occ <- apply(
        ifelse(
          matrix(h, nrow = length(site_ids), ncol = length(h), byrow = TRUE) == 1,
          p,
          1 - p
        ),
        1,
        prod
      )

      # Add probability of being unoccupied
      prob_unocc <- if (all(h == 0)) 1 - psi else rep(0, length(psi))

      expected_counts[hist] <- sum(psi * prob_occ + prob_unocc)
    }
  }

  # Statistic
  obs_all <- table(
    factor(obs_hist, levels = names(expected_counts))
  )

  exp_all <- expected_counts

  T_obs <- sum(
    (obs_all - exp_all)^2 / exp_all,
    na.rm = TRUE
  )

  list(
    statistic = T_obs,
    observed = obs_all,
    expected = exp_all
  )
}
