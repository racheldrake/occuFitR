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
mb_stat <- function(fit){

  y <- unmarked::getY(fit)

  n_sites  <- nrow(y)
  n_visits <- ncol(y)

  # -----------------------------
  # Model predictions
  # -----------------------------

  psi_hat <- predict(fit, type = "state")$Predicted

  det_hat <- predict(fit, type = "det")$Predicted

  p_mat <- matrix(
    det_hat,
    nrow = n_sites,
    ncol = n_visits,
    byrow = TRUE
  )

  # -----------------------------
  # Observed histories
  # -----------------------------

  obs_hist <- apply(
    y,
    1,
    function(x) paste(ifelse(is.na(x), ".", x), collapse = "")
  )

  obs_counts <- table(obs_hist)

  # -----------------------------
  # All possible histories for each
  # missing-data pattern
  # -----------------------------

  patterns <- unique(gsub("[01]", "x", obs_hist))

  history_probability <- function(history, psi, p){

    h <- strsplit(history, "")[[1]]

    keep <- h != "."

    h <- h[keep]
    p <- p[keep]

    valid <- !is.na(p)

    h <- h[valid]
    p <- p[valid]

    if(length(p) == 0)
      return(1)

    prob_occ <-
      psi *
      prod(ifelse(h == "1", p, 1 - p))

    prob_unocc <-
      if(all(h == "0")) 1 - psi else 0

    prob_occ + prob_unocc
  }

  expected_counts <- c()

  for(pattern in patterns){

    fixed <- strsplit(pattern, "")[[1]]

    x_pos <- which(fixed == "x")

    if(length(x_pos) == 0){

      histories <- pattern

    } else {

      combos <- expand.grid(rep(list(c("0","1")), length(x_pos)))

      histories <- apply(combos, 1, function(z){

        h <- fixed
        h[x_pos] <- z
        paste(h, collapse = "")

      })

    }

    site_ids <- which(gsub("[01]", "x", obs_hist) == pattern)

    for(hist in histories){

      expected_counts[hist] <-
        sum(
          sapply(site_ids, function(i){

            history_probability(
              hist,
              psi_hat[i],
              p_mat[i,]
            )

          })
        )

    }

  }

  # -----------------------------
  # Statistic
  # -----------------------------

  obs_all <- table(
    factor(
      obs_hist,
      levels = names(expected_counts)
    )
  )

  exp_all <- expected_counts

  T_obs <-
    sum(
      (obs_all - exp_all)^2 / exp_all,
      na.rm = TRUE
    )

  list(
    statistic = T_obs,
    observed = obs_all,
    expected = exp_all
  )

}
