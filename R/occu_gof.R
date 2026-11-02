#' Parametric bootstrap goodness-of-fit tests for occupancy models
#'
#' Performs parametric bootstrap goodness-of-fit (GOF) tests for a fitted
#' single-season occupancy model. Simulated detection histories are generated
#' using the fitted occupancy and detection probabilities, while retaining
#' the observed pattern of missing visits. The model is refitted to each
#' simulated dataset and the GOF statistics recalculated.
#'
#' Three aggregate GOF statistics are calculated: the MacKenzie-Bailey (MB)
#' statistic, a site-level occupancy and detection probability (\eqn{\psi-p})
#' statistic, and a detection probability (\eqn{p}) statistic. Bootstrap
#' p-values are calculated as the proportion of simulated statistics greater
#' than or equal to the corresponding observed statistic.
#'
#' Optionally, the function also produces histograms of the bootstrap
#' distributions, with the observed statistics indicated, and plots of
#' site-level residuals against the site covariates.
#'
#' @param fit A fitted occupancy model of class \code{unmarkedFitOccu}.
#' @param n_boot Number of parametric bootstrap replicates. Defaults to
#'   \code{1000}.
#' @param seed Optional random seed for reproducibility. If \code{NULL},
#'   no seed is set.
#' @param plot Logical. If \code{TRUE}, produce plots of the bootstrap
#'   distributions and site-level residuals. Defaults to \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{statistic}{A named vector of the observed MB, \eqn{\psi-p}, and
#'   \eqn{p} GOF statistics.}
#'   \item{p_value}{A named vector of bootstrap p-values for the three
#'   GOF statistics.}
#'   \item{bootstrap}{A data frame containing the simulated values of
#'   each GOF statistic across bootstrap replicates.}
#'   \item{true_residuals}{A list containing the observed site-level
#'   residuals for the \eqn{\psi-p} and \eqn{p} statistics.}
#'   \item{plot}{A list of plots. This includes histograms of the three
#'   bootstrap distributions (\code{mb}, \code{psi_p}, and \code{p}) and
#'   site-level residual plots for the \eqn{\psi-p} and \eqn{p} statistics
#'   (\code{site_psi_p} and \code{site_p}). If \code{plot = FALSE},
#'   all elements of this list are \code{NULL}.}
#' }
#'
#' @details
#' For each bootstrap replicate, occupancy states are simulated using the
#' fitted site-level occupancy probabilities. Detection histories are then
#' simulated using the fitted detection probabilities, conditional on
#' occupancy. The observed pattern of missing visits is retained.
#'
#' The occupancy model is refitted to each simulated dataset before the GOF
#' statistics are recalculated. The resulting bootstrap distributions are
#' compared with the observed statistics to obtain upper-tail p-values.
#'
#' The site-level residual plots compare residuals with the site covariates
#' supplied to the original occupancy model. Residuals for the detection
#' probability statistic are calculated only for sites where the species
#' was detected.
#'
#' @seealso \code{\link[unmarked:occu]{unmarked::occu}},
#'   \code{\link[unmarked:getY]{unmarked::getY}}
#'
#' @examples
#' \dontrun{
#'   fit <- unmarked::occu(~ 1 ~ 1, data = umf)
#'
#'   gof <- occu_gof(
#'     fit,
#'     n_boot = 1000,
#'     seed = 123
#'   )
#'
#'   gof$statistic
#'   gof$p_value
#'   gof$plot$mb
#' }
#'
#' @export
occu_gof <- function(fit, n_boot = 1000, seed = NULL, plot = TRUE) {

  if (!is.null(seed)) set.seed(seed)

  y <- unmarked::getY(fit)

  n_sites <- nrow(y)
  n_visits <- ncol(y)

  # --- model predictions ---

  psi_hat <- predict(fit, type = "state")$Predicted

  det_hat <- predict(fit, type = "det")$Predicted

  p_mat <- matrix(
    det_hat,
    nrow = n_sites,
    ncol = n_visits,
    byrow = TRUE
  )


  # --- observed data stats ---

  mb_obs <- mb_stat(fit)
  psi_p_obs <- psi_p_stat(fit)
  p_obs <- p_stat(fit)

  T_obs <- c(
    mb = mb_obs$statistic,
    psi_p = psi_p_obs$statistic,
    p = p_obs$statistic
  )

  # --- save true site-level residuals ---

  true_residuals <- list(
    psi_p = psi_p_obs$residuals,
    p = p_obs$residuals
  )


  # --- site-level covariates, whether included in model or no ---

  plot_covs <- names(fit@data@siteCovs)



  # --- simulation function ---

  sim_once <- function() {

    y_sim <- matrix(NA, n_sites, n_visits)

    for (i in seq_len(n_sites)) {

      # Generate occupancy state from site psi
      z <- rbinom(1, 1, psi_hat[i])

      # Decide if we generate detections
      if (z == 1) {

        # generate detections for each visit
        # where we have an observation
        for (t in seq_len(n_visits)) {

          if (!is.na(y[i, t])) {

            y_sim[i, t] <- rbinom(
              1,
              1,
              p_mat[i, t]
            )

          }

        }

      } else {

        for (t in seq_len(n_visits)) {

          # also need to include missing-ness here
          if (!is.na(y[i, t])) {

            y_sim[i, t] <- 0

          }

        }
      }
    }


    # fit model to generated data

    umf_sim <- unmarkedFrameOccu(
      y = y_sim,
      siteCovs = siteCovs(fit@data),
      obsCovs = obsCovs(fit@data)
    )

    fit_sim <- update(
      fit,
      data = umf_sim
    )


    # calculate simulated data stats

    mb_sim <- mb_stat(fit_sim)
    psi_p_sim <- psi_p_stat(fit_sim)
    p_sim <- p_stat(fit_sim)

    # arrange outputs
    list(

      # aggregate GOF statistics
      statistic = c(
        mb = mb_sim$statistic,
        psi_p = psi_p_sim$statistic,
        p = p_sim$statistic
      )
    )

  }


  # --- bootstrap ---

  boot_results <- future_replicate(
    n_boot,
    sim_once(),
    simplify = FALSE
  )

  # collect aggregate statistics

  boot_stats <- do.call(
    rbind,
    lapply(
      boot_results,
      function(x) x$statistic
    )
  )

  boot_stats <- as.data.frame(boot_stats)

  # --- p-values ---

  p_values <- c(
    mb = mean(
      boot_stats$mb >= T_obs["mb"],
      na.rm = TRUE
    ),

    psi_p = mean(
      boot_stats$psi_p >= T_obs["psi_p"],
      na.rm = TRUE
    ),

    p = mean(
      boot_stats$p >= T_obs["p"],
      na.rm = TRUE
    )
  )


  # --- plots ---

  if (plot) {

    # --- MB plot ---

    plot_mb <- ggplot(
      boot_stats,
      aes(x = mb)
    ) +
      geom_histogram(
        bins = 50,
        fill = "grey70",
        colour = "white"
      ) +
      geom_vline(
        xintercept = T_obs["mb"],
        linetype = 2,
        linewidth = 1
      ) +
      theme_bw() +
      labs(
        x = "MB GOF statistic",
        y = "Frequency",
        title = paste0(
          "MB GOF test: p = ",
          round(p_values["mb"], 4)
        )
      )


    # --- psi-p plot ---

    plot_psi_p <- ggplot(
      boot_stats,
      aes(x = psi_p)
    ) +
      geom_histogram(
        bins = 50,
        fill = "grey70",
        colour = "white"
      ) +
      geom_vline(
        xintercept = T_obs["psi_p"],
        linetype = 2,
        linewidth = 1
      ) +
      theme_bw() +
      labs(
        x = "Psi-p GOF statistic",
        y = "Frequency",
        title = paste0(
          "Psi-p GOF test: p = ",
          round(p_values["psi_p"], 4)
        )
      )


    # --- p plot ---

    plot_p <- ggplot(
      boot_stats,
      aes(x = p)
    ) +
      geom_histogram(
        bins = 50,
        fill = "grey70",
        colour = "white"
      ) +
      geom_vline(
        xintercept = T_obs["p"],
        linetype = 2,
        linewidth = 1
      ) +
      theme_bw() +
      labs(
        x = "P GOF statistic",
        y = "Frequency",
        title = paste0(
          "P GOF test: p = ",
          round(p_values["p"], 4)
        )
      )

    # --- site-level residual plots ---

    plot_site_psi_p <- make_residual_plot(
      fit = fit,
      y = y,
      residuals = psi_p_obs$residuals,
      residual_sites = seq_along(psi_p_obs$residuals),
      covariates = plot_covs,
      test_name = "psi-p"
    )

    # p residuals correspond only to sites where the
    # species was detected
    p_sites <- which(
      rowSums(y, na.rm = TRUE) > 0
    )

    plot_site_p <- make_residual_plot(
      fit = fit,
      y = y,
      residuals = p_obs$residuals,
      residual_sites = p_sites,
      covariates = plot_covs,
      test_name = "p"
    )


  } else {

    plot_mb <- NULL
    plot_psi_p <- NULL
    plot_p <- NULL
    plot_site_psi_p <- NULL
    plot_site_p <- NULL

  }


  list(
    statistic = T_obs,
    p_value = p_values,
    bootstrap = boot_stats,
    true_residuals = true_residuals,

    plot = list(
      mb = plot_mb,
      psi_p = plot_psi_p,
      p = plot_p,
      site_psi_p = plot_site_psi_p,
      site_p = plot_site_p
    )
  )

}
