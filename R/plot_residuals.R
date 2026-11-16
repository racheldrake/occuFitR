#' Plot site-level goodness-of-fit residuals against site covariates
#'
#' Creates diagnostic plots of site-level Pearson residuals against the
#' site covariates stored in an \code{occu_gof} object. Separate panels
#' are produced for each covariate, allowing potential patterns in
#' residuals to be identified.
#'
#' Points are coloured according to the number of detections recorded at
#' each site. Numeric covariates are plotted with a generalised additive
#' model smooth and confidence interval, while categorical covariates
#' are plotted using boxplots and jittered points.
#'
#' @param x An object of class \code{occu_gof}, returned by
#'   \code{\link{occu_gof}}.
#' @param test A character string specifying which residuals to plot.
#'   Either \code{"psi_p"} or \code{"p"}.
#' @param covariates An optional character vector giving the names of the
#'   site covariates to include in the plots. Defaults to all available
#'   site covariates.
#'
#' @return If both numeric and categorical covariates are included,
#'   returns a named list containing two \code{ggplot} objects:
#'   \code{numeric} and \code{categorical}. If only one type of
#'   covariate is included, returns a single \code{ggplot} object.
#'   Returns \code{NULL} if no suitable covariates are available.
#'
#' @details
#' Site-level residuals are joined to the corresponding site covariates
#' and detection counts using site indices. Numeric and categorical
#' covariates are plotted separately to accommodate their different
#' data types.
#'
#' For numeric covariates, a generalised additive model smooth is fitted
#' to the residuals in each panel, with a confidence interval displayed
#' around the fitted trend. For categorical covariates, boxplots show
#' the distribution of residuals within each category, with individual
#' observations displayed as jittered points.
#'
#' Residuals for the \eqn{\psi-p} statistic are available for all sites,
#' whereas residuals for the \eqn{p} statistic are available only for
#' sites where the species was detected.
#'
#' @seealso \code{\link{occu_gof}}
#'
#' @export
plot_residuals <- function(
    x,
    test,
    covariates = NULL
) {

  # Check test selection
  if (!test %in% c("psi_p", "p")) {
    stop(
      "'test' must be either 'psi_p' or 'p'.",
      call. = FALSE
    )
  }

  # Extract residuals and site covariates
  residual_df <- x$residuals[[test]] |>
    tibble::as_tibble() |>
    dplyr::left_join(
      x$site_covs,
      by = "site"
    )

  # Identify available site covariates
  available_covariates <- setdiff(
    names(x$site_covs),
    "site"
  )

  # Use all available covariates by default
  if (is.null(covariates)) {
    covariates <- available_covariates
  }

  # Check that requested covariates are available
  if (!all(covariates %in% available_covariates)) {
    stop(
      "Some specified covariates are not available in the GOF object.",
      call. = FALSE
    )
  }

  # No site-level covariates to plot
  if (length(covariates) == 0) {
    return(NULL)
  }

  # Calculate detection counts
  y <- unmarked::getY(x$fit)

  detection_counts <- tibble::tibble(
    site = seq_len(nrow(y)),
    n_detections = rowSums(y, na.rm = TRUE)
  ) |>
    dplyr::mutate(
      n_detections = factor(
        n_detections,
        levels = sort(unique(n_detections)),
        ordered = TRUE
      )
    )

  # Join detection counts
  residual_df <- residual_df |>
    dplyr::left_join(
      detection_counts,
      by = "site"
    )

  # Identify numeric and categorical covariates
  numeric_covariates <- covariates[
    vapply(
      x$site_covs[covariates],
      is.numeric,
      logical(1)
    )
  ]

  categorical_covariates <- setdiff(
    covariates,
    numeric_covariates
  )

  # Reshape numeric covariates
  if (length(numeric_covariates) > 0) {
    numeric_data <- residual_df |>
      dplyr::select(
        site, residual, n_detections,
        dplyr::all_of(numeric_covariates)
      ) |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(numeric_covariates),
        names_to = "covariate",
        values_to = "value"
      )

    numeric_plot <- ggplot2::ggplot(
      numeric_data,
      ggplot2::aes(
        x = value,
        y = residual,
        colour = n_detections
      )
    ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed"

























      ) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::geom_smooth(
        method = "gam",
        formula = y ~ s(x, bs = "cs"),
        se = TRUE,
        colour = "black"
      ) +
      ggplot2::facet_wrap(
        ~covariate,
        scales = "free_x"
      ) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x = "Covariate",
        y = "Residual",
        colour = "Number of detections"
      ) +
      ggplot2::scale_colour_viridis_d()
  }

  # Reshape categorical covariates
  if (length(categorical_covariates) > 0) {
    categorical_data <- residual_df |>
      dplyr::select(
        site, residual, n_detections,
        dplyr::all_of(categorical_covariates)
      ) |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(categorical_covariates),
        names_to = "covariate",
        values_to = "value"
      ) |>
      dplyr::mutate(
        value = as.factor(value)
      )

    categorical_plot <- ggplot2::ggplot(
      categorical_data,
      ggplot2::aes(
        x = value,
        y = residual,
        colour = n_detections
      )
    ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed"
      ) +
      ggplot2::geom_boxplot(
        ggplot2::aes(group = value),
        colour = "black",
        alpha = 0.2,
        outlier.shape = NA
      ) +
      ggplot2::geom_jitter(
        width = 0.15,
        alpha = 0.6
      ) +
      ggplot2::facet_wrap(
        ~covariate,
        scales = "free_x"
      ) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x = "Covariate",
        y = "Residual",
        colour = "Number of detections"
      ) +
      ggplot2::scale_colour_viridis_d()
  }

  # Return whichever plots are applicable
  if (length(numeric_covariates) > 0 &&
      length(categorical_covariates) > 0) {

    return(list(
      numeric = numeric_plot,
      categorical = categorical_plot
    ))

  } else if (length(numeric_covariates) > 0) {

    return(numeric_plot)

  } else {

    return(categorical_plot)

  }
}
