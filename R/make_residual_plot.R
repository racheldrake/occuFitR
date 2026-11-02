#' Plot site-level goodness-of-fit residuals against site covariates
#'
#' Creates diagnostic plots of site-level Pearson residuals against the
#' site covariates in a fitted occupancy model. Separate panels are produced
#' for each covariate, allowing potential patterns in residuals to be
#' identified.
#'
#' Points are coloured according to the number of detections recorded at
#' each site. A smooth trend is added to each panel to help identify
#' systematic relationships between residuals and covariate values.
#'
#' @param fit A fitted occupancy model of class \code{unmarkedFitOccu}.
#' @param y A matrix of observed detection histories, with sites in rows
#'   and visits in columns. Missing observations should be represented by
#'   \code{NA}.
#' @param residuals A numeric vector of site-level Pearson residuals.
#' @param residual_sites An integer vector identifying the sites
#'   corresponding to the supplied residuals. Site indices should correspond
#'   to the row indices in \code{y} and the site covariates in \code{fit}.
#' @param covariates A character vector giving the names of the site
#'   covariates to include in the plots.
#' @param test_name A character string identifying the GOF test associated
#'   with the residuals. Used to label the y-axis.
#'
#' @return A \code{ggplot} object showing site-level residuals against each
#'   specified covariate, faceted by covariate. Returns \code{NULL} if
#'   \code{covariates} has length zero.
#'
#' @details
#' Site-level residuals are joined to the corresponding site covariates
#' and detection counts using site indices. The data are then reshaped
#' into long format to allow each covariate to be displayed in a separate
#' facet.
#'
#' A generalised additive model smooth is fitted to the residuals in each
#' panel, with a confidence interval displayed around the fitted trend.
#' The x-axis scales are allowed to vary between covariates.
#'
#' @keywords internal
make_residual_plot <- function(
    fit,
    y,
    residuals,
    residual_sites,
    covariates,
    test_name
) {

  # No site-level covariates to plot
  if (length(covariates) == 0) {
    return(NULL)
  }

  site_covs <- fit@data@siteCovs |>
    as_tibble()

  site_covs <- site_covs |>
    mutate(
      site = seq_len(nrow(site_covs))
    )

  detection_counts <- tibble(
    site = seq_len(nrow(y)),
    n_detections = rowSums(y, na.rm = TRUE)
  ) |>
    mutate(
      n_detections = factor(
        n_detections,
        levels = sort(unique(n_detections)),
        ordered = TRUE
      )
    )

  residual_df <- tibble(
    site = residual_sites,
    residual = residuals
  ) |>
    left_join(
      detection_counts,
      by = "site"
    ) |>
    left_join(
      site_covs,
      by = "site"
    )

  residual_long <- residual_df |>
    select(
      site,
      residual,
      n_detections,
      all_of(covariates)
    ) |>
    pivot_longer(
      cols = all_of(covariates),
      names_to = "covariate",
      values_to = "value"
    )

  ggplot(
    residual_long,
    aes(
      x = value,
      y = residual
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = 2,
      colour = "grey50"
    ) +
    geom_point(
      aes(colour = n_detections),
      alpha = 0.6
    ) +
    geom_smooth(
      se = TRUE,
      method = "gam",
      formula = y ~ s(x, bs = "cs", k = 6),
      colour = "black"
    ) +
    facet_wrap(
      ~ covariate,
      scales = "free_x"
    ) +
    theme_bw() +
    labs(
      x = "Covariate value",
      y = paste0(
        "Site-level Pearson ",
        test_name,
        " residual"
      )
    )  +
    scale_colour_viridis_d(
      direction = -1,
      option = "G",
      begin = 0,
      end = 0.9
    )
}
