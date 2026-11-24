
#' Plot spatial patterns in site-level goodness-of-fit residuals
#'
#' Creates two diagnostic plots showing the spatial distribution of
#' site-level Pearson residuals from an \code{occu_gof} object.
#'
#' The first plot displays a two-dimensional generalised additive model
#' smooth of residuals across the convex hull of the study sites, using
#' a divergent colour scale centred on zero. The second displays
#' individual sites as circles, with colour indicating the sign of
#' the residual and size representing its absolute magnitude.
#'
#' @param x An object of class \code{occu_gof}, returned by
#'   \code{\link{occu_gof}}.
#' @param coords A character vector of length two giving the names of the
#'   coordinate columns in the site covariates. The first element
#'   specifies the x-coordinate and the second the y-coordinate.
#' @param test A character string specifying which residuals to plot.
#'   Either \code{"psi_p"} or \code{"p"}.
#' @param crs The coordinate reference system of the supplied coordinates.
#'   Can be a CRS object, EPSG code or PROJ string accepted by
#'   \code{sf::st_crs}.
#' @param boundary Either \code{"convex_hull"} to use the convex hull
#'   of the sampling sites, or an \code{sf} or \code{sfc} polygon object
#'   defining the area over which predictions should be made. If a
#'   supplied boundary has a different CRS, it is transformed to the
#'   CRS specified by \code{crs}.
#' @param ... additional arguments to be added to customise the \code{mgcv::gam} smooth
#'   only. To update plots please use the function output with the usual '+'
#'   syntax for ggplot2. See \code{\link{mgcv::gam}}
#'
#'
#' @return A named list containing two \code{ggplot} objects:
#'   \code{smooth}, showing the spatial smooth of residuals, and
#'   \code{sites}, showing the residuals at individual sites.
#'
#' @details
#' Residuals are joined to their corresponding site coordinates using
#' site indices. The spatial smooth is estimated using a two-dimensional
#' thin-plate regression spline fitted with \code{\link[mgcv]{gam}}. Predictions are
#' restricted to either the convex hull of the sampling
#' sites or a user-supplied polygon boundary.
#'
#' The site-level plot uses a divergent colour scheme to distinguish
#' positive and negative residuals. Circle sizes are scaled according
#' to the absolute residual values.
#'
#' Residuals for the \eqn{\psi-p} statistic are available for all sites,
#' whereas residuals for the \eqn{p} statistic are available only for
#' sites where the species was detected.
#'
#' @seealso \code{\link{occu_gof}}, \code{\link{plot_residuals}}
#'
#' @export
plot_spatial_residuals <- function(
    x,
    coords,
    test,
    crs,
    boundary = "convex_hull",
    ...
) {

  if (!test %in% c("psi_p", "p")) {
    stop(
      "'test' must be either 'psi_p' or 'p'.",
      call. = FALSE
    )
  }

  # Check coordinate specification
  if (length(coords) != 2) {
    stop(
      "'coords' must contain exactly two coordinate column names.",
      call. = FALSE
    )
  }

  available_covariates <- names(x$site_covs)

  # Check CRS
  crs <- sf::st_crs(crs)

  if (is.na(crs)) {
    stop(
      "Please supply a valid coordinate reference system.",
      call. = FALSE
    )
  }

  if (!all(coords %in% available_covariates)) {
    stop(
      "One or both coordinate columns are not available in the GOF object.",
      call. = FALSE
    )
  }

  if (!all(vapply(x$site_covs[coords], is.numeric, logical(1)))) {
    stop(
      "Both coordinate columns must be numeric.",
      call. = FALSE
    )
  }

  # Extract residuals and coordinates
  plot_data <- x$residuals[[test]] |>
    tibble::as_tibble() |>
    dplyr::left_join(
      x$site_covs,
      by = "site"
    ) |>
    dplyr::mutate(
      x_coord = .data[[coords[1]]],
      y_coord = .data[[coords[2]]]
    ) |>
    dplyr::filter(
      is.finite(x_coord),
      is.finite(y_coord),
      is.finite(residual)
    )

  if (nrow(plot_data) < 5) {
    stop(
      "At least five sites with valid coordinates and residuals are required.",
      call. = FALSE
    )
  }

  # Check for sufficient spatial variation
  n_unique_coords <- nrow(
    dplyr::distinct(plot_data, x_coord, y_coord)
  )

  if (n_unique_coords < 5) {
    stop(
      "At least five unique coordinate locations are required.",
      call. = FALSE
    )
  }

  if (length(unique(plot_data$x_coord)) < 2 ||
      length(unique(plot_data$y_coord)) < 2) {
    stop(
      "Coordinates must contain variation in both dimensions.",
      call. = FALSE
    )
  }

  # Identify the sign of each residual
  plot_data <- plot_data |>
    dplyr::mutate(
      residual_sign = factor(
        dplyr::case_when(
          residual > 0 ~ "Positive",
          residual < 0 ~ "Negative",
          TRUE ~ "Zero"
        ),
        levels = c("Negative", "Zero", "Positive")
      )
    )

  # Convert coordinates to sf
  plot_sf <- sf::st_as_sf(
    plot_data,
    coords = c("x_coord", "y_coord"),
    crs = crs,
    remove = FALSE
  )

  # --------------------------------------------------
  # 1. Spatial smooth of residuals
  # --------------------------------------------------

  # Fit a two-dimensional GAM
  spatial_model <- mgcv::gam(
    residual ~ s(
      x_coord,
      y_coord,
      bs = "tp",
      k = min(50, n_unique_coords - 1)
    ),
    data = plot_data,
    method = "REML",
    ...
  )

  # Extract spatial smooth significance
  smooth_summary <- summary(spatial_model)$s.table

  smooth_edf <- smooth_summary[1, "edf"]
  smooth_p <- smooth_summary[1, "p-value"]

  smooth_label <- if (smooth_p < 0.05) {
    "Spatial smooth: p < 0.05"
  } else {
    paste0(
      "Spatial smooth: p = ",
      signif(smooth_p, 3)
    )
  }

  # Generate prediction grid
  prediction_grid <- expand.grid(
    x_coord = seq(
      min(plot_data$x_coord),
      max(plot_data$x_coord),
      length.out = 100
    ),
    y_coord = seq(
      min(plot_data$y_coord),
      max(plot_data$y_coord),
      length.out = 100
    )
  )

  # Predict residuals across the grid
  prediction_grid$residual <- predict(
    spatial_model,
    newdata = prediction_grid
  )

  # --------------------------------------------------
  # Define prediction boundary
  # --------------------------------------------------

  if (identical(boundary, "convex_hull")) {

    # Use convex hull of sampling sites
    prediction_boundary <- sf::st_convex_hull(
      sf::st_union(plot_sf)
    )

  } else {

    # Check supplied boundary
    if (!inherits(boundary, "sf") &&
        !inherits(boundary, "sfc")) {
      stop(
        "'boundary' must be 'convex_hull', an sf object, or an sfc object.",
        call. = FALSE
      )
    }

    # Check boundary CRS
    boundary_crs <- sf::st_crs(boundary)

    if (is.na(boundary_crs)) {
      stop(
        "The supplied boundary must have a defined CRS.",
        call. = FALSE
      )
    }

    # Extract geometry
    boundary_geometry <- if (inherits(boundary, "sf")) {
      sf::st_geometry(boundary)
    } else {
      boundary
    }

    # Check geometry type
    geometry_types <- as.character(
      sf::st_geometry_type(boundary_geometry)
    )

    if (!all(geometry_types %in% c("POLYGON", "MULTIPOLYGON"))) {
      stop(
        "The supplied boundary must contain polygon or multipolygon geometries.",
        call. = FALSE
      )
    }

    # Transform to coordinate CRS if necessary
    boundary_geometry <- sf::st_transform(
      boundary_geometry,
      crs
    )

    # Combine multiple polygons into one geometry
    prediction_boundary <- sf::st_union(
      boundary_geometry
    )
  }

  # Check that the boundary overlaps the sampling sites
  if (!any(lengths(
    sf::st_intersects(plot_sf, prediction_boundary)
  ) > 0)) {
    stop(
      "The prediction boundary does not overlap any sampling sites.",
      call. = FALSE
    )
  }

  # Convert prediction grid to sf
  grid_sf <- sf::st_as_sf(
    prediction_grid,
    coords = c("x_coord", "y_coord"),
    crs = crs
  )

  # Retain predictions within the selected boundary
  inside_boundary <- lengths(
    sf::st_intersects(grid_sf, prediction_boundary)
  ) > 0

  prediction_grid <- prediction_grid[inside_boundary, ]

  if (nrow(prediction_grid) == 0) {
    stop(
      "No prediction grid cells fall within the supplied boundary.",
      call. = FALSE
    )
  }

  smooth_plot <- ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = prediction_grid,
      ggplot2::aes(
        x = x_coord,
        y = y_coord,
        fill = residual
      )
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "Residual"
    ) +
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(
        x = x_coord,
        y = y_coord
      ),
      colour = "black",
      size = 1,
      alpha = 0.4
    ) +
    ggplot2::coord_sf(
      crs = crs,
      default_crs = crs
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = coords[1],
      y = coords[2],
      title = "Spatial smooth of residuals",
      subtitle = smooth_label
    )

  # --------------------------------------------------
  # 2. Site-level residuals
  # --------------------------------------------------

  sites_plot <- ggplot2::ggplot(
    plot_data,
  ) + ggplot2::geom_sf(
    data = sf::st_as_sf(prediction_boundary),
    fill = NA,
    colour = "black",
    linewidth = 0.6
  ) +
    ggplot2::geom_point(ggplot2::aes(
      x = x_coord,
      y = y_coord,
      colour = residual_sign,
      size = abs(residual)
    ),
    alpha = 0.8) +
    ggplot2::scale_colour_manual(
      values = c(
        "Negative" = "#2166AC",
        "Zero" = "grey70",
        "Positive" = "#B2182B"
      ),
      name = "Residual"
    ) +
    ggplot2::scale_size_continuous(
      range = c(0.5, 3),
      name = "Absolute residual"
    ) +
    ggplot2::coord_sf(
      crs = crs,
      default_crs = crs
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = coords[1],
      y = coords[2],
      title = "Site-level residuals"
    )

  # Return both plots
  list(
    smooth = smooth_plot,
    sites = sites_plot
  )
}
