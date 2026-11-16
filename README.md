 <p align="center">
  <img src="images/occuFitR_hex.png" width="200"/>
</p>

# occuFitR

**Goodness-of-fit tests and diagnostic tools for single-season occupancy models in R**

`occuFitR` is an R package providing goodness-of-fit (GOF) tests and diagnostic tools for single-season occupancy models fitted using [`unmarked`](https://rbchan.github.io/unmarked/).

The package implements three complementary GOF statistics within a parametric bootstrap framework, allowing model fit to be assessed while accounting for uncertainty in fitted occupancy and detection probabilities.

## Features

`occuFitR` provides three complementary goodness-of-fit statistics:

- **MacKenzie–Bailey (MB) statistic:** An established goodness-of-fit test for occupancy models (MacKenzie and Bailey, 2004).
- **Occupancy and detection probability (\(\psi-p\)) statistic:** Assesses discrepancies between observed and expected detection frequencies across sites, accounting for variation in occupancy and detection probabilities.
- **Detection probability (\(p\)) statistic:** Assesses discrepancies between observed and expected detection frequencies conditional on a species having been detected at a site.

All three statistics are implemented within a consistent parametric bootstrap framework. This procedure simulates detection histories using fitted probabilities, retains the observed pattern of missing visits, refits the occupancy model to each simulated dataset, and calculates bootstrap distributions and p-values.

The package also provides diagnostic tools to investigate model fit:

- **Bootstrap diagnostic plots:** Visualise bootstrap distributions relative to observed GOF statistics.
- **Residuals against site covariates:** Explore relationships between site-level residuals and numeric or categorical covariates.
- **Spatial residual plots:** Investigate spatial patterns in site-level residuals using two-dimensional GAM smooths and individual site-level residuals.

## Installation

The package is currently under development and can be installed from GitHub:

```r
# install.packages("remotes")
remotes::install_github("racheldrake/occuFitR")
```

## Getting started

For a worked introduction to `occuFitR`, including a simulated occupancy example, running the GOF tests, interpreting results and producing diagnostic plots, see the **Basic usage** vignette:

```r
browseVignettes("occuFitR")
```

The vignette is also available online in the package repository under [`vignettes/basic-usage.Rmd`](vignettes/basic-usage.Rmd).

## Requirements

`occuFitR` builds on the following R packages:

- [`unmarked`](https://rbchan.github.io/unmarked/) – fitting occupancy models and extracting model predictions.
- [`ggplot2`](https://ggplot2.tidyverse.org/) – generating diagnostic plots.
- [`dplyr`](https://dplyr.tidyverse.org/) and [`tidyr`](https://tidyr.tidyverse.org/) – data manipulation.
- [`future.apply`](https://future.apply.futureverse.org/) – bootstrap replication.
- [`mgcv`](https://cran.r-project.org/package=mgcv) – fitting GAM smooths for diagnostic plots.
- [`sf`](https://r-spatial.github.io/sf/) – handling spatial coordinates and polygon boundaries.

## Development status

`occuFitR` is under active development. Functions and documentation may change as the package develops.

## Authors

**Rachel Drake (Maintainer)**  
University of St Andrews  
[rld3@st-andrews.ac.uk](mailto:rld3@st-andrews.ac.uk)

**Alison Johnston (Author)**  
University of St Andrews

**Ben Stevenson (Author)**  
University of St Andrews

**David Borchers (Author)**  
University of St Andrews

## Licence

This package is licensed under the **GNU General Public License v3.0 (GPL-3.0)**. You may redistribute and/or modify it under the terms of this licence. See the [LICENSE](LICENSE) file for details.
