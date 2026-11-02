<p align="center">
  <img src="images/hex-sticker.png" width="200">
</p>

# occuFitR

**Goodness-of-fit tools for single-season occupancy models in R**

`occuFitR` is an R package providing goodness-of-fit (GOF) tests and diagnostic tools for single-season occupancy models fitted using [`unmarked`](https://rbchan.github.io/unmarked/).

The package implements three complementary GOF statistics to assess discrepancies between observed detection histories and those expected under a fitted occupancy model. Parametric bootstrap procedures are used to generate reference distributions for the test statistics, allowing model fit to be assessed while accounting for the fitted occupancy and detection probabilities.

## Features

`occuFitR` provides three complementary goodness-of-fit (GOF) statistics for single-season occupancy models:

- **MacKenzie–Bailey (MB) statistic:** An established goodness-of-fit test for occupancy models, originally described by MacKenzie and Bailey (2004), with an existing implementation in AICcmodavg. The implementation in occuFitR produces identical results to AICcmodavg and is included to provide this established test alongside the additional GOF statistics within a common parametric bootstrap framework.
- **Occupancy and detection probability (\(\psi-p\)) statistic:** Compares observed and expected numbers of detections across sites, accounting for variation in both occupancy and detection probabilities.
- **Detection probability (\(p\)) statistic:** Compares observed and expected numbers of detections conditional on a species having been detected at a site.

The package also provides a parametric bootstrap procedure that applies all three statistics within a consistent framework. This procedure:

- Simulates detection histories using fitted occupancy and detection probabilities.
- Retains the observed pattern of missing visits.
- Refits the occupancy model to each simulated dataset.
- Calculates bootstrap distributions and p-values for all three GOF statistics.
- Produces diagnostic plots of bootstrap distributions and site-level residuals against site covariates.

The main contribution of `occuFitR` is therefore the implementation of the additional \(\psi-p\) and \(p\) statistics, together with a unified bootstrap and diagnostic framework for these tests.

## Installation

The package is currently under development. Once available on GitHub, it can be installed using:

```r
# install.packages("remotes")
remotes::install_github("racheldrake/occuFitR")
```

## Usage

`occuFitR` is designed to work with single-season occupancy models fitted using `unmarked`.

First, fit an occupancy model:

```r
library(unmarked)
library(occuFitR)

# Example occupancy model
fit <- occu(
  formula = ~ 1 ~ 1,
  data = umf
)
```

Run the goodness-of-fit procedure:

```r
gof <- occu_gof(
  fit,
  n_boot = 1000,
  seed = 123
)
```

The number of bootstrap replicates and random seed can be adjusted. Set `plot = FALSE` to return the numerical results without generating diagnostic plots.

## Results

The `occu_gof()` function returns a list containing:

| Component | Description |
|---|---|
| `statistic` | Observed values of the three GOF statistics |
| `p_value` | Parametric bootstrap p-values for each statistic |
| `bootstrap` | Bootstrap distributions of the three statistics |
| `true_residuals` | Observed site-level residuals for the \(\psi-p\) and \(p\) tests |
| `plot` | Diagnostic plots for the bootstrap distributions and site-level residuals |

For example:

```r
# Observed GOF statistics
gof$statistic

# Bootstrap p-values
gof$p_value

# Bootstrap distributions
head(gof$bootstrap)

# Individual diagnostic plots
gof$plot$mb
gof$plot$psi_p
gof$plot$p

# Site-level residual plots
gof$plot$site_psi_p
gof$plot$site_p
```

Small bootstrap p-values indicate that the observed GOF statistic is unusually large relative to the distribution expected under the fitted model, providing evidence of lack of fit.

## Requirements

The package builds on the following R packages:

- [`unmarked`](https://rbchan.github.io/unmarked/) – fitting occupancy models and extracting model predictions.
- [`ggplot2`](https://ggplot2.tidyverse.org/) – generating diagnostic plots.
- [`dplyr`](https://dplyr.tidyverse.org/) and [`tidyr`](https://tidyr.tidyverse.org/) – data manipulation.
- [`future.apply`](https://future.apply.futureverse.org/) – bootstrap replication.

Additional dependencies may be required as development continues.

## Development status

`occuFitR` is under active development. Functions and documentation may change as the package develops.

## Authors

**Rachel Drake (Maintainer)**  
University of St Andrews  
rld3@st-andrews.ac.uk

**Alison Johnston (Author)**  
University of St Andrews  

**Ben Stevenson (Author)**  
University of St Andrews 

**David Borchers (Author)**  
University of St Andrews 

## Licence

This package is licensed under the **GNU General Public License v3.0 (GPL-3.0)**. You may redistribute and/or modify it under the terms of this licence. See the [LICENSE](LICENSE) file for details.
