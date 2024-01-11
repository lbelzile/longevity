
<!-- README.md is generated from README.Rmd. Please edit that file -->

# longevity <img src="man/figures/longevity_sticker.png" align="right"/>

<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version-last-release/longevity)](https://cran.r-project.org/package=longevity)
[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%203%29-blue.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![Downloads](http://cranlogs.r-pkg.org/badges/longevity?color=brightgreen)](http://www.r-pkg.org/pkg/longevity)
<!-- badges: end -->

The `longevity` package proposes estimation routines for modeling excess
lifetime. Core functionalities include maximum likelihood estimation for
parametric models (exponential, Gompertz, Weibull, generalized Pareto,
extended generalized Pareto, piecewise generalized Pareto), threshold
selection plots for survival data, nonparametric maximum likelihood
estimation, profile likelihood estimation for the endpoint of the
distribution of exceedances.

## Installation

<!-- You can install the released version of longevity from [CRAN](https://CRAN.R-project.org) with: -->
<!-- ``` r -->
<!-- install.packages("longevity") -->
<!-- ``` -->

You can install the development version of longevity from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("lbelzile/longevity")
```

<!-- `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/master/examples>. -->

## Features

- maximum likelihood estimation routines
- simulation of left-truncated and right-truncated/right-censored data
- hazard plots with profile-likelihood based confidence intervals
- threshold selection diagnostics with profile and Wald pointwise
  confidence intervals
- quantile-quantile plots
- likelihood ratio tests for nested models
- likelihood ratio tests for a categorical explanatory
- score and likelihood ratio tests for piecewise generalized Pareto
  distribution, extending Northrop and Coleman (2014), with $p$-value
  paths
- nonparametric maximum likelihood estimate of the distribution function
  with arbitrary truncation and censoring using the EM algorithm of
  Turnbull (1976) - C++ implementation.
- profile likelihood for generalized Pareto (endpoint)
- hazard functions for all parametric models, using `helife`
