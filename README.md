
<!-- README.md is generated from README.Rmd. Please edit that file -->

# longevity <img src="tools/longevity_sticker.png" align="right" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

The `longevity` packages proposes estimation routines for modelling
excess lifetime of old people. Functionalities include

Parametric models for excess lifetime (exponential, Gompertz, Weibull,
generalized Pareto, extended generalized Pareto, piecewise generalized
Pareto), with

-   [x] maximum likelihood estimation routines
-   [x] simulation of left-truncated and right-truncated/right-censored
    data
-   [x] hazard plots with profile-likelihood based confidence intervals
-   [x] threshold selection diagnostics with profile and Wald pointwise
    confidence intervals
-   [x] quantile-quantile plots
-   [x] likelihood ratio tests for nested models
-   [x] likelihood ratio tests for a categorical explanatory
-   [x] score tests for piecewise generalized Pareto distribution,
    extending Northrop and Coleman (2014)
-   [x] nonparametric maximum likelihood estimate of the distribution
    function with arbitrary truncation and censoring using the EM
    algorithm of Turnbull (1976).

**TODO** list:

-   C++ implementation of the nonparametric MLE for the distribution
    function of Turnbull (1976)
-   profile likelihood for generalized Pareto (endpoint)
-   uncertainty for diagnostic plots (via bootstrap)
-   local hazard models (with associated plots)
-   bootstrap *p*-values for nested models (`anova`)
-   unit tests
-   vignettes

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
