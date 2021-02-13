
<!-- README.md is generated from README.Rmd. Please edit that file -->

# longevity <img src="tools/longevity_sticker.png" align="right" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

The `longevity` packages proposes estimation routines for modelling
excess lifetime of old people. Functionalities include

Parametric models for excess lifetime (exponential, Gompertz, Weibull,
generalized Pareto, extended generalized Pareto), with

-   maximum likelihood estimation
-   simulation for doubly truncated data
-   hazard plots with profile-likelihood based confidence intervals
-   threshold selection diagnostics
-   likelihood ratio tests for features
-   quantile-quantile plots
-   bootstrap goodness-of-fit tests
-   likelihood ratio tests for nested models

The package will also include tools for esitmation of the nonparametric
hazard estimates of Turnbull (1978) with pointwise confidence intervals.

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
