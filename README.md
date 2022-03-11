
<!-- README.md is generated from README.Rmd. Please edit that file -->

# longevity <img src="man/figures/longevity_sticker.png" align="right"/>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
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

## Current features

-   [x] maximum likelihood estimation routines
-   [x] simulation of left-truncated and right-truncated/right-censored
    data
-   [x] hazard plots with profile-likelihood based confidence intervals
-   [x] threshold selection diagnostics with profile and Wald pointwise
    confidence intervals
-   [x] quantile-quantile plots
-   [x] likelihood ratio tests for nested models
-   [x] likelihood ratio tests for a categorical explanatory
-   [x] score and likelihood ratio tests for piecewise generalized
    Pareto distribution, extending Northrop and Coleman (2014), with
    *p*-value plots
-   [x] nonparametric maximum likelihood estimate of the distribution
    function with arbitrary truncation and censoring using the EM
    algorithm of Turnbull (1976) - C++ implementation.
-   [x] profile likelihood for generalized Pareto (endpoint)

## Improvements

-   [ ] Bootstrap *p*-values for nested model (`anova`)
-   [ ] Remove/keep Kolmogorov-Smirnov test (depending on whether it
    makes sense)
-   [ ] Adapt `npsurv` for interval censoring
-   [ ] Add plots of (local) hazards with delta-method based confidence
    intervals
-   [ ] Check whether bootstrap scheme (conditional) should be adapted
    to sample truncation patterns too (using Shen’s IPW estimator for
    the nonparametric MLE of the distribution of the truncation bounds.

## Testing

-   [ ] Check `npsurv` for the case of
    1)  left-truncation (Bell-Lynden estimator) and
    2)  double truncation (is the problem in `DTDA`?)
-   [ ] Check functions fail when packages listed in ‘Suggests’ are
    absent
-   [ ] Check that all ANOVA nesting works as expected (with null
    distribution)
-   [ ] Verify fitting procedure in multiple instances, including
    interval censoring, left and right truncation, etc.
-   [ ] Check all plots type are produced with both base **R** and
    `ggplot2`

## Package on CRAN

-   [ ] Add tests and examples for each function
-   [ ] Add vignettes
-   [x] Use `pkgdown` to create a webpage
-   [ ] Root out data sets that cannot go on CRAN (with accompanying
    examples)
-   [ ] Submit to the CRAN
