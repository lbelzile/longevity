% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parametric.R
\name{fit_elife}
\alias{fit_elife}
\title{Fit excess lifetime models}
\usage{
fit_elife(
  dat,
  thresh,
  ltrunc = NULL,
  rtrunc = NULL,
  rcens = NULL,
  type = c("none", "ltrc", "ltrt"),
  family = c("exp", "gp", "weibull", "gomp", "extgp", "gppiece"),
  weights = rep(1, length(dat)),
  export = FALSE
)
}
\arguments{
\item{dat}{vector of threshold exceedances}

\item{thresh}{threshold}

\item{ltrunc}{lower truncation limit, possibly zero}

\item{rtrunc}{upper truncation limit}

\item{rcens}{logical indicating right-censoring (\code{TRUE} for censored)}

\item{type}{string, either \code{ltrt} for left- and right-truncated data or \code{ltrc} for left-truncated right-censored data}

\item{family}{string; choice of parametric family, either exponential (\code{exp}), Weibull (\code{weibull}), generalized Pareto (\code{gp}), Gompertz (\code{gomp}) or extended generalized Pareto (\code{extgp}).}

\item{weights}{weights for observations}

\item{export}{logical; should data be included in the returned object to produce diagnostic plots? Default to \code{FALSE}.}
}
\value{
an object of class \code{elife_par}
}
\description{
This function is a wrapper around constrained optimization
routines for different models with interval truncation and
left-truncation and right-censoring.
}
\note{
The extended generalized Pareto model is constrained
to avoid regions where the likelihood is flat so \eqn{\xi \in [-1, 10]} in
the optimization algorithm.
}
