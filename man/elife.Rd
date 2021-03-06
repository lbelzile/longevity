% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parametric.R
\name{elife}
\alias{elife}
\alias{qelife}
\alias{pelife}
\alias{relife}
\title{Excess lifetime distributions}
\usage{
qelife(
  p,
  scale,
  shape,
  family = c("exp", "gp", "weibull", "gomp", "extgp"),
  lower.tail = TRUE
)

pelife(
  q,
  scale,
  shape,
  family = c("exp", "gp", "weibull", "gomp", "extgp"),
  lower.tail = TRUE,
  log.p = FALSE
)

relife(n, scale, shape, family = c("exp", "gp", "weibull", "gomp", "extgp"))
}
\arguments{
\item{p}{vector of probabilities}

\item{scale}{scale parameter, strictly positive.}

\item{shape}{vector of shape parameter(s).}

\item{family}{string indicating the parametric model, one of \code{exp}, \code{gp}, \code{gomp}, \code{weibull} and \code{extgp}}

\item{lower.tail}{logical; if \code{TRUE} (default), the lower tail probability \eqn{\Pr(X \leq x)} is returned.}

\item{q}{vector of quantiles.}

\item{log.p}{logical; if \code{FALSE} (default), values are returned on the probability scale.}
}
\description{
Quantile and distribution function of excess lifetime distribution
for threshold exceedances.
}
\keyword{internal}
