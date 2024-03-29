% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{.turnbull_em_dtrunc}
\alias{.turnbull_em_dtrunc}
\title{Turnbull EM algorithm (low storage implementation)}
\usage{
.turnbull_em_dtrunc(
  tsets,
  lcens,
  rcens,
  ltrunc,
  rtrunc,
  weights,
  cens = TRUE,
  trunc = TRUE,
  tol = 1e-12,
  zerotol = 1e-10,
  maxiter = 100000L
)
}
\arguments{
\item{tsets}{Turnbull's sets}

\item{lcens}{numeric vector of left censoring}

\item{rcens}{numeric vector of right censoring}

\item{ltrunc}{numeric vector of left truncation}

\item{rtrunc}{numeric vector of right truncation}

\item{weights}{vector of weights for observations}

\item{cens}{logical; if \code{FALSE}, then \code{censUpp = censLow} and a particular update can be avoided in the EM algorithm}

\item{tol}{tolerance level for terminating the EM algorithm}

\item{maxiter}{maximum number of iteration for the EM algorithm}
}
\value{
a list with the probabilities and the standard errors
}
\description{
Turnbull EM algorithm (low storage implementation)
}
\keyword{internal}
