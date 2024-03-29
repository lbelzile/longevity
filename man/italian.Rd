% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{italian}
\alias{italian}
\title{Italian semi-supercentenarian}
\format{
A data frame with 3836 rows and 6 variables:
\describe{
\item{ndays}{survival time (in days)}
\item{bdate}{birth date (\code{Date})}
\item{ddate}{death date (\code{Date}), or \code{NA_Date_}
if the person is alive at the end of the sampling}
\item{ltrunc}{minimum age (in days); the maximum of 38351 days
(approximately 105 years) or the number of days reached in 2009}
\item{event}{integer indicating the censoring pattern; \code{0} for
right-censored records, \code{1} for fully observed}
\item{gender}{factor indicating gender of individual,
either \code{female} or \code{male}}
}
}
\usage{
italian
}
\description{
This data frame contains information about 3836 Italians
individually validated survival lifetimes times in days
of all persons in Italy who were at least 105 years old
at some point in the period from 1 January 2009
to 31 December 2015.
Observations are left-truncated and right-censored.
These data are not publicly available, but can be purchased
from the Italian National Institute of Statistics by
registering at the Contact Center and mentioning the
Semi-supercentenarian Survey and Marco Marsili
as contact person.
}
\references{
Istituto Nazionale di Statistica
}
\keyword{datasets}
