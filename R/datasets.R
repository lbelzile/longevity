#' Dutch survival data
#'
#' This data frame contains information about all Dutch who died
#' above age 92 years between 1986 and 2015. Observations are
#' doubly truncated and such bounds are calculated based on the
#' range of plausible values for these variables.
#' There are 226 records that are interval-censored and interval-truncated
#' for which \code{bdate}, \code{ddate} and \code{ndays} is missing (\code{NA}).
#'
#' @format A data frame with 305143 rows and 11 variables:
#' \describe{
#' \item{ndays}{survival time (in days)}
#' \item{bdate}{the smallest plausible birth date given information about month of birth and death and survival (in days)}
#' \item{bmonth}{month of birth}
#' \item{byear}{year of birth}
#' \item{ddate}{the largest plausible death date given information about month of birth and death and survival (in days)}
#' \item{dmonth}{month of death}
#' \item{dyear}{year of death}
#' \item{ltrunc}{minimum age (in days); the maximum of either 92 years or the number of days reached in 1986}
#' \item{rtrunc}{maximum age (in days) an individual could have reached by the end of 2015}
#' \item{gender}{factor indicating gender of individual, either \code{female} or \code{male}}
#' \item{valid}{quality flag; \code{A} for individuals born in the Netherlands, \code{B} for individuals born abroad who died in the Netherlands}
#' }
#' @references Einmahl, J.J., J.H.J. Einmahl and L. de Haan (2019). \emph{Limits to Human Life Span Through Extreme Value Theory}, Journal of the American Statistical Association, \bold{114}(527), 1075-1080.
#'
#' @source https://doi.org/10.1080/01621459.2018.1537912
"dutch"