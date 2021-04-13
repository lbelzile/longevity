#' Simulate excess lifetime with truncation or right-censoring
#'
#' This function dispatches simulations accounting for potential left-truncation (remove by setting lower to zero).
#' If \code{type2=ltrt}, simulated observations will be lower than the upper bounds \code{upper}.
#' If \code{type2=ltrc}, simulated observations are capped at \code{upper} and the observation is right-censored (\code{rcens=TRUE}).
#'
#' @param n sample size
#' @param scale scale parameter(s)
#' @param shape shape parameter(s)
#' @param family string; choice of parametric family, either exponential (\code{exp}), Weibull (\code{weibull}), generalized Pareto (\code{gp}), Gompertz (\code{gomp}), Gompertz-Makeham (\code{gompmake}) or extended generalized Pareto (\code{extgp}).
#' @param lower vector of lower bounds
#' @param upper vector of upper bounds
#' @param type2 string, either \code{none}, \code{ltrt} for left- and right-truncated data or \code{ltrc} for left-truncated right-censored data
#' @export
#' @return either a vector of observations or, if \code{type2=ltrc}, a list with \code{n} observations \code{dat} and a logical vector of the same length with \code{TRUE} for right-censored observations and \code{FALSE} otherwise.
samp_elife <- function(n,
                       scale,
                       shape = NULL,
                       lower = 0,
                       upper = Inf,
                       family = c("exp","gp","gomp","gompmake","weibull","extgp","gppiece"),
                       type2 = c("none","ltrt","ltrc","ditrunc")){
  family <- match.arg(family)
  if(family != "exp" & is.null(shape)){
    stop("Invalid shape parameter.")
  }
  type2 <- match.arg(type2)
  if(type2 == "none"){
    relife(n = n,
           scale = scale,
           shape = shape,
           family = family)
  } else if(type2 == "ltrt"){
    r_dtrunc_elife(n = n,
                   scale = scale,
                   shape = shape,
                   lower = lower,
                   upper = upper,
                   family = family)
  } else if(type2 == "ltrc"){
    r_ltrc_elife(n = n,
                 scale = scale,
                 shape = shape,
                 lower = lower,
                 upper = upper,
                 family = family)
  } else if(type2 == "ditrunc"){
    r_ditrunc_elife(n = n,
                    scale = scale,
                    shape = shape,
                    lower = lower,
                    upper = upper,
                    family = family)
  }
}



#' Sample observations from an interval truncated excess lifetime distribution
#'
#' @inheritParams samp_elife
#' @return a vector of \code{n} observations
#' @export
#' @keywords internal
#' @examples
#' n <- 100L
#' # the lower and upper bound are either scalar,
#' # or else vectors of length n
#' r_dtrunc_elife(n = n, scale = 1, shape = -0.1,
#'                lower = pmax(0, runif(n, -0.5, 1)),
#'                upper = runif(n, 6, 10),
#'                family = "gp")
r_dtrunc_elife <- function(n,
                          scale,
                          shape,
                          lower,
                          upper,
                          family = c("exp","gp","gomp","gompmake","weibull","extgp")
                          ){
  family <- match.arg(family)
  if(length(lower) > 1){
    stopifnot("Sample size should match length of 'lower'" = length(lower) == n)
  }
  if(length(upper) > 1){
    stopifnot("Sample size should match length of 'upper'" = length(upper) == n)
  }
  if(family == "gompmake"){
    stopifnot("Scale and shape parameters must be positive" = isTRUE(all(c(scale[1] > 0, scale[2] >= 0, shape[1] >= 0))))
    if(isTRUE(all.equal(shape[1], 0, ignore.attributes = TRUE))){
      family == "gp"
      shape[1] <- 0
    } else if(isTRUE(all.equal(scale[2], 0, ignore.attributes = TRUE))){
      family == "gomp"
      scale <- scale[1]
    } else{
      return(qgompmake(pgompmake(q = lower, scale = scale[1], lambda = scale[2], shape = shape[1]) +
                     runif(n)*(pgompmake(upper,  scale = scale[1], lambda = scale[2], shape = shape[1]) - pgompmake(lower, scale = scale[1], lambda = scale[2], shape = shape[1])),
                     scale = scale[1], lambda = scale[2], shape = shape[1])
      )
    }
  }
  if(family == "gomp"){
    stopifnot("Scale and shape parameters must be positive" = isTRUE(scale > 0 && shape[1] >= 0))
    if(isTRUE(all.equal(shape[1], 0, ignore.attributes = TRUE))){
      family == "gp"
      shape[1] <- 0
    } else{
      return(qgomp(pgomp(lower, scale = scale, shape = shape[1]) +
                     runif(n)*(pgomp(upper, scale = scale, shape = shape[1]) - pgomp(lower, scale = scale, shape = shape[1])),
                   scale = scale, shape = shape[1])
              )
    }
  }
  if(family == "exp"){
    stopifnot("Scale parameter must be positive" = isTRUE(scale > 0))
    shape <- 0
    family <- "gp"
  }
  if(family == "gp"){
    stopifnot("Scale and shape parameters must be in range of admissible values" = isTRUE(scale > 0 && shape[1] >= -1)
              # The following warning is discontinued because the distribution function evaluates to 1 beyond the support
              # "Upper bound must be lower than the maximum value of the support." = ifelse(shape[1] < 0, max(upper) < -scale/shape[1], TRUE)
    )
    return(qgpd(pgpd(lower, scale = scale, shape = shape[1]) +
         runif(n)*(pgpd(upper, scale = scale, shape = shape[1]) - pgpd(lower, scale = scale, shape = shape[1])),
       scale = scale, shape = shape[1])
  )
  } else if(family == "extgp"){
    stopifnot("Scale and shape parameters must be in range of admissible values" = isTRUE(scale > 0 && shape[1] >= 0),
              "Shape parameters must be in the range of admissible values" = 1 - shape[1] / shape[2] > 0,
              "Upper bound must be lower than the maximum value of the support." = ifelse(shape[2] < 0, max(upper) < scale / shape[1] * log(1 - shape[1] / shape[2]), TRUE)
    )
    return(qextgp(pextgp(lower, scale = scale, shape1 = shape[1], shape2 = shape[2]) +
                  runif(n)*(pextgp(upper, scale = scale, shape1 = shape[1], shape2 = shape[2]) - pextgp(lower, scale = scale, shape1 = shape[1], shape2 = shape[2])),
                scale = scale, shape1 = shape[1], shape2 = shape[2])
    )
  } else if(family == "weibull"){
    stopifnot("Scale and shape parameters must be in range of admissible values" = isTRUE(scale > 0 && shape[1] > 0))
    return(qweibull(pweibull(lower, scale = scale, shape = shape[1]) +
                  runif(n)*(pweibull(upper, scale = scale, shape = shape[1]) - pweibull(lower, scale = scale, shape = shape[1])),
                scale = scale, shape = shape[1])
    )
  }
}



#' Sample observations from an interval truncated excess lifetime distribution
#'
#' @inheritParams samp_elife
#' @return a vector of \code{n} observations
#' @export
#' @keywords internal
#' @examples
#' n <- 100L
#' # the lower and upper bound are either scalar,
#' # or else vectors of length n
#' r_dtrunc_elife(n = n, scale = 1, shape = -0.1,
#'                lower = pmax(0, runif(n, -0.5, 1)),
#'                upper = runif(n, 6, 10),
#'                family = "gp")
r_ditrunc_elife <- function(n,
                           scale,
                           shape,
                           lower,
                           upper,
                           family = c("exp","gp","gomp","gompmake","weibull","extgp")
){
  stopifnot("`lower` should be a matrix." = is.matrix(lower),
            "`upper` should be a matrix." = is.matrix(upper),
            "`lower` must be a matrix with 2 columns." = ncol(lower) == 2,
            "`upper` must be a matrix with 2 columns." = ncol(upper) == 2,
            "`lower` and `upper` should have the same number of rows." = nrow(lower) == nrow(upper),
            "`lower` and `upper` should have 1 or n rows." = nrow(upper) == 1 | nrow(upper) == n)
  family <- match.arg(family)
  if(family == "exp"){
    quantf <- function(par, dat){
      qexp(p = dat, rate = 1/par[1])
    }
    cdf <- function(par, dat, lower.tail = TRUE, log.p = FALSE){
      pexp(q = dat, rate = 1/par[1], lower.tail = lower.tail, log.p = log.p)
    }
  } else if(family == "gp"){
    quantf <- function(par, dat){
      qgpd(p = dat, loc = 0, scale = par[1], shape = par[2])
      }
    cdf <- function(par, dat, lower.tail = TRUE, log.p = FALSE){
      pgpd(q = dat, loc = 0, scale = par[1], shape = par[2], lower.tail = lower.tail, log.p = log.p)
      }
  } else if(family == "weibull"){
    quantf <- function(par, dat){
      qweibull(p = dat, scale = par[1], shape = par[2])
      }
    cdf <- function(par, dat, lower.tail = TRUE, log.p = FALSE){
      pweibull(q = dat, scale = par[1], shape = par[2], lower.tail = lower.tail, log.p = log.p)
      }
  } else if(family == "extgp"){
    quantf <- function(dat, par){
      qextgp(p = dat, scale = par[1], shape1 = par[2], shape2 = par[3])
    }
    cdf <- function(dat, par, lower.tail = TRUE, log.p = FALSE){
      pextgp(q = dat, scale = par[1], shape1 = par[2], shape2 = par[3], lower.tail = lower.tail, log.p = log.p)
    }
  } else if(family == "gompmake"){
    quantf <- function(dat, par){
      qgompmake(p = dat, scale = par[1], shape = par[2], lambda = par[3])
    }
    cdf <- function(dat, par, lower.tail = TRUE, log.p = FALSE){
      pgompmake(q = dat, scale = par[1], shape = par[2], lambda = par[3], lower.tail = lower.tail, log.p = log.p)
    }
  }
  # Not all records are doubly truncated records...
  # If there are NAs, then second interval has prob zero

  if(family == "gompmake"){
    par <- c(scale[1], shape, scale[2])
  } else{
    par <- c(scale, shape)
  }
  probint1 <- cdf(dat = upper[,1], par = par) - cdf(dat = lower[,1], par = par)
  probint2 <- ifelse(is.na(lower[,2]), 0,
                     cdf(upper[,2], par = par) - cdf(lower[,2], par = par))
  if(isTRUE(any((probint1 + probint2) > 1))){
    stop("Invalid input")
  }

  cutoff <- probint1 / (probint1 + probint2)
  unif <- runif(n)
  quantf(ifelse(unif < cutoff,
         unif*(probint1 + probint2) + cdf(lower[,1], par = par),
         unif*(probint1 + probint2) + cdf(lower[,2], par = par) - probint1),
    par = par)
}


#' Sample observations from a left-truncated right-censored excess lifetime distribution
#'
#' @param n sample size
#' @param scale scale parameter
#' @param shape shape parameter(s)
#' @param family string; choice of parametric family, either exponential (\code{exp}), Weibull (\code{weibull}), generalized Pareto (\code{gp}), Gompertz (\code{gomp}) or extended generalized Pareto (\code{extgp}).
#' @param lower vector of lower bounds
#' @param upper vector of upper bounds above which data are right-truncated
#' @return a list with \code{n} observations \code{dat} and a logical vector of the same length with \code{TRUE} for right-censored observations and \code{FALSE} otherwise.
#' @export
#' @keywords internal
#' @examples
#' n <- 100L
#' # the lower and upper bound are either scalar,
#' # or else vectors of length n
#' r_ltrc_elife(n = n, scale = 5, shape = -0.1,
#'                lower = pmax(0, runif(n, -0.5, 1)),
#'                upper = 5,
#'                family = "gp")
r_ltrc_elife <- function(n,
                         scale,
                         shape,
                         lower,
                         upper,
                         family = c("exp","gp","gomp","gompmake","weibull","extgp")
){
  family <- match.arg(family)
  if(length(lower) > 1){
    stopifnot("Sample size should match length of 'lower'" = length(lower) == n)
  }
  if(length(upper) > 1){
    stopifnot("Sample size should match length of 'upper'" = length(upper) == n)
  }
  if(family == "gompmake"){
    stopifnot("Scale and shape parameters must be positive" = isTRUE(all(c(scale[1] > 0, scale[2] >=0, shape >= 0))))
    if(isTRUE(all.equal(shape[1], 0, ignore.attributes = TRUE))){
      family == "gp"
      shape[1] <- 0
    } else if(isTRUE(all.equal(scale[2], 0, ignore.attributes = TRUE))){
      family == "gomp"
      scale <- scale[1]
    } else{
      samp <- qgompmake(pgompmake(lower, scale = scale[1], lambda = scale[2], shape = shape[1]) +
                      runif(n)*(1 - pgompmake(lower, scale = scale[1], lambda = scale[2], shape = shape[1])),
                    scale = scale[1], lambda = scale[2], shape = shape[1])

    }
  }
  if(family == "gomp"){
    stopifnot("Scale and shape parameters must be positive" = isTRUE(scale > 0 && shape[1] >= 0))
    if(isTRUE(all.equal(shape[1], 0, ignore.attributes = TRUE))){
      family == "gp"
      shape[1] <- 0
    } else{
      samp <- qgomp(pgomp(lower, scale = scale, shape = shape[1]) +
                     runif(n)*(1 - pgomp(lower, scale = scale, shape = shape[1])),
                   scale = scale, shape = shape[1])

    }
  }

  if(family == "exp"){
    stopifnot("Scale parameter must be positive" = isTRUE(scale > 0))
    shape <- 0
    family <- "gp"
  }
  if(family == "gp"){
    stopifnot("Scale and shape parameters must be in range of admissible values" = isTRUE(scale > 0 && shape[1] >= -1)
              # The following warning is discontinued because the distribution function evaluates to 1 beyond the support
              # "Upper bound must be lower than the maximum value of the support." = ifelse(shape[1] < 0, max(upper) < -scale/shape[1], TRUE)
    )
    samp <- qgpd(pgpd(lower, scale = scale, shape = shape[1]) +
                  runif(n)*(1 - pgpd(lower, scale = scale, shape = shape[1])),
                scale = scale, shape = shape[1])
  } else if(family == "extgp"){
    stopifnot("Scale and shape parameters must be in range of admissible values" = isTRUE(scale > 0 && shape[1] >= 0),
              "Shape parameters must be in the range of admissible values" = 1 - shape[1] / shape[2] > 0,
              "Upper bound must be lower than the maximum value of the support." = ifelse(shape[2] < 0, max(upper) < scale / shape[1] * log(1 - shape[1] / shape[2]), TRUE)
    )
    samp <- qextgp(pextgp(lower, scale = scale, shape1 = shape[1], shape2 = shape[2]) +
                    runif(n)*(1 - pextgp(lower, scale = scale, shape1 = shape[1], shape2 = shape[2])),
                  scale = scale, shape1 = shape[1], shape2 = shape[2])
  } else if(family == "weibull"){
    stopifnot("Scale and shape parameters must be in range of admissible values" = isTRUE(scale > 0 && shape[1] > 0))
    samp <- qweibull(pweibull(lower, scale = scale, shape = shape[1]) +
                      runif(n)*(1 - pweibull(lower, scale = scale, shape = shape[1])),
                    scale = scale, shape = shape[1])
  }
  rcens <- ifelse(samp < upper, FALSE, TRUE)
  samp <- pmin(samp, upper)
  return(list(dat = samp,
              rcens = rcens)
         )
}

#' Sample lifetime from excess lifetime model
#'
#' Given parameters of a \code{elife} distribution, sampling window and
#' birth dates with excess lifetimes, sample new observations; excess lifetime
#' at \code{c1} are sampled from an exponential distribution, whereas
#' the birth dates are sampled from a jittered histogram-based distribution
#' The new excess lifetime above the threshold are right-censored if they exceed
#' \code{c2}.
#'
#' @inheritParams r_dtrunc_elife
#' @param xcal date at which individual reaches \code{u} years
#' @param c1 date, first day of the sampling frame
#' @param c2 date, last day of the sampling frame
#' @return list with new birthdates (\code{xcal}), excess lifetime at \code{c1} (\code{ltrunc}),
#' excess lifetime above \code{u} (\code{dat}) and right-censoring indicator (\code{rightcens}).
#' @export
#' @keywords internal
samp2_elife <- function(n,
                         scale,
                         shape,
                         family = c("exp","gp","gomp","gompmake","weibull","extgp"),
                         xcal,
                         c1,
                         c2){
  stopifnot("Date `c1` is missing" = !missing(c1),
            "Date `c2` is missing" = !missing(c2),
            "`c1` should be a single date" = is.Date(c1) && length(c1) == 1L,
            "`c2` should be a single date" = is.Date(c2) && length(c2) == 1L)
  ltrunc <- as.numeric(pmax(0, c1 - xcal))
  sample_dates <- function(n, xcal, c1, c2, ltrunc){
    sample_ltrunc <- function(n, ltrunc){
      sort(round(rexp(n, rate = 1/mean(365.25*ltrunc[ltrunc>0])))/365.25, decreasing = TRUE)
    }
    nday <- as.numeric(xcal[ind]-c1)
    nday <- nday[nday>0]
    nmax <- as.numeric(c2-c1)
    ssltrunc <- sample_ltrunc(round(sum(ltrunc>0)*n/length(ltrunc)), ltrunc = ltrunc)
    xhist <- hist(nday, plot = FALSE)
    bins <- with(xhist, sample(length(mids), n-length(ssltrunc), p=density, replace=TRUE)) # choose a bin
    result <- round(runif(length(bins), xhist$breaks[bins], pmin(xhist$breaks[bins+1], nmax-1)))
    list(xcal = as.Date(round(c(-ssltrunc*365.25, sort(result))), origin = c1),
         ltrunc = c(ssltrunc, rep(0, n-length(ssltrunc))))
  }

  traject <- rdtrunc_elife(n = n,
                           scale = scale,
                           shape = shape,
                           lower = 0,
                           upper = Inf,
                           family = family)
  sdates <- sample_dates(n = n, xcal, c1, c2, ltrunc)
  lifeah <- pmax(1,pmin(round(365.25*(traject + sdates$ltrunc)), as.numeric(c2-sdates$xcal)))
  rcens_new <- lifeah == as.numeric(c2-sdates$xcal)
  list(xcal = sdates$xcal, ltrunc = sdates$ltrunc, dat = lifeah/365.25, rightcens = rcens_new)
}
