#' Distribution function of the generalized Pareto distribution
#'
#' @param q vector of quantiles.
#' @param loc location parameter.
#' @param scale scale parameter, strictly positive.
#' @param shape shape parameter.
#' @param lower.tail logical; if \code{TRUE} (default), the lower tail probability \eqn{\Pr(X \leq x)} is returned.
#' @param log.p logical; if \code{FALSE} (default), values are returned on the probability scale.
#' @return a vector of (log)-probabilities of the same length as \code{q}
#' @export
#' @keywords internal
pgpd <- function(q,
                 loc = 0,
                 scale = 1,
                 shape = 0,
                 lower.tail = TRUE,
                 log.p = FALSE){
  if (min(scale) <= 0){
    stop("invalid scale")
  }
  if (length(shape) != 1){
    stop("invalid shape")
  }
  q <- pmax(q - loc, 0)/scale
  if (shape == 0){
    p <- 1 - exp(-q)
  } else {
    p <- 1 - exp((-1/shape)*log1p(pmax(-1, shape * q)))
  }
  if (!lower.tail){
    p <- 1 - p
  }
  if (log.p){
    p <- log(p)
  }
  return(p)
}

#' Density function of the generalized Pareto distribution
#'
#' @inheritParams pgpd
#' @param x vector of quantiles.
#' @param log logical; if \code{FALSE} (default), return the density, else the log likelihood of the individual observations.
#' @return a vector of (log)-density.
#' @export
#' @keywords internal
dgpd <- function (x,
                  loc = 0,
                  scale = 1,
                  shape = 0,
                  log = FALSE){
  if (length(scale) != 1 || scale <= 0){
    stop("invalid scale parameter")
  }
  if (length(shape) != 1){
    stop("invalid shape")
  }
  d <- (x - loc)/scale
  index <- (d > 0 & ((1 + shape * d) > 0)) | is.na(d)
  if (shape == 0) {
    d[index] <- log(1/scale) - d[index]
    d[!index] <- -Inf
  }
  else {
    d[index] <- log(1/scale) - (1/shape + 1) * log1p(shape * d[index])
    d[!index] <- -Inf
  }
  if (!log){
    d <- exp(d)
  }
  return(d)
}

#' Quantile function of the generalized Pareto distribution
#'
#' @param p vector of probabilities.
#' @inheritParams pgpd
#' @return vector of quantiles
#' @export
#' @keywords internal
qgpd <- function(p,
                 loc = 0,
                 scale = 1,
                 shape = 0,
                 lower.tail = TRUE){
  if (min(p, na.rm = TRUE) < 0 || max(p, na.rm = TRUE) > 1)
    stop("`p' must contain probabilities in (0,1)")
  if (length(scale) != 1L || scale <= 0) {
    stop("invalid scale")
  }
  if (length(shape) != 1){
    stop("invalid shape")
  }
  if (lower.tail){
    p <- 1 - p
  }
  if (shape == 0){
    return(loc - scale * log(p))
  } else {
    return(loc + scale * (p^(-shape) - 1)/shape)
  }
}


#' Distribution function of the Gompertz distribution
#'
#' @inheritParams pgpd
#' @export
#' @return a vector of (log)-probabilities of the same length as \code{q}
#' @keywords internal
pgomp <- function(q,
                  scale = 1,
                  shape = 0,
                  lower.tail = TRUE,
                  log.p = FALSE){
  if (min(scale) <= 0){
    stop("invalid scale")
  }
  if (length(shape) != 1 || shape[1] < 0){
    stop("invalid shape")
  }
  if(isTRUE(all.equal(shape, 0, check.attributes = FALSE))){
    # Exponential data
    return(pexp(q = q,
                rate = 1/scale,
                lower.tail = lower.tail,
                log.p = log.p)
    )
  } else {
    p <- 1 - exp(-(exp(shape*q/scale)-1)/shape)
  }
  if (!lower.tail){
    p <- 1 - p
  }
  if (log.p){
    p <- log(p)
  }
  return(p)
}

#' Quantile function of the Gompertz distribution
#'
#' @param p vector of probabilities.
#' @inheritParams pgpd
#' @export
#' @return vector of quantiles
#' @keywords internal
qgomp <- function(p,
                  scale = 1,
                  shape = 0,
                  lower.tail = TRUE){
  if (min(p, na.rm = TRUE) < 0 || max(p, na.rm = TRUE) > 1)
    stop("`p' must contain probabilities in (0,1)")
  if (length(scale) != 1L || scale <= 0) {
    stop("invalid scale")
  }
  if (length(shape) != 1 || shape < 0){
    stop("invalid shape")
  }
  if(isTRUE(all.equal(shape, 0, check.attributes = FALSE))){
    # Exponential data
    return(qexp(p = p,
                rate = 1/scale,
                lower.tail = lower.tail)
    )
  }
  if(!lower.tail){
    p <- 1 - p
  }
  return(scale / shape * log(1 - shape * log(1 - p)))
}

#' Distribution function of the Gompertz-Makeham distribution
#'
#' @inheritParams pgpd
#' @param lambda exponential rate
#' @export
#' @return a vector of (log)-probabilities of the same length as \code{q}
#' @keywords internal
pgompmake <- function(q,
                      scale = 1,
                      shape = 1,
                      lambda = 1,
                      lower.tail = TRUE,
                      log.p = FALSE){
  stopifnot("`scale` should be a vector of length 1." = length(scale) == 1L,
            "`shape` should be a vector of length 1." = length(shape) == 1L,
            "`lambda` should be a vector of length 1." = length(lambda) == 1L,
            "`scale` must be positive." = scale > 0,
            "`shape` must be non-negative." = shape >= 0,
            "`lambda` must be non-negative." = lambda >= 0
  )

  if(isTRUE(all.equal(shape, 0, check.attributes = FALSE))){
    # Exponential data
    return(pexp(q = q,
                rate = lambda + 1/scale,
                lower.tail = lower.tail,
                log.p = log.p)
    )
  } else {
    p <- 1 - exp(-lambda*q - (exp(shape*q/scale)-1)/shape)
  }
  if (!lower.tail){
    p <- 1 - p
  }
  if (log.p){
    p <- log(p)
  }
  return(p)
}

#' Quantile function of the Gompertz-Makeham distribution
#'
#' @param p vector of probabilities.
#' @inheritParams pgpd
#' @param lambda exponential rate
#' @export
#' @return vector of quantiles
#' @keywords internal
qgompmake <- function(p,
                      scale = 1,
                      shape = 1,
                      lambda = 1,
                      lower.tail = TRUE){
  stopifnot("`scale` should be a vector of length 1." = length(scale) == 1L,
            "`shape` should be a vector of length 1." = length(shape) == 1L,
            "`lambda` should be a vector of length 1." = length(lambda) == 1L,
            "`scale` must be positive." = scale > 0,
            "`shape` must be non-negative." = shape >= 0,
            "`lambda` must be non-negative." = lambda >= 0
  )
      stopifnot("Install package \"gsl\" to use \"qgompmake\" with the Gompertz model.\n Try `install.packages(\"gsl\")`" = requireNamespace("gsl", quietly = TRUE))
  if (min(p, na.rm = TRUE) < 0 || max(p, na.rm = TRUE) > 1)
    stop("`p' must contain probabilities in (0,1)")
  if(isTRUE(all.equal(shape, 0, check.attributes = FALSE))){
    # Exponential data - but parameter not identifiable
    return(qexp(p = p,
                rate = lambda + 1/scale,
                lower.tail = lower.tail)
    )
  }
  if(isTRUE(all.equal(lambda, 0, check.attributes = FALSE))){
    # Exponential data - but parameter not identifiable
    return(qgomp(p = p,
                scale = scale,
                shape = shape,
                lower.tail = lower.tail)
    )
  }
  if(!lower.tail){
    p <- 1 - p
  }
  return(1/(shape*lambda) - log(1 - p) / lambda - scale / shape * gsl::lambert_W0(exp(1/(scale*lambda))*(1-p)^(-shape/(scale*lambda)) / (scale*lambda)))
}

#' Density function of the Gompertz-Makeham distribution
#'
#' @param x vector of quantiles.
#' @inheritParams pgpd
#' @param lambda exponential rate
#' @export
#' @return vector of density
#' @keywords internal
dgompmake <- function(x,
                      scale = 1,
                      shape = 1,
                      lambda = 1,
                      log = FALSE){
  stopifnot("`scale` should be a vector of length 1." = length(scale) == 1L,
            "`shape` should be a vector of length 1." = length(shape) == 1L,
            "`lambda` should be a vector of length 1." = length(lambda) == 1L,
            "`scale` must be positive." = scale > 0,
            "`shape` must be non-negative." = shape >= 0,
            "`lambda` must be non-negative." = lambda >= 0
  )
  if(isTRUE(all.equal(shape, 0, check.attributes = FALSE))){
    # Exponential data - but parameter not identifiable
    return(dexp(x = x,
                rate = lambda + 1/scale,
                log = log)
    )
  }
  x <- pmax(0, x)
  ld1 <- log(lambda + exp(shape * x / scale) / scale)
  ld1 <- ifelse(is.finite(ld1), ld1, 0)
  ldens <- ld1 - lambda*x - (exp(shape * x / scale) - 1) / shape
  if(log){
    return(ldens)
  } else{
    return(exp(ldens))
  }
}

#' Distribution function of the extended generalized Pareto distribution
#'
#' @inheritParams pgpd
#' @param shape1 positive shape parameter \eqn{\beta}; model defaults to generalized Pareto when it equals zero.
#' @param shape2 shape parameter \eqn{\gamma}; model reduces to Gompertz when \code{shape2=0}.
#' @return a vector of (log)-probabilities of the same length as \code{q}
#' @export
#' @keywords internal
pextgp <- function(q,
                   scale = 1,
                   shape1,
                   shape2,
                   lower.tail = TRUE,
                   log.p = FALSE){
  if (length(scale) != 1 || min(scale) <= 0){
    stop("invalid scale")
  }
  if (length(shape1) != 1 || shape1 < 0){
    stop("invalid shape1")
  }
  if (length(shape2) != 1){
    stop("invalid shape2")
  }
  if(isTRUE(all.equal(shape1, 0, check.attributes = FALSE))){
    # Exponential data
    return(pgpd(q = q,
                scale = scale,
                shape = shape2,
                lower.tail = lower.tail,
                log.p = log.p)
    )
  } else if(isTRUE(all.equal(shape2, 0, check.attributes = FALSE))){
    return(pgomp(q = q,
                 scale = scale,
                 shape = shape1,
                 lower.tail = lower.tail,
                 log.p = log.p)
    )
  }
  p <- 1 - pmax(0,(1+shape2*(exp(shape1*q/scale)-1)/shape1)^(-1/shape2))
  if (!lower.tail){
    p <- 1 - p
  }
  if (log.p){
    p <- log(p)
  }
  return(p)
}
#' Density function of the Gompertz distribution
#'
#' @param x vector of quantiles
#' @param scale positive scale parameter
#' @param shape non-negative shape parameter
#' @param log logical; if \code{TRUE}, return the log density
#' @keywords internal
dgomp <- function(x,
                  scale = 1,
                  shape,
                  log = FALSE){
  stopifnot("`shape` must be non-negative." = shape >= 0,
            "`scale` must be positive." = scale > 0)
  if(shape < 1e-8){
    ldens <-  -log(scale) + -x/scale
  } else{
    ldens <-  -log(scale) + (shape*x/scale - exp(shape*x/scale)/shape + 1/shape)
  }
  if(log){
    return(ldens)
  } else{
    exp(ldens)
  }
}
#' Density function of the extended generalized Pareto distribution
#'
#' @inheritParams pgpd
#' @param shape1 positive shape parameter \eqn{\beta}; model defaults to generalized Pareto when it equals zero.
#' @param shape2 shape parameter \eqn{\gamma}; model reduces to Gompertz when \code{shape2=0}.
#' @param log logical; if \code{TRUE}, return the log-density
#' @return a vector of (log)-density of the same length as \code{x}
#' @export
#' @keywords internal
dextgp <- function(x,
                   scale = 1,
                   shape1,
                   shape2,
                   log = FALSE){
  stopifnot("`shape1` must be non-negative." = shape1 >= 0,
            "`shape2` must be larger than -1." = shape2 >= -1,
            "`scale` must be positive." = scale > 0)
  if(abs(shape2) < 1e-8 && abs(shape1) < 1e-8){
    ldens <-  -log(scale) + -x/scale
  } else if(abs(shape2) < 1e-8 && abs(shape1) > 1e-8){ #Gompertz
    ldens <-  -log(scale) + (shape1*x/scale - exp(shape1*x/scale)/shape1 + 1/shape1)
  } else if(abs(shape2) >= 1e-8 && abs(shape1) < 1e-8){ #generalized Pareto
    ldens <-  - log(scale) - (1/shape2 + 1)*log(pmax(0, 1+x*shape2/scale))
  } else { #extended
    ldens <- -log(scale) + (-1/shape2 - 1)*log(shape2*(exp(shape1*x/scale) - 1)/shape1 + 1) + shape1*x/scale
  }
  if(log){
    return(ldens)
  } else{
    exp(ldens)
  }
}

#' Quantile function of the extended generalized Pareto distribution
#'
#' @param p vector of probabilities.
#' @inheritParams pgpd
#' @inheritParams pextgp
#' @return vector of quantiles
#' @export
#' @keywords internal
qextgp <- function(p,
                   scale = 1,
                   shape1,
                   shape2,
                   lower.tail = TRUE){
  if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 1)
    stop("`p' must contain probabilities in (0,1)")
  if (length(scale) != 1 || min(scale) <= 0){
    stop("invalid scale")
  }
  if (length(shape1) != 1 || shape1 < 0){
    stop("invalid shape1")
  }
  if (length(shape2) != 1 || shape2 < -1){
    stop("invalid shape2")
  }
  if(isTRUE(all.equal(shape1, 0, check.attributes = FALSE))){
    # Exponential data
    return(qgpd(p = p,
                scale = scale,
                shape = shape2,
                lower.tail = lower.tail)
    )
  } else if(isTRUE(all.equal(shape2, 0, check.attributes = FALSE))){
    return(qgomp(p = p,
                 scale = scale,
                 shape = shape1,
                 lower.tail = lower.tail)
    )
  }
  if (lower.tail){
    p <- 1 - p
  }
  return(scale / shape1 * log(shape2 / shape1 * ((1 - p)^(-shape2) - 1) + 1))
}



#' Excess lifetime distributions
#'
#' Quantile and distribution function of excess lifetime distribution
#' for threshold exceedances.
#' @param q vector of quantiles.
#' @param p vector of probabilities
#' @param scale scale parameter(s)
#' @param shape vector of shape parameter(s).
#' @param lower.tail logical; if \code{TRUE} (default), the lower tail probability \eqn{\Pr(X \leq x)} is returned.
#' @param log.p logical; if \code{FALSE} (default), values are returned on the probability scale.
#' @param family string indicating the parametric model, one of \code{exp}, \code{gp}, \code{gomp}, \code{gompmake}, \code{weibull} and \code{extgp}
#' @name elife

#' @rdname elife
#' @export
#' @keywords internal
qelife <- function(p,
                   scale,
                   shape,
                   family = c("exp", "gp", "weibull", "gomp", "gompmake", "extgp"),
                   lower.tail = TRUE){
  family <- match.arg(family)
  switch(family,
         exp = qexp(p, rate = 1/scale, lower.tail = lower.tail),
         gp = qgpd(p = p, loc = 0, scale = scale, shape = shape, lower.tail = lower.tail),
         weibull = qweibull(p = p, shape = shape, scale = scale, lower.tail = lower.tail),
         gomp = qgomp(p = p, scale = scale, shape = shape, lower.tail = lower.tail),
         gompmake = qgompmake(p = p, scale = scale[1], lambda = scale[2], shape = shape, lower.tail = lower.tail),
         extgp = qextgp(p = p, scale = scale, shape1 = shape[1], shape2 = shape[2], lower.tail = lower.tail))
}

#' @rdname elife
#' @export
#' @keywords internal
pelife <- function(q,
                   scale,
                   shape,
                   family = c("exp", "gp", "weibull", "gomp", "gompmake","extgp"),
                   lower.tail = TRUE,
                   log.p = FALSE){
  family <- match.arg(family)
  switch(family,
         exp = pexp(q = q, rate = 1/scale, lower.tail = lower.tail, log.p = log.p),
         gp = pgpd(q = q, loc = 0, scale = scale, shape = shape, lower.tail = lower.tail, log.p = log.p),
         weibull = pweibull(q = q, shape = shape, scale = scale, lower.tail = lower.tail, log.p = log.p),
         gomp = pgomp(q = q, scale = scale, shape = shape, lower.tail = lower.tail, log.p = log.p),
         gompmake = pgompmake(q = q, scale = scale[1], shape = shape, lambda = scale[2], lower.tail = lower.tail, log.p = log.p),
         extgp = pextgp(q = q, scale = scale, shape1 = shape[1], shape2 = shape[2], lower.tail = lower.tail, log.p = log.p))
}

#' @rdname elife
#' @export
#' @keywords internal
relife <- function(n,
                   scale,
                   shape,
                   family = c("exp", "gp", "weibull", "gomp", "gompmake", "extgp")
){
  family <- match.arg(family)
  switch(family,
         exp = rexp(n = n, rate = 1/scale),
         gp = qgpd(p = runif(n), loc = 0, scale = scale, shape = shape),
         weibull = rweibull(n = n, shape = shape, scale = scale),
         gomp = qgomp(p = runif(n), scale = scale, shape = shape),
         gompmake = qgompmake(p = runif(n), scale = scale[1], lambda = scale[2], shape = shape),
         extgp = qextgp(p = runif(n), scale = scale, shape1 = shape[1], shape2 = shape[2])
  )
}
