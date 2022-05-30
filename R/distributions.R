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
  stopifnot("\"loc\" must be a vector of length 1." = length(loc) == 1L,
            "\"loc\" must be finite" = isTRUE(is.finite(loc)))
  check_elife_dist(scale = scale,
                   shape = shape,
                   family = "gp")
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
  stopifnot("\"loc\" must be a vector of length 1." = length(loc) == 1L,
            "\"loc\" must be finite" = isTRUE(is.finite(loc)))
  check_elife_dist(scale = scale,
                   shape = shape,
                   family = "gp")
  if(isTRUE(all.equal(shape, 0, check.attributes = FALSE))){
    return(dexp(x = x - loc, rate = 1/scale, log = log))
  }
  d <- (x - loc)/scale
  index <- (d >= 0 & ((1 + shape * d) >= 0)) | is.na(d)
  d[index] <- log(1/scale) -
    (1/shape + 1) * log(1 + shape * d[index])
  d[!index & !is.na(d)] <- -Inf
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
    stop("\"p' must contain probabilities in (0,1)")
  stopifnot("\"loc\" must be a vector of length 1." = length(loc) == 1L,
            "\"loc\" must be finite" = isTRUE(is.finite(loc)))
  check_elife_dist(scale = scale,
                   shape = shape,
                   family = "gp")
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
  check_elife_dist(scale = scale,
                   shape = shape,
                   family = "gomp")
  if(isTRUE(all.equal(shape, 0, check.attributes = FALSE))){
    # Exponential data
    return(pexp(q = q,
                rate = 1/scale,
                lower.tail = lower.tail,
                log.p = log.p)
    )
  } else {
    p <- pmax(0, 1 - exp(-expm1(exp(log(shape) + log(q) - log(scale)))/shape))
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
  if (min(p, na.rm = TRUE) < 0 || max(p, na.rm = TRUE) > 1){
    stop("\"p' must contain probabilities in (0,1)")
  }
  check_elife_dist(scale = scale,
                   shape = shape,
                   family = "gomp")
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
                      shape = 0,
                      lambda = 0,
                      lower.tail = TRUE,
                      log.p = FALSE){
  check_elife_dist(scale = c(scale, lambda),
                   shape = shape,
                   family = "gompmake")
  if(isTRUE(all.equal(shape, 0, check.attributes = FALSE))){
    # Exponential data
    return(pexp(q = q,
                rate = lambda + 1/scale,
                lower.tail = lower.tail,
                log.p = log.p)
    )
  } else {
    p <- pmax(0, 1 - exp(-lambda*q - (exp(shape*q/scale)-1)/shape))
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
#' @note The quantile function is defined in terms of Lambert's W function. Particular parameter combinations (small values of \code{lambda} lead to numerical overflow; the function throws a warning when this happens.
#'
#' @param p vector of probabilities.
#' @inheritParams pgpd
#' @param lambda exponential rate
#' @export
#' @return vector of quantiles
#' @keywords internal
qgompmake <- function(p,
                      scale = 1,
                      shape = 0,
                      lambda = 0,
                      lower.tail = TRUE){
  check_elife_dist(scale = c(scale, lambda),
                   shape = shape,
                   family = "gompmake")
      # stopifnot("Install package \"gsl\" to use \"qgompmake\" with the Gompertz model.\n Try \"install.packages(\"gsl\")\"" = requireNamespace("gsl", quietly = TRUE))
      if (min(p, na.rm = TRUE) < 0 || max(p, na.rm = TRUE) > 1){
        stop("\"p' must contain probabilities in (0,1)")
      }
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
  q <- rep(NA, length(p))
  y <- exp(1/(scale*lambda))*(1-p)^(-shape/(scale*lambda)) / (scale*lambda)
 index <- is.finite(y)
 if(sum(index) != sum(is.finite(p) & p < 1)){
   warning("Numerical overflow: some quantiles map to infinity.")
 }
  q[index] <-
       1/(shape*lambda) - log(1 - p[index]) / lambda - scale / shape * .LambertW0(y[index])
  q[is.infinite(y)] <- Inf
  return(q)
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
                      shape = 0,
                      lambda = 0,
                      log = FALSE){
  check_elife_dist(scale = c(scale, lambda),
                   shape = shape,
                   family = "gompmake")
  if(isTRUE(all.equal(shape, 0, check.attributes = FALSE))){
    # Exponential data - but parameter not identifiable
    return(dexp(x = x,
                rate = lambda + 1/scale,
                log = log)
    )
  }
  ld1 <- log(lambda + exp(shape * x / scale) / scale)
  ld1 <- ifelse(is.finite(ld1), ld1, 0)
  ldens <- ld1 - lambda*x - (exp(shape * x / scale) - 1) / shape
  ldens <- ifelse(x < 0 | is.infinite(x), -Inf, ldens)
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
                   shape1 = 0,
                   shape2 = 0,
                   lower.tail = TRUE,
                   log.p = FALSE){
  check_elife_dist(scale = scale,
                   shape = c(shape1, shape2),
                   family = "extgp")
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
  p <-  1 - pmax(0,(1+shape2*(exp(shape1*pmax(0,q)/scale)-1)/shape1))^(-1/shape2)
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
#' @export
#' @keywords internal
dgomp <- function(x,
                  scale = 1,
                  shape = 0,
                  log = FALSE){
  check_elife_dist(scale = scale,
                   shape = shape,
                   family = "gomp")
  if(shape < 1e-8){
    return(dexp(x = x,
                rate = 1/scale,
                log = log))
  } else{
    ldens <-  ifelse(x < 0 |  is.infinite(x), -Inf, -log(scale) + (shape*x/scale - exp(shape*x/scale)/shape + 1/shape))
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
                   shape1 = 0,
                   shape2 = 0,
                   log = FALSE){
  check_elife_dist(scale = scale,
                   shape = c(shape1, shape2),
                   family = "extgp")
  if(abs(shape2) < 1e-8 && abs(shape1) < 1e-8){
    return(dexp(x = x, rate = 1/scale, log = log))
  } else if(abs(shape2) < 1e-8 && abs(shape1) > 1e-8){ #Gompertz
    ldens <-  -log(scale) + (shape1*x/scale - exp(shape1*x/scale)/shape1 + 1/shape1)
    ldens <- ifelse(x < 0, -Inf, ldens)
  } else if(abs(shape2) >= 1e-8 && abs(shape1) < 1e-8){ #generalized Pareto
   return(dgpd(x = x, loc = 0, scale = scale,
               shape = shape2, log = log))
  } else { #extended
    ldens <- suppressWarnings(-log(scale) + (-1/shape2 - 1)*log(shape2*(exp(shape1*x/scale) - 1)/shape1 + 1) + shape1*x/scale)

    ldens <- ifelse(x < 0 | is.infinite(x), -Inf, ldens)
    if(shape2 < 0){
      ldens <- ifelse(x > scale * log(1-shape1/shape2) / shape1,
                      -Inf,
                      ldens)
    }
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
                   shape1 = 0,
                   shape2 = 0,
                   lower.tail = TRUE){
  if (min(p, na.rm = TRUE) < 0 || max(p, na.rm = TRUE) > 1){
    stop("\"p\" must contain probabilities in [0,1]")
  }
  check_elife_dist(scale = scale,
                   shape = c(shape1, shape2),
                   family = "extgp")
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
  return(scale / shape1 * log(shape1 / shape2 * (p^(-shape2) - 1) + 1))
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
                   scale = 1,
                   shape,
                   family = c("exp", "gp", "weibull", "gomp", "gompmake", "extgp"),
                   lower.tail = TRUE){
  family <- match.arg(family)
  if(missing(shape) & family != "exp"){
    stop("Missing \"shape\" parameter.")
  }
  check_elife_dist(scale = scale,
                   shape = shape,
                   family = family)
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
                   scale = 1,
                   shape,
                   family = c("exp", "gp", "weibull", "gomp", "gompmake","extgp"),
                   lower.tail = TRUE,
                   log.p = FALSE){
  family <- match.arg(family)
  if(missing(shape) & family != "exp"){
    stop("Missing \"shape\" parameter.")
  }
  check_elife_dist(scale = scale,
                   shape = shape,
                   family = family)
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
                   scale = 1,
                   shape,
                   family = c("exp", "gp", "weibull", "gomp", "gompmake", "extgp")
){
  family <- match.arg(family)
  if(missing(shape) & family != "exp"){
    stop("Missing \"shape\" parameter.")
  }
  check_elife_dist(scale = scale,
                   shape = shape,
                   family = family)
  switch(family,
         exp = rexp(n = n, rate = 1/scale),
         gp = qgpd(p = runif(n), loc = 0, scale = scale, shape = shape),
         weibull = rweibull(n = n, shape = shape, scale = scale),
         gomp = qgomp(p = runif(n), scale = scale, shape = shape),
         gompmake = qgompmake(p = runif(n), scale = scale[1], lambda = scale[2], shape = shape),
         extgp = qextgp(p = runif(n), scale = scale, shape1 = shape[1], shape2 = shape[2])
  )
}


#' @rdname elife
#' @export
#' @keywords internal
delife <- function(x,
                   scale = 1,
                   shape,
                   family = c("exp", "gp", "weibull", "gomp", "gompmake", "extgp"),
                   log = FALSE
){

  family <- match.arg(family)
  if(missing(shape) & family != "exp"){
    stop("Missing \"shape\" parameter.")
  }
  check_elife_dist(scale = scale,
                   shape = shape,
                   family = family)
  switch(family,
         exp = dexp(x = x, rate = 1/scale, log = log),
         gp = dgpd(x = x, loc = 0, scale = scale, shape = shape, log = log),
         weibull = dweibull(x = x, shape = shape, scale = scale, log = log),
         gomp = dgomp(x = x, scale = scale, shape = shape, log = log),
         gompmake = dgompmake(x = x, scale = scale[1], lambda = scale[2], shape = shape, log = log),
         extgp = dextgp(x = x, scale = scale, shape1 = shape[1], shape2 = shape[2], log = log)
  )
}

#' @keywords internal
check_elife_dist <- function(scale,
                             shape,
                             family = c("exp", "gp", "weibull", "gomp", "gompmake", "extgp")){
  family <- match.arg(family)
  if(family != "gompmake"){
  stopifnot("Invalid scale parameter: must be a positive scalar." =
              isTRUE(all(length(scale) == 1L,
                         is.finite(scale),
                         scale > 0)))
  }
  if(family == "gp"){
    stopifnot("\"shape\" should be a vector of length 1." = length(shape) == 1L,
              "\"shape\" must be finite." = isTRUE(is.finite(shape)))
  } else if(family == "weibull"){
    stopifnot("\"shape\" should be a vector of length 1." = length(shape) == 1L,
              "\"shape\" must be positive." = isTRUE(is.finite(shape) & shape > 0)
    )
  } else if(family == "gomp"){
    stopifnot("\"shape\" should be a vector of length 1." = length(shape) == 1L,
    "\"shape\" must be non-negative." = isTRUE(is.finite(shape) & shape >= 0)
    )
  } else if(family == "gompmake"){
    stopifnot("Invalid scale parameter." =
                isTRUE(all(length(scale) == 2L,
                           is.finite(scale[1]),
                           scale[1] > 0)),
              "\"shape\" should be a vector of length 1." = length(shape) == 1L,
               "\"shape\" must be non-negative." = isTRUE(all(is.finite(shape), shape >= 0)),
                                                          "\"lambda\" must be non-negative." = isTRUE(all(is.finite(scale[2]), scale[2] >= 0))
    )
  } else if(family == "extgp"){
    stopifnot("\"shape\" should be a vector of length 2." = length(shape) == 2L,
              "\"shape1\" must be non-negative." = isTRUE(is.finite(shape[1]) & shape[1] >= 0),
    "\"shape2\" must be finite." = isTRUE(is.finite(shape[2]))
    )
  }
}


# The Lambert W function
#
# This is extracted from VGAM
# @author Thomas W. Yee
# @keywords internal
.LambertW0 <- function (x, tolerance = 1e-10, maxit = 50)
{
  ans <- x
  ans[!is.na(x) & x < -exp(-1)] <- NA
  ans[!is.na(x) & x >= -exp(-1)] <- log1p(x[!is.na(x) & x >=
                                              -exp(-1)])
  ans[!is.na(x) & x >= 0] <- sqrt(x[!is.na(x) & x >= 0])/2
  cutpt <- 3
  if (any(myTF <- !is.na(x) & x > cutpt)) {
    L1 <- log(x[!is.na(x) & x > cutpt])
    L2 <- log(L1)
    wzinit <- L1 - L2 + (L2 + (L2 * (-2 + L2)/(2) + (L2 *
                                                       (6 + L2 * (-9 + L2 * 2))/(6) + L2 * (-12 + L2 * (36 +
                                                                                                          L2 * (-22 + L2 * 3)))/(12 * L1))/L1)/L1)/L1
    ans[myTF] <- wzinit
  }
  for (ii in 1:maxit) {
    exp1 <- exp(ans)
    exp2 <- ans * exp1
    delta <- (exp2 - x)/(exp2 + exp1 - ((ans + 2) * (exp2 -
                                                       x)/(2 * (ans + 1))))
    ans <- ans - delta
    if (all(is.na(delta)) || max(abs(delta), na.rm = TRUE) <
        tolerance)
      break
    if (ii == maxit)
      warning("did not converge")
  }
  ans[is.infinite(NA)] <- Inf
  ans
}
