#' Distribution function of the generalized Pareto distribution
#'
#' @param q vector of quantiles.
#' @param loc location parameter.
#' @param scale scale parameter, strictly positive.
#' @param shape shape parameter.
#' @param lower.tail logical; if \code{TRUE} (default), the lower tail probability \eqn{\Pr(X \leq x)} is returned.
#' @param log.p logical; if \code{FALSE} (default), values are returned on the probability scale.
#' @return a vector of (log)-probabilities of the same length as \code{q}
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
  nn <- x
  index <- (d > 0 & ((1 + shape * d) > 0)) | is.na(d)
  if (shape == 0) {
    d[index] <- log(1/scale[index]) - d[index]
    d[!index] <- -Inf
  }
  else {
    d[index] <- log(1/scale[index]) - (1/shape + 1) * log1p(shape * d[index])
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
#' @keywords internal
qgpd <- function(p,
                 loc = 0,
                 scale = 1,
                 shape = 0,
                 lower.tail = TRUE){
  if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 1)
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


#' Distribution function of the Gompertz-Makeham distribution
#'
#' @inheritParams pgpd
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

#' Quantile function of the Gompertz-Makeham distribution
#'
#' @param p vector of probabilities.
#' @inheritParams pgpd
#' @return vector of quantiles
#' @keywords internal
qgomp <- function(p,
                  scale = 1,
                  shape = 0,
                  lower.tail = TRUE){
  if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 1)
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
  if (lower.tail){
    p <- 1 - p
  }
  return(scale / shape * log(1 - shape * log(1 - p)))
}

#' Distribution function of the extended generalized Pareto distribution
#'
#' @inheritParams pgpd
#' @param shape1 positive shape parameter \eqn{\beta}; model defaults to generalized Pareto when it equals zero.
#' @param shape2 shape parameter \eqn{\gamma}; model reduces to Gompertz when \code{shape2=0}.
#' @return a vector of (log)-probabilities of the same length as \code{q}
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
  p <- 1 - pmax(0,(1+shape[2]*(exp(shape1*q/scale)-1)/shape1)^(-1/shape2))
  if (!lower.tail){
    p <- 1 - p
  }
  if (log.p){
    p <- log(p)
  }
  return(p)
}

#' Quantile function of the extended generalized Pareto distribution
#'
#' @param p vector of probabilities.
#' @inheritParams pgpd
#' @inheritParams pextgp
#' @return vector of quantiles
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
  if (length(shape2) != 1){
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

#' Likelihood for left-truncated and right-censored or interval-truncated data
#'
#' Computes the log-likelihood for various parametric models suitable for threshold
#' exceedances
#' @param par vector of parameters
#' @param thresh threshold
#' @param type string, either \code{ltrt} for left- and right-truncated data or \code{ltrc} for left-truncated right-censored data
#' @param dat vector of threshold exceedances
#' @param rcens logical indicating right-censoring (\code{TRUE} for censored)
#' @param ltrunc lower truncation limit, possibly zero
#' @param ltrunc upper truncation limit
#' @param weights weights for observations
#' @param family string; choice of parametric family, either exponential (\code{exp}), Weibull (\code{weibull}), generalized Pareto (\code{gp}), Gompertz (\code{gomp}) or extended generalized Pareto (\code{extgp}).
#' @return log-likelihood value
nll_elife <- function(par,
                      dat,
                      thresh = 0,
                      type = c("none","ltrt","ltrc"),
                      rcens,
                      ltrunc,
                      rtrunc,
                      family = c("exp","gp","gomp","weibull","extgp"),
                      weights = rep(1, length(dat))){
  family <- match.arg(family)
  type <- match.arg(type)
  stopifnot("Threshold must be non-negative" = thresh >= 0,
            "Only a single threshold value is allowed" = length(thresh) == 1L
            )
  if(thresh > 0){
    # Keep only exceedances, shift observations
    ind <- dat > thresh
    dat <- dat[ind] - thresh
    if(type != "none"){
      ltrunc <- pmax(0, ltrunc[ind] - thresh)
    }
    weights <- weights[ind]
    if(type == "ltrc"){
      if(missing(rcens) || missing(ltrunc)){
        stop("Missing inputs for left-truncated and right-censored data (`rcens` or `ltrunc`).")
      }
      rcens <- rcens[ind]
    } else if(type == "ltrt"){
      if(missing(rtrunc) || missing(ltrunc)){
        stop("Missing inputs for left- and right-truncated data (`ltrunc` or `rtrunc`).")
      }
      rtrunc <- rtrunc[ind] - thresh
    }
  }
  if(family == "gomp"){
    family = "extgp"
    par <- c(par, 0)
  }
  # Define density and survival function for each of the parametric families
  if(family == "exp"){
    stopifnot("Length of parameters for exponential family is one" = length(par) == 1)
    if(par < 0){
      return(1e10)
    }
    ldensf <- function(par, dat){ dexp(x = dat, rate = 1/par[1], log = TRUE)}
    lsurvf <- function(par, dat){ pexp(q = dat, rate = 1/par[1], lower.tail = FALSE, log.p = TRUE)}
  } else if(family == "gp"){
    stopifnot("Length of parameters for generalized Pareto family is two." = length(par) == 2)
    if(par[1] < 0 || par[2] < (-1+1e-8)){
      return(1e10)
    }
    if(par[2] < 0 && (-par[2]/par[1] < max(dat))){
         return(1e10)
    }
    ldensf <- function(par, dat){ dgpd(x = dat, loc = 0, scale = par[1], shape = par[2], log = TRUE)}
    lsurvf <- function(par, dat){ pgpd(q = dat, loc = 0, scale = par[1], shape = par[2], lower.tail = FALSE, log.p = TRUE)}
  } else if(family == "weibull"){
    if(par[1] <= 0 || par[2] <= 0){
      return(1e10)
    }
      ldensf <- function(par, dat){dweibull(x = dat, scale = par[1], shape = par[2], log = TRUE)}
      lsurvf <- function(par, dat){pweibull(q = dat, scale = par[1], shape = par[2], lower.tail = FALSE, log.p = TRUE)}
  } else if(family == "extgp"){
    # scale par[1]
    # beta = par[2]; if zero, recover generalized Pareto
    # xi = par[3]; if zero, recover Gompertz
    maxdat <- max(dat)
    # both shape parameters can be negative, but the data
    # must lie inside the support of the distribution
    if(abs(par[3]) <= 1e-8 && abs(par[2]) <= 1e-8){
      bounds <- par[1]
    } else if(abs(par[3]) > 1e-8 && abs(par[2]) < 1e-8){
      bounds <- c(par[1], ifelse(par[3] < 0, 1+par[3]*maxdat/par[1], Inf))
    } else if(abs(par[3]) < 1e-8 && abs(par[2]) > 1e-8){
      bounds <- par[1:2]
    } else {
      bounds <- c(par[1:2], 1+par[3]*(exp(par[2]*maxdat/par[1])-1)/par[2])
    }
    if(isTRUE(any(bounds <= 0))){
      return(1e10)
    }
    lsurvf <- function(dat, par){
      if(abs(par[3]) < 1e-8 && abs(par[2]) < 1e-8){ #exponential
        -dat/par[1]
      } else if(abs(par[3]) < 1e-8 && abs(par[2]) > 1e-8){ #Gompertz
        -exp(par[2]*dat/par[1])/par[2] + 1/par[2]
      } else if(abs(par[3]) >= 1e-8 && abs(par[2]) < 1e-8){ #generalized Pareto
        (-1/par[3])*log(1+dat*par[3]/par[1])
      } else if(abs(par[3]) >= 1e-8 && abs(par[2]) >= 1e-8){ #extended
        (-1/par[3])*log(1+par[3]*(exp(par[2]*dat/par[1])-1)/par[2])
      }
    }
    ldensf <- function(dat, par){
      if(abs(par[3]) < 1e-8 && abs(par[2]) < 1e-8){
        -log(par[1]) + -dat/par[1]
      } else if(abs(par[3]) < 1e-8 && abs(par[2]) > 1e-8){ #Gompertz
        -log(par[1]) + (par[2]*dat/par[1] - exp(par[2]*dat/par[1])/par[2] + 1/par[2])
      } else if(abs(par[3]) >= 1e-8 && abs(par[2]) < 1e-8){ #generalized Pareto
        - log(par[1]) - (1/par[3] + 1)*log(pmax(0, 1+dat*par[3]/par[1]))
      } else if(abs(par[3]) >= 1e-8 && abs(par[2]) >= 1e-8){ #extended
       -log(par[1]) + (-1/par[3] - 1)*log(par[3]*(exp(par[2]*dat/par[1]) - 1)/par[2] + 1) + par[2]*dat/par[1]
      }
    }
  }
  # Likelihoods, depending on the sampling scheme (left- and right-truncated, left-truncated and right-censored, none)
  if(type == "ltrc"){
    # Check number of parameters
    g1 <- intersect(which(!rcens), which(ltrunc > 0))
    g2 <- intersect(which(rcens), which(ltrunc > 0))
    ll <- 0
    #Contribution from observations in sampling frame (density divided by truncation interval)
    if(sum(!rcens)>0){
      ll <- sum(weights[!rcens]*ldensf(dat = dat[!rcens], par = par))
      if(length(g1) > 0){
        ll <- ll - sum(weights[g1]*lsurvf(dat = ltrunc[g1], par = par))  #right truncated individuals
      }
    }
    if(sum(rcens)>0){
      ll <- ll +  sum(weights[rcens]*lsurvf(dat  =dat[rcens], par = par))
      if(length(g2) > 0){
        ll <- ll - sum(weights[g2]*lsurvf(dat = ltrunc[g2], par = par))
      }
    }
  } else if(type == "ltrt"){
    ll <- sum(weights*(ldensf(dat = dat, par = par) -
                         log(exp(lsurvf(dat = ltrunc, par = par)) - exp(lsurvf(dat = rtrunc, par = par)))))
  } else if(type == "none"){
    ll <- sum(weights*ldensf(dat = dat, par = par))
  }
  if (!is.finite(ll)) {
    return(1e10)
  }  else {
    return(-ll)
  }
}

#' Optimization routine for excess lifetime models
#'
#' This function is a wrapper around constrained optimization
#' routines for different models with interval truncation and
#' left-truncation and right-censoring.
#'
#' @importFrom Rsolnp solnp
optim_elife <- function(dat,
                        thresh,
                        ltrunc,
                        rtrunc,
                        rcens,
                        type = c("none", "ltrc", "ltrt"),
                        family = c("exp", "gp", "weibull", "gomp", "extgp"),
                        weights){
  n <- length(dat)
  if(is.null(weights)){
    weights <- rep(1, length(dat))
  } else{
    stopifnot("All weights should be positive" = isTRUE(all(weights >=0 )))
  }
  family <- match.arg(family, several.ok = TRUE)
  type <- match.arg(type)
  if(type == "ltrc" && isTRUE(any(missing(ltrunc), missing(rcens)))){
    stop("User must provide `ltrunc` and `rcens` for left-truncated and right-censored data.")
  }
  if(type == "ltrt" && isTRUE(any(missing(ltrunc), missing(rtrunc)))){
    stop("User must provide `ltrunc` and `rtrunc` for left- and right-truncated data.")
  }
  stopifnot("Threshold must be positive" = thresh > 0,
            "Only a single threshold is allowed" = length(thresh) == 1L,
            "Only a single sampling scheme in `type` is allowed" = length(type) == 1L
  )
  if(!missing(ltrunc)){
    stopifnot("`dat` and `ltrunc` must be of the same length." = n == length(ltrunc),
              "`ltrunc` must be lower than `dat`" = isTRUE(all(ltrunc <= dat))
    )
  }
  if(!missing(rtrunc)){
    stopifnot("Argument `rtrunc` is not needed for chosen sampling scheme" = type == "ltrc",
              "`dat` and `rtrunc` must be of the same length." = n == length(rtrunc),
              "`rtrunc` must be higher than `dat`" = isTRUE(all(rtrunc >= dat)),
              "`ltrunc` must be larger than `ltrunc" = isTRUE(all(rtrunc > ltrunc))
    )
  }
  if(!missing(rcens)){
    stopifnot("Argument `rcens` is not needed for chosen sampling scheme" = type == "ltrc",
              "`dat` and `rcens` must be of the same length." = n == length(rcens),
              "`rcens` must be a vector of logicals.`" = islogical(rcens)
    )
  }
  # Perform constrained optimization for the different routines
  if(family == "exp"){
    # closed-form solution for left-truncated and right-censored data
    if(type == "none"){
      exc_i <- dat>thresh
      mle <- weighted.mean(x = (dat - thresh)[exc_i], w = weights[exc_i], na.rm = TRUE)
      se_mle <- mle / sqrt(sum(weights[exc_i]))
      ll <- sum(weights[exc_i]*dexp((dat - thresh)[exc_i], rate = 1/mle, log = TRUE))
      conv <- TRUE
    } else if(type == "ltrc"){
      exp_mle_ltrc <- function(dat, ltrunc, rcens, weights){
        sum(weights*(dat - ltrunc)) / sum(weights[!rcens])
      }
      exc_i <- dat > thresh
      mle <- exp_mle_ltrc(dat = (dat - thresh)[exc_i],
                          ltrunc = pmax(0, ltrunc[exc_i] - thresh),
                          rcens = rcens[exc_i],
                          weights = weights[exc_i])
      se_mle <- mle/sqrt(sum(weights[exc_i][!rcens[exc_i]]))
      ll <- - nll_elife(par = mle,
                        dat = dat,
                        ltrunc = ltrunc,
                        rcens = rcens,
                        type = "ltrc",
                        thresh = thresh,
                        family = "exp",
                        weights = weights)
      conv <- TRUE
    } else if(type == "ltrt"){
      opt_mle <- optim(par = mean(dat) - thresh,
                       fn = nll_elife,
                       method = "Brent",
                       lower = 1e-8,
                       upper = 3 * (max(dat) - thresh),
                       type = type,
                       thresh = thresh,
                       dat = dat,
                       ltrunc = ltrunc,
                       rtrunc = rtrunc,
                       weights = weights,
                       family = "exp",
                       hessian = TRUE
      )
      mle <- opt_mle$par
      se_mle <- sqrt(diag(solve(opt_mle$hessian)))
      ll <- -opt_mle$value
      conv <- opt_mle$convergence == 0
    }
  } else {
    if(family == "gp"){
      hin <- function(par, dat, thresh, ...){
        # scale > 0, xi > -1, xdat < -xi/sigma if xi < 0
        c(par[1], par[2], ifelse(par[2] < 0, thresh - par[1]/par[2] - max(dat)), 1e-5)
      }
      ineqLB <- c(0, -1, 0)
      ineqUB <- rep(Inf, 3)
      start <- mev::fit.gpd(xdat = dat, thresh = thresh)$par
    } else if(family == "weibull"){
      hin <- function(par, dat, thresh, ...){
        c(par[1], par[2])
      }
      ineqLB <- c(0,0)
      ineqUB <- rep(Inf, 2)
      start <- c(mean(dat), 1)
    } else if(family == "gomp"){
      hin <- function(par, dat, thresh, ...){
        par[1:2]
      }
      ineqLB <- c(0,0)
      ineqUB <- rep(Inf, 2)
      start <- c(mean(dat), 0.5)
    } else if(family == "extgp"){
      #parameters are (1) scale > 0, (2) beta >= 0 (3) gamma
      hin <- function(par, dat, thresh, ...){
        c(par[1:2],
          ifelse(par[3] < 0, 1-par[2]/par[3], 1e-5),
          ifelse(par[3] < 0, thresh + par[1]/par[2]*log(1-par[2]/par[3]) - max(dat), 1e-5)
        )
      }
      ineqLB <- rep(0, 4)
      ineqUB <- rep(Inf, 4)
      start <- c(mean(dat), 1, 0.1)
      #TODO try also fitting the GP/EXP/Gompertz and see which is best?
    }
    opt_mle <- Rsolnp::solnp(pars = start,
                             family = family,
                             ltrunc = ltrunc,
                             rtrunc = rtrunc,
                             rcens = rcens,
                             fun = nll_elife,
                             ineqfun = hin,
                             dat = dat,
                             thresh = thresh,
                             type = type,
                             weights = weights,
                             ineqLB = ineqLB,
                             ineqUB = ineqUB,
                             control = list(trace = 0))
    mle <- opt_mle$pars
    vcov <- try(solve(opt_mle$hessian[-(1:length(ineqLB)),-(1:length(ineqLB))]))
    if(is.character(vcov)){
      vcov <- NULL
    } else{
      se_mle <- try(sqrt(diag(vcov)))
    }
    if(is.character(se_mle)){
      se_mle <- rep(NA, length(mle))
    }
    ll <- -opt_mle$values[length(opt_mle$values)]
    conv <- opt_mle$convergence == 0
  }
  structure(list(par = mle,
                 std.error = se_mle,
                 loglik = ll,
                 nexc = sum(dat>thresh),
                 vcov = vcov,
                 convergence = conv),
            class = "longevity_parmod")
}

#' Maximum likelihood estimation of parametric models for excess lifetime
#'
#' This function performs constrained optimization to find the maximum likelihood estimator
#' for a range of parametric models suitable for the analysis of longevity data.
#' These can be estimated simultaneously for a range of thresholds.
#'
#' @inheritParams np.surv.interv
#' @param family string, one of \code{exp}, \code{gp}, \code{gomp} or \code{ext} for exponential, generalized Pareto, Gompertz or extended family, respectively.
#' @param thresh a vector of thresholds
#' @param weights vector of weights, defaults to one for each observation.
#' @return a list containing
fit_elife <- function(dat,
                      thresh,
                      ltrunc,
                      rtrunc,
                      rcens,
                      type = c("none", "ltrc", "ltrt"),
                      family = c("exp", "gp", "weibull", "gomp", "extgp"),
                      weights = NULL) {

# Return a list of tables with parameter estimates as
# a function of the different thresholds
results <- list()
results$par <- results$std.error <-
  matrix(NA,
        nrow = length(thresh),
        ncol = switch(family,
                      "exp" = 1L,
                      "gp" = 2L,
                      "gomp" = 2L,
                      "extgp" = 3L,
                      "weibull" = 2L)
          )
results$convergence <- rep(FALSE, length(thresh))
results$nexc <- rep(0, length(thresh))
results$thresh <- thresh
results$family <- family
results$type <- results$type
for(i in 1:length(thresh)){
  results_i <- optim_elife(dat = dat,
                           thresh = thresh[i],
                           ltrunc = ltrunc,
                           rtrunc = rtrunc,
                           rcens = rcens,
                           type = type,
                           family = family,
                           weights = weights
                          )
  results$nexc[i] <- sum(dat > thresh[i])
  results$par[i,] <- results_i$mle
  #TODO write tests for cases where standard errors do not exist
  results$std.error[i,] <- results_i$se
  results$convergence[i] <- results_i$convergence
  results$loglik[i] <- results_i$loglik
  }
  return(results)
}


#' @export
print.elife_par <-
  function(x,
          digits = getOption("digits"),
          na.print = "", ...){

      cat("Model:", x$family, "distribution.", "\n")
    if(x$type != "none"){
      cat("Sampling:",
          switch(x$type,
                 ltrt = "left- and right-truncated data",
                 ltrc = "left-truncated, right-censored data"),"\n")
    }
      cat("Log-likelihood:", round(x$loglik, digits), "\n")

      cat("\nThreshold:", round(x$thresh, digits), "\n")
      cat("Number of exceedances:", x$nexc, "\n")
      cat("\nEstimates\n")
      print.default(format(x$par, digits = digits), print.gap = 2, quote = FALSE, ...)
      if (!is.null(x$se) && x$estimate[1] > -0.5) {
        cat("\nStandard Errors\n")
        print.default(format(x$std.err, digits = digits), print.gap = 2, quote = FALSE, ...)
      }
      cat("\nOptimization Information\n")
      cat("  Convergence:", x$convergence, "\n")
  invisible(x)
}


#' @export
logLik.elife_par <- function(object, ...) {
  val <- object$loglik
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}



# Methods for class "mev_gpd", returned by mev::fit.gpd()

#' @export
nobs.elife_par <- function(object, ...) {
  return(object$nexc)
}

#' @export
coef.elife_par <- function(object, ...) {
  return(object$par)
}

#' @export
vcov.elife_par <- function(object, ...) {
  return(object$vcov)
}
