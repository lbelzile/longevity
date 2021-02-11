#' @keywords internal
hazard_fn_elife <- function(x,
                         par,
                         family = c("exp","gp","gomp","weibull","extgp")){
  family <- match.arg(family)
  stopifnot("`x` must be numeric" = is.numeric(x),
            "`x` must be positive" = isTRUE(all(x > 0)),
    "Length of parameter vector does not match model definition" =length(par) == switch(family, exp = 1, gp = 2, gomp = 2, weibull = 2, extgp = 3)
    )
  if(family == 'extgp' && isTRUE(all.equal(par[3], 0, check.attributes = FALSE))){
    family <- "gomp"
  }
  if(family == 'extgp' && isTRUE(all.equal(par[2], 0, check.attributes = FALSE))){
    family <- "gp"
    par <- par[-2]
  }
  # Define hazard functions for various parametric models
switch(family,
       exp = rep(1/par[1], length.out = length(x)),
       gp =  1/(par[1] + par[2] * x),
       weibull = par[2]*par[1]^(-par[2])*x^(par[2]-1), #par[1] = scale, par[2] = shape
       extgp =  (1/par[1]) * exp(par[2]*x/par[1])/(1+par[3]*(exp(par[2]*x/par[1])-1)/par[2]),
       gomp = exp(par[2]*x/par[1])/par[1]
       )
}

#' Profile likelihood for hazard
#'
#' This function computes the hazard for different \code{elife} parametric
#' model with profile-likelihood based confidence intervals.
#' It is also used to provide local hazard plots at varying thresholds.
#'
#' @param x value of the threshold exceedance at which to estimate the hazard
#' @param plot logical; if true, display the profile log-likelihood. Default to \code{FALSE}.
#' @inheritParams nll_elife
#' @return an invisible object of class \code{elife_hazard} containing information about the profile likelihood
#' @examples
#' n <- 2500
#' dat <- rdtrunc_elife(n = n, scale = 2,
#' family = "gp", shape = 0.1,
#' lower = ltrunc <- runif(n),
#' upper = rtrunc <- (5 + runif(n)))
#' hazard_elife(x = 2, dat = dat, ltrunc = ltrunc, rtrunc = rtrunc, type = "ltrt", family = "gp")
hazard_elife <- function(x,
                        dat,
                        thresh = 0,
                        type = c("none","ltrt","ltrc"),
                        rcens,
                        ltrunc,
                        rtrunc,
                        family = c("exp","gp","gomp","weibull","extgp"),
                        weights = rep(1, length(dat)),
                        level = 0.95,
                        psi = NULL,
                        plot = FALSE){
  stopifnot("The level of the confidence interval must be in [0,1]" = level > 0 && level < 1,
            "Level should be a numeric of length one" = length(level) == 1L,
            "The value of `x` should be numeric." = is.numeric(x)
            )
   # This function should
  # 1) create a wrapper function to reparametrize the likelihood in terms of hazard
  # 2) compute the restricted maximum likelihood at each point
  # 3) return an object with hazard and the log-likelihood value
  # 4) include a @confint and a @plot method to display the profile
  # Compute maximum likelihood estimator
  mle <- optim_elife(dat = dat,
                     thresh = thresh,
                     ltrunc = ltrunc,
                     rtrunc = rtrunc,
                     rcens = rcens,
                     type = type,
                     family = family,
                     weights = weights)
  # Compute MLE of hazard
  mle_haz <- as.numeric(hazard_fn_elife(par = mle$par, x = x, family = family))

  if(family == "exp"){
    # constant hazard function
    if(is.null(psi)){
      psi <- 1/(mle$par + seq(to = -mle$std.error*3.5, from = mle$std.error*7, length.out = 101))
    }
    psi <- psi[psi>0]
    npll <- sapply(psi, function(par){
      nll_elife(par = 1/par,
                dat = dat,
                thresh = thresh,
                type = type,
                rcens = rcens,
                ltrunc = ltrunc,
                rtrunc = rtrunc,
                family = family,
                weights = weights)})
    # plot(psi, shifted_npll, type = "l");
  } else if(family == "gp"){
    if(mle$par[2] < 0 && x > -mle$par[1]/mle$par[2]){
      stop("Value of x is outside of the range of the distribution evaluated at the maximum likelihood estimate.")
    }
    inv_haz <- function(hazard, xi, x){
        as.numeric((1/hazard)-xi*x)
    }
    if(is.null(psi)){
      haz_stderror <- sqrt(diag(solve(numDeriv::hessian(func = function(par){
        nll_elife(par = c(inv_haz(par, xi = mle$par[2], x = x), mle$par[1]), dat = dat, thresh = thresh, type = type, rcens = rcens, ltrunc = ltrunc, rtrunc = rtrunc, family = family, weights = weights)},
        x = mle_haz))))
      psi <- mle_haz + seq(-4*haz_stderror, 4*haz_stderror, length.out = 101L)
    }
      psi <- psi[psi > 0]
    mdat <- max(dat) - thresh
    ubound <- ifelse(x > mdat, x, mdat)
    npll <- sapply(psi, function(haz_i){
      opt <- optimize(f = function(xi){
          nll_elife(par = c(1/haz_i-xi*x, xi),
                    dat = dat,
                    rcens = rcens,
                    ltrunc = ltrunc,
                    rtrunc = rtrunc,
                    weights = weights,
                    family = family,
                    type = type,
                    thresh = thresh)},
        upper = 1/(haz_i*x),
        lower = -1)
      opt$objective
    })


  } else if(family == "weibull"){

    inv_haz <- function(hazard, alpha, x){
      as.numeric((hazard*x^(1-alpha)/alpha)^(-1/alpha))
    }
    if(is.null(psi)){
    haz_stderror <- sqrt(diag(solve(numDeriv::hessian(func = function(par){
      nll_elife(par = c(inv_haz(par[1], mle$par[2], x = x), mle$par[2]), dat = dat, thresh = thresh, type = type, rcens = rcens, ltrunc = ltrunc, rtrunc = rtrunc, family = family, weights = weights)},
      x = mle_haz)))[1])
    psi <- mle_haz + seq(-5*haz_stderror, 5*haz_stderror, length.out = 101L)
    }
    psi <- psi[psi > 0]
    npll <- matrix(0, nrow = 2, ncol = length(psi))
    mid <- which.min(abs(psi-mle_haz))
    for(i in c(mid:length(psi), (mid-1):1)){
      opt <- Rsolnp::solnp(pars = ifelse(i >= mid, ifelse(i==mid, mle$par[2], npll[1,i-1]), npll[1,i+1]),
                           f = function(alpha){
        nll_elife(par = c(inv_haz(hazard = psi[i], alpha = alpha, x = x), alpha),
                  dat = dat,
                  rcens = rcens,
                  ltrunc = ltrunc,
                  rtrunc = rtrunc,
                  weights = weights,
                  family = family,
                  type = type,
                  thresh = thresh)},
        ineqfun = function(alpha){alpha},
        ineqLB = 0,
        ineqUB = Inf, control = list(trace = 0))
      npll[,i] <- c(opt$pars, opt$values[length(opt$values)])
    }
    npll <- npll[2,]
  } else if(family == "gomp"){

    inv_haz <- function(hazard, sigma, x){
      as.numeric(sigma*log(sigma*hazard)/x)
    }
    if(is.null(psi)){
    # if(mle$par[2] < 1e-4){
    #   psi <- seq(from = 0.01,
    #              to = 1/mle$par[1]+4/mle$par[1],
    #              length.out = 101)
    # } else{
      haz_stderror <- sqrt(diag(solve(numDeriv::hessian(func = function(par){
        nll_elife(par = c(mle$par[1],inv_haz(par, mle$par[1], x = x)),
                  dat = dat, thresh = thresh,
                  type = type, rcens = rcens,
                  ltrunc = ltrunc, rtrunc = rtrunc,
                  family = family, weights = weights)},
      x = mle_haz,
      method.args=list(eps=1e-7, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE))
      )))[1]
      psi <- mle_haz + seq(-5*haz_stderror, 5*haz_stderror, length.out = 101L)
    }
      psi <- psi[psi > 0]
    # }
    npll <- matrix(0, nrow = 2, ncol = length(psi))
    mid <- which.min(abs(psi-mle_haz)) + 1L
    for(i in c(mid:length(psi), (mid-1):1)){
      # Ensure a feasible solution
      opt <- Rsolnp::solnp(pars = pmax(1/psi[i]+1e-4, ifelse(i >= mid, ifelse(i==mid, mle$par[1], npll[1,i-1]), npll[1,i+1])),
                           f = function(sigma){
                             nll_elife(par = c(sigma, inv_haz(hazard = psi[i], sigma = sigma, x = x)),
                                       dat = dat,
                                       rcens = rcens,
                                       ltrunc = ltrunc,
                                       rtrunc = rtrunc,
                                       weights = weights,
                                       family = family,
                                       type = type,
                                       thresh = thresh)},
                           ineqfun = function(alpha){alpha},
                           ineqLB = 1/psi[i],
                           ineqUB = 1e8, control = list(trace = 0))
      npll[,i] <- c(opt$pars, opt$values[length(opt$values)])
    }
    npll <- npll[2,]

  } else if(family == "extgp"){

  inv_haz <- function(hazard, beta, sigma, x){
      if(abs(beta) > 1e-6){
        as.numeric(beta*(exp(beta*x/sigma)/(sigma*hazard)-1)/(exp(beta*x/sigma)-1))
      } else{
        as.numeric((1/hazard-sigma)/x)
      }
    }
    # Compute grid of values at which to evaluate the hazard
  if(is.null(psi)){
    haz_stderror <- sqrt(diag(solve(numDeriv::hessian(func = function(par){
      nll_elife(par = c(mle$par[1:2], inv_haz(par, sigma = mle$par[1], beta = mle$par[2], x = x)),
                dat = dat, thresh = thresh, type = type, rcens = rcens, ltrunc = ltrunc, rtrunc = rtrunc, family = family, weights = weights)},
      x = mle_haz))))
      psi <- seq(from = mle_haz - 5*haz_stderror,
                 to = mle_haz + 5*haz_stderror,
                 length.out = 101)
    }
    psi <- psi[psi > 0]

    npll <- matrix(0, nrow = 3, ncol = length(psi))
    mid <- which.min(abs(psi-mle_haz)) + 1L
    for(i in c(mid:length(psi), (mid-1):1)){
      if(i == mid){
        start <- mle$par[1:2]
      } else if(i > mid){
        start <- npll[1:2, i-1]
      } else{
        start <- npll[1:2, i+1]
      }
      # Ensure a feasible solution
      opt <- Rsolnp::solnp(pars = start,
                           f = function(par){
                             sigma = par[1];
                             beta = par[2];
                             nll_elife(par = c(sigma, beta, inv_haz(hazard = psi[i], sigma = sigma, beta = beta, x = x)),
                                       dat = dat,
                                       rcens = rcens,
                                       ltrunc = ltrunc,
                                       rtrunc = rtrunc,
                                       weights = weights,
                                       family = family,
                                       type = type,
                                       thresh = thresh)},
                           ineqfun = function(par){
                             sigma = par[1];
                             beta = par[2]
                             xi = inv_haz(hazard = psi[i], sigma = sigma, beta = beta, x = x)
                             c(sigma, beta, xi,
                                 ifelse(xi < 0, -beta/xi, 1e-5),
                                 ifelse(xi < 0, thresh + sigma/beta*log(1-beta/xi) - max(dat), 1e-5)
                               )
                             },
                             ineqLB = c(0, 0, -1, 0, 0),
                             ineqUB = c(rep(Inf, 2), 10, rep(Inf, 2)),
                           control = list(trace = 0))
      npll[,i] <- c(opt$pars, opt$values[length(opt$values)])
    }
    npll <- npll[3,]
  }

  profile_confint <- function(psi, npll, psi_mle, nll_mle, level = 0.95){
    stopifnot("Maximum likelihood estimate should be numeric." = is.numeric(psi_mle),
              "Length of mle must be one" = length(psi_mle) == 1L,
              "Grid of psi values must contain maximum likelihood estimate for the hazard." = min(psi) < psi_mle & max(psi) > psi_mle,
              "Log likelihood values must be the same length as points over which profiling is done." = length(psi) == length(npll),
              "The log likelihood function must be largest at the MLE." =  nll_mle <= min(npll))
    pred <- c(predict(smooth.spline(x = npll[psi < psi_mle] - nll_mle, y = psi[psi < psi_mle]), qchisq(level, 1)/2)$y,
              predict(smooth.spline(x = npll[psi > psi_mle] - nll_mle, y = psi[psi > psi_mle]), qchisq(level, 1)/2)$y)
    return(pred)
  }
  mnpll <- which.min(npll)
  if(npll[mnpll] < -mle$loglik){
    mle_haz <- psi[mnpll]
    nll_mle <- npll[mnpll]
  } else{
    nll_mle <- -mle$loglik
  }
  pconfint <- profile_confint(psi = psi,
                              psi_mle = mle_haz,
                              npll = npll,
                              nll_mle = nll_mle,
                              level = level)
  shifted_npll <- - (npll + mle$loglik)
  keep <- shifted_npll > -7
  retv <- structure(
    list(hazards = psi[keep],
         pll = shifted_npll[keep],
         confint = pconfint,
         par = as.numeric(mle_haz),
         level = level),
    class = "elife_hazard"
  )
  if(plot){
    plot(retv)
  }
  return(invisible(retv))
}

#' @export
plot.elife_hazard <- function(x, ...){
  args <- list(...)
  plot(x = x$hazards,
       y = x$pll,
       type = "l",
       bty = "l",
       xlab = "hazard",
       ylab = "profile log-likelihood",
       ylim = c(-5,0),
       panel.first = {
  abline(h = -qchisq(level, 1)/2, col = "gray")}, ...)
  rug(c(x$confint, x$par), ticksize = 0.05)
}

confint.elife_hazard <- function(x, ...){
x$confint

}
#
# hazard_plot_elife <- function(dat, thresh){
#     # TODO This function should take as argument the result
#     # of @elife_hazard and loop over thresh
#     # + return a
# }
#
# local_hazard_plot_elife <- function(dat, thresh){
#   # TODO This function should take as argument
#   #
#   # of @elife_hazard and loop over thresh
#   # + return a
# }
#

