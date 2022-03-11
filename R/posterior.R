# Log of the maximal data information prior
lprior_mdi_elife <- function(par,
                           family = c("exp","gp","gomp")){
  family <- match.arg(family)
  obound <- switch(family,
                   exp = par[1] <= 0,
                   gp = par[1] <= 0 | par[2] < -1,
                   gomp = par[1] <=0 | par[2] <= 0)
  if(obound){
    return(-Inf)
  }
  if(family == "exp"){
    return(-log(par[1]))
  } else if(family == "gp"){
    return(-log(par[1]) - par[2] - 1)
  } else if(family == "gomp"){
    stopifnot("Install package \"gsl\" to use \"lprior_mdi_elife\" with the Gompertz model.\n Try `install.packages(\"gsl\")`" = requireNamespace("gsl", quietly = TRUE))
    e1 <- exp(1/par[2])*gsl::expint_E1(x = 1/par[2])
    #The limit exp(1/b)*E1(1/b) = 0 as b -> 0
    return(-log(par[1]) + ifelse(is.na(e1),0,e1))
  }
}

#' Box-Cox transformation function
#'
#' Given a vector of parameters, apply the Box-Cox transformation.
#'
#' @export
#'@keywords internal
boxcox_transfo <- function(par, lambda = rep(1, length(par))){
  stopifnot(length(par) == length(lambda),
            is.numeric(par),
            is.numeric(lambda))
  ifelse(lambda == 0, log(par), (par^lambda-1)/lambda)
}
# Posterior density for selected model using
# rust and the maximal data information prior

#' Log posterior distribution with MDI priors
#'
#' Log of the posterior distribution for excess lifetime
#' distribution with maximal data information priors.
#' @export
#' @inheritParams nll_elife
lpost_elife <- function(par,
                        time,
                        time2 = NULL,
                        event = NULL,
                        type = c("right","left","interval","interval2"),
                        ltrunc = NULL,
                        rtrunc = NULL,
                        family = c("exp","gp","gomp"),
                        thresh = 0,
                        weights = rep(1, length(time)),
                        status = NULL,
                        ...){
  type <- match.arg(type)
  family <- match.arg(family)
  obound <- switch(family,
                   exp = par[1] <= 0,
                   gp = par[1] <= 0 | par[2] < -1,
                   gomp = par[1] <=0 | par[2] <= 0)
  if(obound){
    return(-1e20)
  }
loglik  <- -nll_elife(par = par,
          time = time,
          time2 = time2,
          event = event,
          type = type,
          ltrunc = ltrunc,
          rtrunc = rtrunc,
          thresh = thresh,
          family = family,
          weights = weights,
          status = status)
  logprior <- lprior_mdi_elife(family = family, par = par)
  lpost <- loglik + logprior
  ifelse(is.finite(lpost), lpost, -1e20)
}

