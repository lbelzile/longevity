#' Profile likelihood for the endpoint of the generalized Pareto distribution
#'
#' This function returns
#'
#' @export
#' @inheritParams hazard_elife
#' @param psi mandatory vector of endpoints at which to compute the profile
#' @param confint logical; if \code{TRUE}, return a \code{level} confidence interval instead of a list with the profile log-likelihood components
#' @param level numeric; the level for the confidence intervals
#' @return a list with the mle of the endpoint and the profile log-likelihood
prof_gp_endpt <- function(time,
                          time2 = NULL,
                          event = NULL,
                          thresh = 0,
                          type = c("right","left","interval","interval2"),
                          ltrunc = NULL,
                          rtrunc = NULL,
                          weights = rep(1, length(time)),
                          psi = NULL,
                          confint = FALSE,
                          level = 0.95){
  stopifnot("Endpoints must be positive" = all(psi > thresh))
  psi <- psi - thresh[1]
  type <- match.arg(type)
  # Compute the maximum log-likelihood
  mle <- fit_elife(time = time,
                   time2 = time2,
                   event = event,
                   thresh = thresh,
                   ltrunc = ltrunc,
                   rtrunc = rtrunc,
                   type = type,
                   family = "gp",
                   weights = weights)
  np <- length(psi)
  param <- matrix(nrow = np, ncol = 2L)
  param[,1] <- psi
  colnames(param) <- c("endpt","shape")
  ll <- vector(mode = "numeric", np)
  opt_fun <- function(xi, endpoint){
    sigma <- -endpoint*xi
    nll_elife(par = c(sigma, xi),
              time = time,
              time2 = time2,
              event = event,
              weights = weights,
              thresh = thresh,
              type = type,
              ltrunc = ltrunc,
              rtrunc = rtrunc,
              family = "gp")
  }
  for(i in seq_along(psi)){
    opt <- optim(fn = opt_fun,
      method = "Brent",
      par = -0.1,
      lower = -1,
      upper = -1e-8,
      control = list(reltol=1e-12),
      endpoint = psi[i])
    ll[i] <- -opt$value
    param[i,2] <- opt$par
  }
  mle_endpt <- as.numeric(thresh + ifelse(mle$par[2] >= 0,
                                          Inf,
                                          -mle$par[1]/mle$par[2]))
  prof <- structure(
          list(psi = psi + thresh,
               lambda = param[,2],
               pll = ll,
               maxpll = mle$loglik,
               mle = mle_endpt,
               psi.max = mle_endpt,
               nexc = mle$nexc),
          class = "elife_profile")
  if(is.finite(prof$mle) & confint){
    confint <- try(conf_interv(prof, level = level))
    return(confint)
  } else{
    return(prof)
  }
}
#
# plot.elife_profile <- function()
