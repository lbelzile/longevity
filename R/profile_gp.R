#' #' Profile likelihood for the endpoint of the generalized Pareto distribution
#' #'
#' #' @inheritParams hazard_elife
#' #' @param psi optional vector of endpoints at which to compute the profile
#' #' @return a vector of length three containing twice the negative log-likelihood value, the endpoint value and the maximum of the nuisance lambda (i.e., the shape parameter).
#' prof_gp_endpoint <- function(dat,
#'                              thresh = 0,
#'                              type = c("none","ltrt","ltrc"),
#'                              rcens,
#'                              ltrunc,
#'                              rtrunc,
#'                              weights = rep(1, length(dat)),
#'                              level = 0.95,
#'                              psi = NULL,
#'                              plot = FALSE){
#'   stopifnot("The level of the confidence interval must be in [0,1]" = level > 0 && level < 1,
#'             "Level should be a numeric of length one" = length(level) == 1L,
#'   )
#'   type <- match.arg(type)
#'   for(thresh )
#'   mle <- optim_elife(dat = dat,
#'                      thresh = thresh,
#'                      ltrunc = ltrunc,
#'                      rtrunc = rtrunc,
#'                      rcens = rcens,
#'                      type = type,
#'                      family = "gp",
#'                      weights = weights)
#'   if(mle$par['shape'] > 0){
#'     stop("The maximum likelihood estimator for the endpoint is infinite.")
#'   }
#'   opt <- optim(fn = function(xi){
#'     sigma <- -endpoint*xi
#'     gpd_cens(par = c(sigma, xi), dat = dat, rightcens = rightcens, slow = slow)},
#'     method = "Brent", par = -0.1,
#'     lower = -1,
#'     upper = -1e-8,
#'     control = list(reltol=1e-12))
#'   res <- -2*opt$value
#'   return(c(res, -opt$par*endpoint, opt$par))
#' }
#'
#' #' Profile likelihood for the endpoint of the generalized Pareto distribution
#' #'
#' #' Profile likelihood for the endpoint of the
#' #' generalized Pareto distribution
#' #' for left-truncated and right-censored observations.
#' #'
#' #' @inheritsParam gpd_dtrunc
#' #' @param endpoint value of the endpoint at which to compute the profile
#' #' @return a vector of length three containing twice the negative log-likelihood value, the endpoint value and the maximum of the nuisance lambda (i.e., the shape parameter).
#' prof_gpd_dtrunc_endpoint <- function(endpoint, dat, slow, supp){
#'   stopifnot(length(endpoint) == 1L)
#'   opt <- optim(fn = function(xi){
#'     sigma <- -endpoint*xi
#'     gpd_dtrunc(par = c(sigma, xi), dat = dat, supp = supp, slow = slow)},
#'     method = "Brent", par = -0.1,
#'     lower = -1,
#'     upper = -1e-8,
#'     control = list(reltol=1e-12))
#'   res <- -2*opt$value
#'   return(c(res, -opt$par*endpoint, opt$par))
#' }
#'
#' #' Profile likelihood for shape parameter of the generalized Pareto distribution
#' #'
#' #' #' This function implements the profile likelihood with left-truncated and right-censored observations
#' #' and must be applied repeatedly for each threshold
#' #' @inheritParams gpd_cens
#' #' @return point estimate and confidence limits for the shape parameter
#' prof_gpd_cens_xi_confint <- function(dat, rightcens, slow){
#'   xis <- seq(-0.5, 0.5, length = 200)
#'   mdat <- max(dat)
#'   dev <- sapply(xis, function(xi){
#'     opt <- optimize(f = function(lambda){gpd_cens(par = c(lambda, xi),
#'                                                   dat = dat,  rightcens = rightcens, slow = slow)},
#'                     interval = c(ifelse(xi < 0, mdat*abs(xi), 1e-8), 1e9), tol = 1e-10)
#'     c(-2*opt$objective, opt$minimum)
#'   })
#'   ms <- which.max(dev[1,])
#'   opt_mle <- optim(f = function(theta){gpd_cens(par = theta, dat = dat,  rightcens = rightcens, slow = slow)},
#'                    par = c(dev[2,ms],xis[ms]), method = "N", hessian = TRUE,
#'                    control = list(parscale = c(500, 0.01), reltol = 1e-10, maxit = 1000))
#'
#'   mle <- opt_mle$par
#'   maxll <- -2*opt_mle$value
#'   prof <- list(psi = xis, pll = dev[1,]-maxll, maxpll = 0, mle = mle,
#'                psi.max = mle[2], std.error = sqrt(solve(opt_mle$hessian)[2,2]))
#'   confint <- confint_int(prof, parm = "profile")
#'   confint
#'   # return(c(res, psi, opt$minimum))
#' }
#'
#'
#' #' Likelihood root function for shape parameter of the generalized Pareto distribution
#' #'
#' #' This function returns the likelihood root of the profile log-likelihood for the shape of the generalized Pareto distribution with left-truncated and right-censored data.
#' #' Specifically, \eqn{-r^2/2} is the profile likelihood and the two-sided p-value is\code{qchisq(p, 1)/2}.
#' #'
#' #' @param psi value of the shape at which to compute the p-value
#' #' @param thetahat maximum likelihood estimates of the scale and shape parameters
#' #' @inheritParams gpd_cens
#' #' @return a p-value
#' rfun_xi_gpd_cens <- function(psi, thetahat, dat, rightcens, slow){
#'   if(abs(psi) > 1e-6){
#'     llp <- prof_gpd_cens_xi(xi = psi, dat = dat, rightcens = rightcens, slow = slow)
#'   } else{
#'     rate <- exp_mle_lt_rc(dat = dat, rightcens = rightcens, slow = slow)
#'     llp <- c(-2*gpd_cens(par = c(rate, 0), dat = dat,  rightcens = rightcens, slow = slow), 0)
#'   }
#'   -2*gpd_cens(par = thetahat, dat = dat, rightcens = rightcens, slow = slow) - llp[1]
#' }
#'
#'
#'
#' #' Profile deviance of exponential sub-model
#' #'
#' #' This function computes the profile when the shape is zero and returns
#' #' the value of twice the negative log-likelihood at that value
#' #' @inheritParams gpd_cens
#' #' @return twice the negative profile log-likelihood
#' nllhxizero <- function(dat, rightcens, slow){
#'   opt <- optimize(f = function(lambda){gpd_cens(par = c(lambda, 0),
#'                                                 dat = dat,  rightcens = rightcens, slow = slow)},
#'                   interval = c(1e-8, 1e9), tol = 1e-10)
#'   c(-2*opt$objective)
#' }
