#' Nonparametric estimation of the survival function
#'
#' The survival function is obtained through the EM algorithm
#' described in Turnbull (1978) for left-truncated right-censored
#' or else doubly truncated observations; censoring is assumed to be
#' non-informative. The exceedances \eqn{y_i-u}
#' are separated into discretized intervals of the form
#' \deqn{p_j = \Pr[(j-1)\delta \leq Y-u < j\delta], \quad j=1, \ldots, J}
#' where \eqn{J} is chosen such that \eqn{\max(Y_j) - u < J \delta}.
#'
#' The unknown parameters of the model are \eqn{p_j (j=1, \ldots, J)}
#' subject to the constraint that \eqn{\sum_{j=1}^J p_j=1}.
#'
#' @param dat vector of exceedances of length \code{n}, i.e., observations that are larger than the thresh
#' @param thresh double thresh
#' @param rcens logical vector of length \code{n}, where \code{rcens} is \code{TRUE} for right-censored observation and \code{FALSE} otherwise.
#' @param ltrunc vector of size \code{n} of lower truncation time
#' @param rtrunc vector of size \code{n} of upper truncation time
#' @param delta width of the intervals
#' @param tol double, relative tolerance for convergence of the EM algorithm
#' @return an object of class \code{stepfun}
#' @examples
#' set.seed(2021)
#' n <- 100L
#' dat <- evd::rgpd(n = n, scale = 15, shape = -0.1)
#' ltrunc <- rep(0, n)
#' rtrunc <- rpois(n = n, lambda = 50) + dat
#' rcens <- NULL
#' thresh <- 0
#' delta <- 2
#' npi <- np.interv.mle(dat = dat, thresh = thresh,
#' delta = delta, rtrunc = rtrunc, ltrunc = ltrunc)
np.interv.mle <- function(dat,
                               thresh,
                               delta,
                               rcens = NULL,
                               ltrunc = NULL,
                               rtrunc = NULL,
                               tol = 1e-12) {
  #TODO include interval censoring?
  stopifnot(
    "Data must be either right-censored or right-truncated." = is.null(rcens) |
      is.null(rtrunc),
    "Argument `delta` is missing." = !missing(delta),
    "Argument `delta` should be of length one." = length(delta) == 1L,
    "Argument `delta` should be positive." = delta > 0,
    "Argument `thresh` is missing." = !missing(thresh),
    "Argument `thresh` should be positive." = thresh >= 0,
    "Argument `thresh` should be of length one." = length(thresh) == 1L,
    "Argument `dat` missing." = !missing(dat),
    "Argument `dat` should be a vector" = is.vector(dat),
    "Argument `tol` should be of length one." = length(tol) == 1L,
    "Argument `tol` should be positive." = isTRUE(tol > 0),
    "Argument `tol` should be smaller than 0.01." = isTRUE(tol < 0.01)
  )

  n <- length(dat)
  if (!is.null(rcens)) {
    stopifnot(
      "`rcens` should be a vector" = is.vector(rcens),
      "`rcens` should be the same length as `dat`" =  length(rcens) == n,
      "`rcens` should be of type `logical`" = is.logical(rcens)
    )
    rcens <- as.logical(rcens)
  }
  if (!is.null(rtrunc)) {
    stopifnot(
      "`rtrunc`` should be a vector" = is.vector(rtrunc),
      "`rtrunc` should be the same length as `dat`" =  length(rtrunc) == n,
      "`rtrunc` should be smaller or equal to `dat`" = isTRUE(all(rtrunc > dat))
    )
  }
  if (!is.null(ltrunc)) {
    stopifnot(
      "`ltrunc` should be a vector" = is.vector(ltrunc),
      "`ltrunc` should be the same length as `dat`" =  length(ltrunc) == n,
      "`ltrunc` should be greater or equal to `dat`" = isTRUE(all(ltrunc <= dat)),
      "`ltrunc` must be larger than the thresh" = isTRUE(min(ltrunc) >= thresh)
    )
  }
  if (!is.null(ltrunc) & !is.null(rtrunc)) {
    "`ltrunc` should be smaller than `rtrunc`" = isTRUE(all(ltrunc < rtrunc))
  }
  # Two scenarios:
  # (a) left-truncated and right-censored data
  # (b) left-truncated and right-truncated data
  J <- ceiling((max(dat) - thresh) / delta)
  interv <- thresh + seq(0, J * delta, by = delta)
  if (length(interv) < 2) {
    stop("Only one interval: consider increasing the size of `delta`.")
  }
  # Find in which interval the observation lies
  omega_lb <- findInterval(x = dat,
                           vec = interv,
                           all.inside = TRUE)
  if (!is.null(rcens)) {
    # For right-censored data,
    # this should be 1 for all subsequent intervals
    omega_ub <- ifelse(rcens, J, omega_lb)
  } else{
    # otherwise, it is the interval in which death occurs
    omega_ub <- omega_lb
  }
  if (!is.null(ltrunc)) {
    eta_lb <-
      findInterval(x = ltrunc - thresh,
                   vec = interv,
                   all.inside = TRUE)
  } else{
    eta_lb <- rep(1L, n)
  }
  if (!is.null(rtrunc)) {
    eta_ub <-
      findInterval(x = rtrunc - thresh,
                   vec = interv,
                   all.inside = TRUE)
  } else{
    eta_ub <- rep(J, n)
  }

  start <- rep(1 / J, J - 1)
  mb <- min(c(omega_lb, eta_lb))
  # Local optimum: if we have truncated data, the
  # EM algorithm won't update the components below
  # the lowest threshold
  if(!missing(ltrunc)){
    start[seq(length.out = mb - 1)] <- 1e-25
  }
  copt <- optim(
    par = logit(start),
    fn = np_nll,
    omega_lb = omega_lb,
    omega_ub = omega_ub,
    eta_lb = eta_lb,
    eta_ub = eta_ub,
    method = "BFGS",
    control = list(maxit = 1e6, reltol = 1e-10)
  )
  sol <- c(expit(copt$par), 1 - sum(expit(copt$par)))
  ## EM algorithm
  # Initialize EM
  p <- sol
  convergence <- FALSE


  C_mat <- D_mat <- matrix(0, nrow = n, ncol = J)
  if (!is.null(rcens)) {
    which_rcens <- which(rcens)
  }
  n_iter <- 0L
  for (i in 1:n) {
    # this doesn't need to be updated
    C_mat[i, omega_lb[i]] <- 1
  }
  while (!convergence) {
    n_iter <- n_iter + 1L
    # E-step
    # Update C_mat (only relevant for right-censored data)
    if (!is.null(rcens)) {
      for (i in 1:n) {
        if (rcens[i]) {
          ij <- omega_lb[i]:omega_ub[i]
          C_mat[i, ij] <- p[ij] / sum(p[ij])
        }
      }
    }
    # Update D_mat
    for (i in 1:n) {
      ij <- eta_lb[i]:eta_ub[i]
      D_mat[i, -ij] <- p[-ij] / sum(p[ij])
    }

    # M-step: maximize the log-likelihood
    pnew <- colSums(C_mat + D_mat)
    pnew <- pnew / sum(pnew)
    pnew[pnew < 1e-12 * n] <- 0
    if (max(p - pnew) < tol) {
      convergence <- TRUE
    }
    p <- pnew
  }
  # Compute the hessian matrix of the log-likelihood
  # but without the last parameter
  # and excluding zero weights
  #   np_nll_hessian(p = p[-length(p)],
  #                omega_lb = omega_lb,
  #                omega_ub = omega_ub,
  #                eta_lb = eta_lb,
  #                eta_ub = eta_ub)
  # np_nll_hessian_finite_diff(pf = p,
  #                            omega_lb = omega_lb,
  #                            omega_ub = omega_ub,
  #                            eta_lb = eta_lb,
  #                            eta_ub = eta_ub)
   hess <- try(np_nll_hessian(p = p[-length(p)],
                               omega_lb = omega_lb,
                               omega_ub = omega_ub,
                               eta_lb = eta_lb,
                               eta_ub = eta_ub))
  covmat <- matrix(0, length(p)-1, length(p)-1)
  wpos <- p[-length(p)] > 1e-8
  obs_info <- try(solve(hess[wpos, wpos]))
  if(!is.character(hess)){
    covmat[wpos, wpos] <- obs_info
  } else{
    covmat <- NULL
  }
  list(interval = interv,
       par = p,
       vcov = covmat)

  #TODO
  # determine how to obtain confidence intervals
  #  Delta-method gives variance for the sum - can be used for normal intervals
  #  (a) Hall-Wellner (1980) type bands
  #  (b) Some methods based on Monte Carlo methods for Brownian motion (Kendall, Marin, Robert)
  #  (c) empirical likelihood?
}

#' Create a survival function using weights and intervals
#' @param npi output from np.interv.mle, a list with components \code{interval}, \code{par} and \code{vcov}
#' @param conf.type transformation for the pointwise confidence intervals;
#' @param cumhaz plot the cumulative hazard rather than the probability in state or survival
#' @param ... additional arguments passed to plot
plot.np.surv <- function(npi,
                         conf.type = c("log", "log-log", "plain", "logit", "arcsin"),
                         cumhaz = FALSE,
                         ...){
  if(isTRUE(any(is.null(interval), is.null(par)))){
    stop("Intervals or vector of probability missing in `npi`.")
  }
  S_fn <- stepfun(x = npi$interval, y = 1 - c(0, cumsum(npi$par), 1), right = FALSE)
  np <- length(npi$par) - 1
  if(!is.null(npi$vcov)){
    var_S <- sapply(1:np, function(i){sum(npi$vcov[1:i, 1:i])})
  }

}


# Compute the hessian matrix using finite difference
np_nll_hessian_finite_diff <- function(pf,
                                       omega_lb,
                                       omega_ub,
                                       eta_lb,
                                       eta_ub) {
  # Unconstrained negative log likelihood function
  # This function takes as argument the full vector p
  # and rescales it so that the vector of probability
  # continues summing to one
  np_nll_unconst <- function(pf){
    stopifnot(isTRUE(all(pf >= 0)))
    pf <- pf/sum(pf)
    llp <- 0
    for (i in 1:length(omega_lb)) {
      llp <- llp + log(sum(pf[omega_lb[i]:omega_ub[i]])) - log(sum(pf[eta_lb[i]:eta_ub[i]]))
    }
    return(-llp)
  }
  covmat <- matrix(0, length(pf), length(pf))
  grad <- rep(0, length(pf))
  evec <- rep(0, length(pf))
  nllm <- np_nll_unconst(pf)
  for(i in 1:length(pf)){
    evec[i] <- 1e-7
    grad[i] <- np_nll_unconst(pf + evec)
    evec[i] <- 0
  }
  covmat[1,1] <- (np_nll_unconst(pf + 2*evec) - 2*grad[1] + nllm)/1e-14
  for(i in 2:length(pf)){
    evec[i] <- 1e-7
    covmat[i,i] <- (np_nll_unconst(pf + 2*evec) - 2*grad[i] + nllm)/1e-14
    for(j in 1:(i-1)){
      evec[j] <- 1e-7
      covmat[i,j] <-  covmat[j,i] <- (np_nll_unconst(pf + evec)-grad[j] - grad[i] + nllm)/1e-14
      evec[j] <- 0
    }
    evec[i] <- 0
  }
  return(covmat)
}
#' Hessian of (marginal), aka incomplete log likelihood function
#' @param p vector of D-1 parameters
np_nll_hessian <- function(p,
                           omega_lb,
                           omega_ub,
                           eta_lb,
                           eta_ub){
  p <- c(p, 1-sum(p))
  np <- length(p)
  n <- length(omega_lb)
  hessian <- matrix(0, np, np)
  for(i in 1:n){
    oseq <- omega_lb[i]:omega_ub[i]
    hessian[oseq, oseq] <- hessian[oseq, oseq] - 1/sum(p[oseq])^2
    eseq <- eta_lb[i]:eta_ub[i]
    hessian[eseq, eseq] <- hessian[eseq, eseq] + 1/sum(p[eseq])^2
  }
  return(-hessian[-np,-np])
}

#' Marginal log likelihood function of the nonparametric multinomial with censoring and truncation
#' @param p vector of \code{D-1} parameters
#' @param omega_lb index of interval in which death occurs
#' @param omega_ub index of interval in which death occurs (if death is observed), or else the largest interval.
#' @param eta_lb vector of largest index for the lower truncation
#' @param eta_ub vector of smallest index for the upper truncation
#' @param transform logical; are parameters on logit scale? Default to \code{TRUE}
np_nll <- function(p,
                   omega_lb,
                   omega_ub,
                   eta_lb,
                   eta_ub,
                   transform = TRUE) {
  if (transform) {
    p <- expit(p)
  }
  pf <- c(p, 1 - sum(p))
  if (isTRUE(any(c(pf > 1, pf < 0)))) {
    return(1e8)
  }
  llp <- 0
  for (i in 1:length(omega_lb)) {
    llp <- llp + log(sum(pf[omega_lb[i]:omega_ub[i]])) - log(sum(pf[eta_lb[i]:eta_ub[i]]))
  }
  return(-llp)
}

expit <- function(x) {
  1 / (1 + exp(-x))
}
logit <- function(x) {
  log(x) - log(1 - x)
}
