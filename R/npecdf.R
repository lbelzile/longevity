#' Nonparametric estimation of the survival function
#'
#' The survival function is obtained through the EM algorithm
#' described in Turnbull (1976) for left-truncated right-censored
#' or else doubly truncated observations; censoring is assumed to be
#' non-informative. The survival function changes only
#' at the J distinct exceedances \eqn{y_i-u} and truncation points.
#'
#' The unknown parameters of the model are \eqn{p_j (j=1, \ldots, J)}
#' subject to the constraint that \eqn{\sum_{j=1}^J p_j=1}.
#'
#' @param dat vector of exceedances of length \code{n}, i.e., observations that are larger than the thresh
#' @param thresh double thresh
#' @param rcens logical vector of length \code{n}, where \code{rcens} is \code{TRUE} for right-censored observation and \code{FALSE} otherwise.
#' @param ltrunc vector of size \code{n} of lower truncation time
#' @param rtrunc vector of size \code{n} of upper truncation time
#' @param tol double, relative tolerance for convergence of the EM algorithm
#' @param vcov logical; should the observed information matrix be computed? Default to \code{FALSE}
#' @return an object of class \code{stepfun}
#' @export
#' @examples
#' set.seed(2021)
#' n <- 20L
#' # Create fake data
#' ltrunc <- pmax(0, runif(n, -0.5, 1))
#' rtrunc <- runif(n, 6, 10)
#' dat <- r_dtrunc_elife(n = n, scale = 1,
#'                shape = -0.1,
#'                lower = ltrunc,
#'                upper = rtrunc,
#'                family = "gp")
#' thresh <- 0
#' npi <- np.ecdf(dat = dat, thresh = thresh,
#' rtrunc = rtrunc, ltrunc = ltrunc)
np.ecdf <- function(dat,
                    thresh,
                    rcens = NULL,
                    ltrunc = NULL,
                    rtrunc = NULL,
                    tol = 1e-12,
                    vcov = FALSE) {
  stopifnot(
    "Data must be either right-censored or right-truncated." = is.null(rcens) |
      is.null(rtrunc),
    "Argument `thresh` is missing." = !missing(thresh),
    "Argument `thresh` should be positive." = min(thresh) >= 0,
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
      "`ltrunc` should be smaller or equal to `dat`" = isTRUE(all(ltrunc <= dat)),
      "`ltrunc` must be larger than the thresh" = isTRUE(min(ltrunc) >= thresh[1])
    )
  }
  if (!is.null(ltrunc) & !is.null(rtrunc)) {
    "`ltrunc` should be smaller than `rtrunc`" = isTRUE(all(ltrunc < rtrunc))
  }
  # Two scenarios:
  # (a) (left-truncated) and right-censored data
  # handled via KM estimator
  # TODO
  # (b) interval truncated
  # mdat <- max(dat)
  unex <- sort(unique(dat) - thresh[1])
  J <- length(unex)
  if (length(unex) < 2) {
    stop("Only one unexal: consider increasing the size of `delta`.")
  }
  # Find in which interval the observation lies
  omega_lb <- findInterval(x = dat - thresh[1],
                           vec = unex,
                           all.inside = FALSE)
  if (!is.null(rcens)) {
    # For right-censored data,
    # this should be 1 for all subsequent intervals
    omega_ub <- ifelse(rcens, J, omega_lb)
  } else{
    # otherwise, it is the unexal in which death occurs
    omega_ub <- omega_lb
  }
  if (!is.null(ltrunc)) {
    eta_lb <-
      findInterval(x = pmax(0, ltrunc - thresh[1]),
                   vec = unex,
                   all.inside = TRUE)
  } else{
    eta_lb <- rep(1L, n)
  }
  if (!is.null(rtrunc)) {
    eta_ub <-
      findInterval(x = rtrunc - thresh[1],
                   vec = unex,
                   all.inside = FALSE)
  } else{
    eta_ub <- rep(J, n)
  }
#
#   start <- rep(1 / J, J - 1)
#   copt <- optim(
#     par = logit(start),
#     fn = np_nll,
#     omega_lb = omega_lb,
#     omega_ub = omega_ub,
#     eta_lb = eta_lb,
#     eta_ub = eta_ub,
#     method = "BFGS",
#     control = list(maxit = 1e6, reltol = 1e-10)
#   )
#   sol <- c(expit(copt$par), 1 - sum(expit(copt$par)))
  ## EM algorithm
  # Initialize EM
  # p <- sol
  p <- rep(1/J, J)
  pnew <- rep(0, J)
  C_mat <- D_mat <- matrix(0, nrow = n, ncol = J)
  if (!is.null(rcens)) {
    which_rcens <- which(rcens)
  }
  n_iter <- 0L
  n_iter_max <- 1e4L
  for (i in 1:n) {
    # this doesn't need to be updated
    C_mat[i, omega_lb[i]] <- 1
  }
  while (n_iter < n_iter_max) {
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
    pnew[pnew < 1e-14 * n] <- 0
  if(max(abs(p - pnew)) < tol && n_iter_max != n_iter) {
    n_iter_max <- n_iter + 1L
  }
  p <- pnew
  }
  if(vcov){
  # Compute the hessian matrix of the log-likelihood
  uij <- C_mat + D_mat
  Omega_mat <- Eta_mat <- matrix(FALSE, nrow = n, ncol = J)
  for(i in 1:n){
    Omega_mat[i,omega_lb[i]:omega_ub[i]] <- TRUE
    Eta_mat[i,eta_lb[i]:eta_ub[i]] <- TRUE
  }
  sum_Omega_i_sq <- as.numeric(Omega_mat %*% p)^2
  sum_Eta_i_sq <- as.numeric(Eta_mat %*% p)^2
  infomat <- matrix(0, J-1, J-1)
  for(j in 1:(J-1)){
    for(k in j:(J-1)){
    infomat[j,k] <- -sum((Omega_mat[,k]-Omega_mat[,J])*(Omega_mat[,j]-Omega_mat[,J])/sum_Omega_i_sq) +
      sum((Eta_mat[,k]-Eta_mat[,J])*(Eta_mat[,j]-Eta_mat[,J])/sum_Eta_i_sq)
    infomat[k,j] <- infomat[j,k]
    }
  }
  covmat <- try(solve(-infomat), silent = TRUE)
  if(is.character(covmat)){
    covmat <- NULL
  }
  } else{
    covmat <- NULL
  }
  list(x = unex,
       par = p,
       vcov = covmat)
  # TODO
  # determine how to obtain confidence intervals
  #  Delta-method gives variance for the sum - can be used for normal intervals
  #  (a) Hall-Wellner (1980) type bands
  #  (b) Some methods based on Monte Carlo methods for Brownian motion (Kendall, Marin, Robert)
  #  (c) empirical likelihood?
}

#' Create a survival function using weights and intervals
#' @param npi output from np.surv, a list with components \code{interval}, \code{par} and \code{vcov}
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


#' Marginal log likelihood function of the nonparametric multinomial with censoring and truncation
#' @param p vector of \code{D-1} parameters
#' @param omega_lb index of interval in which death occurs
#' @param omega_ub index of interval in which death occurs (if death is observed), or else the largest interval.
#' @param eta_lb vector of largest index for the lower truncation
#' @param eta_ub vector of smallest index for the upper truncation
#' @param transform logical; are parameters on logit scale? Default to \code{TRUE}
#' @keywords internal
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

#' Inverse logistic link function
#'
#' @param x vector of real values
#' @return probability vector
#' @keywords internal
expit <- function(x) {
  1 / (1 + exp(-x))
}
#' Logistic link function
#'
#' @param x vector of probabilities
#' @return vector of real values
#' @keywords internal
logit <- function(x) {
  log(x) - log(1 - x)
}
