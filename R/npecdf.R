#' Nonparametric estimation of the survival function
#'
#' The survival function is obtained through the EM algorithm
#' described in Turnbull (1976) for left-truncated right-censored
#' or else doubly truncated observations; censoring is assumed to be
#' non-informative. The survival function changes only
#' at the \code{J} distinct exceedances \eqn{y_i-u} and truncation points.
#'
#' @note This function is a vanilla R implementation of the EM algorithm;
#' the function \link{npsurv} uses a backbone Cpp implementation
#' and will be faster for most settings (but it doesn't account for thresholds).
#'The function can return  the variance covariance matrix of the \code{J-1}
#' estimated probabilities, which is computed by inverting the negative
#' hessian of the incomplete log likelihood.
#'
#' The method is currently limited to 2500 unique failure time, since
#' the implementation has a heavy memory footprint of the order O(\eqn{J^2})
#'
#' The unknown parameters of the model are \eqn{p_j (j=1, \ldots, J)}
#' subject to the constraint that \eqn{\sum_{j=1}^J p_j=1}.
#' @inheritParams npsurv
#' @param thresh double thresh
#' @param tol double, relative tolerance for convergence of the EM algorithm
#' @param vcov logical; should the observed information matrix be computed? Default to \code{FALSE}
#' @return a list with elements
#' \itemize{
#' \item{\code{cdf}: }{right-continuous \code{stepfun} object defined by probabilities}
#' \item{\code{time}: }{unique failure exceedance times defining intervals}
#' \item{\code{prob}: }{\code{J} vector of probability of failure in interval}
#' \item{\code{vcov}: }{either \code{NULL} or an \code{J-1} by \code{J-1} covariance matrix for the first components}
#' \item{\code{niter}: }{number of iterations before EM algorithm convergence from equiprobable}
#' }
#' @useDynLib longevity, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @export
#' @examples
#' set.seed(2021)
#' n <- 20L
#' # Create fake data
#' ltrunc <- pmax(0, runif(n, -0.5, 1))
#' rtrunc <- runif(n, 6, 10)
#' dat <- samp_elife(n = n, scale = 1,
#'                shape = -0.1,
#'                lower = ltrunc,
#'                upper = rtrunc,
#'                family = "gp",
#'                type2 = "ltrt")
#' npi <- np_elife(time = dat,
#'                 rtrunc = rtrunc,
#'                 ltrunc = ltrunc,
#'                 vcov = TRUE)
np_elife <- function(time,
                     time2 = NULL,
                     event = NULL,
                     type = c("right","left","interval","interval2"),
                     thresh = 0,
                     ltrunc = NULL,
                     rtrunc = NULL,
                     tol = 1e-12,
                     vcov = FALSE) {
  stopifnot(
    "Argument `thresh` should be positive." = min(thresh) >= 0,
    "Argument `thresh` should be a vector of length 1." = length(thresh) == 1L,
    "Argument `time` missing." = !missing(time),
    "Argument `time` should be a vector" = is.vector(time),
    "Argument `tol` should be of length one." = length(tol) == 1L,
    "Argument `tol` should be positive." = isTRUE(tol > 0),
    "Argument `tol` should be smaller." = isTRUE(tol < 1e-4)
  )
  survout <- .check_surv(time = time,
                         time2 = time2,
                         event = event,
                         type = type)
  time <- survout$time
  time2 <- survout$time2
  status <- survout$status
  if(!is.null(ltrunc) && length(ltrunc) == 1L){
    ltrunc <- rep(ltrunc, length(time))
  }
  if(!is.null(rtrunc) && length(rtrunc) == 1L){
    rtrunc <- rep(rtrunc, length(time))
  }
  if(thresh[1] > 0){
    # Keep only exceedances, shift observations
    # We discard left truncated observations and interval censored
    # if we are unsure whether there is an exceedance
    ind <- ifelse(status == 2, FALSE, time > thresh[1])
    weights <- weights[ind] # in this order
    time <- time[ind] - thresh[1]
    time2 <- time2[ind] - thresh[1]
    status <- status[ind]
    if(!is.null(ltrunc)){
      ltrunc <- pmax(0, ltrunc[ind] - thresh[1])
    }
    if(!is.null(rtrunc)){
      rtrunc <- rtrunc[ind] - thresh[1]
    }
  }
  if(!is.null(ltrunc)){
    if(isTRUE(any(ltrunc > time, ltrunc > time2, na.rm = TRUE))){
      stop("Left-truncation must be lower than observation times.")
    }
  }
  if(!is.null(rtrunc)){
    if(isTRUE(any(rtrunc < time, rtrunc < time2, na.rm = TRUE))){
      stop("Right-truncation must be lower than observation times.")
    }
  }
  n <- length(time)
  cens <- !isTRUE(all(status == 1L, na.rm = TRUE))
  if (!is.null(rtrunc)) {
    stopifnot(
      "`rtrunc`` should be a vector" = is.vector(rtrunc),
      "`rtrunc` should be the same length as `time`" =  length(rtrunc) == n
    )
  }
  if (!is.null(ltrunc)) {
      stopifnot(
      "`ltrunc` should be a vector" = is.vector(ltrunc),
      "`ltrunc` should be the same length as `dat`" =  length(ltrunc) == n
    )
  }
  if (!is.null(ltrunc) & !is.null(rtrunc)) {
    stopifnot("`ltrunc` should be smaller than `rtrunc`" = isTRUE(all(ltrunc < rtrunc)))
  }
  if(is.null(ltrunc) && is.null(rtrunc)){
    trunc <- FALSE
  } else{
    trunc <- TRUE
  }
  # Two scenarios:
  if(cens){
    unex <- sort(unique(time[status == 1L]))
  } else{
    unex <- sort(unique(time))
  }
  J <- length(unex)
  if (length(unex) < 2) {
    stop("Only one interval: consider increasing the size of `delta`.")
  }
  if(length(J) > 2500L){
    stop("The R implementation is currently limited to 2500 unique failure times for memory reasons.")
  }
  # Find in which interval the observation lies
  cens_lb <- pmin(J-1, findInterval(
    x = time + (status %in% c(0L,3L))*1e-10,
    vec = unex,
    left.open = TRUE,
    all.inside = FALSE)) + 1L
  # get the first observation
  cens_lb[is.na(cens_lb)] <- 1L
  # values before first observation
  if(!cens){
    cens_ub <- cens_lb
  } else{
    cens_ub <- pmax(0, findInterval(x = time2 - (status %in% c(2,3))*1e-10,
                                    vec = unex,
                                    left.open = TRUE,
                                    all.inside = FALSE)) + 1L
    cens_ub[is.na(cens_ub)] <- J
  }

  if(is.null(ltrunc) & is.null(rtrunc)){
    trunc <- FALSE
    trunc_lb <- rep(1L, n)
    trunc_ub <- rep(J, n)
  } else{
    trunc <- TRUE
    if(!is.null(ltrunc)){
      trunc_lb <- findInterval(x = ltrunc,
                               vec = unex,
                               left.open = TRUE,
                               all.inside = FALSE) + 1L
    } else{
      trunc_lb <- rep(1L, n)
    }
    if(!is.null(rtrunc)){
      trunc_ub <- pmin(J, findInterval(x = rtrunc,
                               vec = unex,
                               left.open = TRUE,
                               all.inside = FALSE) + 1L)
    } else{
      trunc_ub <- rep(J, n)
    }
  }
  p <- rep(1/J, J)
  pnew <- rep(0, J)
  C_mat <- D_mat <- matrix(0, nrow = n, ncol = J)
  n_iter <- 0L
  n_iter_max <- 1e4L
  for (i in 1:n) {
    # this doesn't need to be updated
    C_mat[i, cens_lb[i]] <- 1
  }
  while (n_iter < n_iter_max) {
    n_iter <- n_iter + 1L
    # E-step
    # Update C_mat (only relevant for right-censored data)
    if (!cens){
      for (i in 1:n) {
          ij <- cens_lb[i]:cens_ub[i]
          C_mat[i, ij] <- p[ij] / sum(p[ij])
      }
    }
    # Update D_mat
    for (i in 1:n) {
      ij <- trunc_lb[i]:trunc_ub[i]
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
    Omega_mat[i,cens_lb[i]:cens_ub[i]] <- TRUE
    Eta_mat[i,trunc_lb[i]:trunc_ub[i]] <- TRUE
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
  survfun <- .wecdf(x = unex, w = p, type = "surv")
  std.err <- NULL
  if(vcov && !is.null(covmat)){
    if(is.matrix(covmat) & nrow(covmat) >= 1){
    std.err <- vapply(1:nrow(covmat), function(j){
      sum(covmat[j:nrow(covmat),j:nrow(covmat)])
      }, numeric(1))
    }
  }
  invisible(
    list(
       survfun = survfun,
       surv = 1-cumsum(p)[-length(p)],
       stdsurv = std.err,
       xval = unex,
       prob = p,
       vcov = covmat,
       niter = n_iter))
}

#' Marginal log likelihood function of the nonparametric multinomial with censoring and truncation
#' @param p vector of \code{D-1} parameters
#' @param cens_lb index of interval in which death occurs
#' @param cens_ub index of interval in which death occurs (if death is observed), or else the largest interval.
#' @param trunc_lb vector of largest index for the lower truncation
#' @param trunc_ub vector of smallest index for the upper truncation
#' @param transform logical; are parameters on logit scale? Default to \code{TRUE}
#' @keywords internal
np_nll <- function(p,
                   cens_lb,
                   cens_ub,
                   trunc_lb,
                   trunc_ub,
                   transform = TRUE) {
  if (transform) {
    p <- expit(p)
  }
  pf <- c(p, 1 - sum(p))
  if (isTRUE(any(c(pf > 1, pf < 0)))) {
    return(1e8)
  }
  llp <- 0
  for (i in seq_along(cens_lb)) {
    llp <- llp + log(sum(pf[cens_lb[i]:cens_ub[i]])) - log(sum(pf[trunc_lb[i]:trunc_ub[i]]))
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

#' Weighted empirical distribution function
#'
#' @param x vector of length \code{n} of values
#' @param w vector of weights of length \code{n}
#' @param type string, one of distribution function (\code{dist}) or survival function (\code{surv})
#' @author Adapted from spatstat (c) Adrian Baddeley and Rolf Turner
#' @keywords internal
.wecdf <- function(x, w, type = c("dist","surv")){
  # Adapted from spatstat (c) Adrian Baddeley and Rolf Turner
  type <- match.arg(type)
  stopifnot(length(x) == length(w))  #also returns error if x is multidimensional
  nbg <- is.na(x)
  x <- x[!nbg]
  w <- w[!nbg]
  n <- length(x)
  w <- w / sum(w)
  if (n < 1)
    stop("'x' must have 1 or more non-missing values")
  ox <- order(x)
  x <- x[ox]
  w <- w[ox]
  vals <- sort(unique(x))
  xmatch <- factor(match(x, vals), levels = seq_along(vals))
  wmatch <- tapply(w, xmatch, sum)
  wmatch[is.na(wmatch)] <- 0
  if(type == "dist"){
    rval <- approxfun(vals, cumsum(wmatch),
                      method = "constant",
                      yleft = 0,
                      yright = 1,
                      f = 0,
                      ties = "ordered")
  } else{
    rval <- approxfun(vals, c(1,1-cumsum(wmatch)[-length(wmatch)]),
                      method = "constant",
                      yleft = 1,
                      yright = 0,
                      f = 1,
                      ties = "ordered")
  }
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  invisible(rval)
}

#' @export
#' @keywords internal
.check_surv <- function(time,
                        time2 = NULL,
                        event = NULL,
                        type = c("right","left","interval","interval2")){
  stopifnot("User must provide a time argument" = !missing(time))
   time <- as.numeric(time)
   n <- length(time)
   type <- match.arg(type)
   if(!is.null(time2) && type %in% c("interval","interval2")){
     time2 <- as.numeric(time2)
     stopifnot("Both `time` and `time2` should be numeric vectors of the same length." = length(time2) == length(time))
   } else if(!is.null(time2)){
     stop("`time2` should only be provided for `type=interval`.")
   }
   if(type == "right"){
     #warning("Event should be provided, else all events are assumed to be observed.")
     if(is.null(event)){
       time2 <- time
       status <- rep(1L, length(time))
     } else{
       time2 <- ifelse(as.logical(event), time, NA)
       status <- ifelse(as.logical(event), 1L, 0L)
     }
   } else if(type == "left"){
     time2 <- time
     time <- ifelse(as.logical(event), time, NA)
     status <- ifelse(as.logical(event), 1L, 2L)
   } else if(type == "interval"){
     stopifnot("Event vector is missing" = !is.null(event))
     event <- as.integer(event)
     stopifnot("Event should be a vector of integers between 0 and 3" = isTRUE(all(event %in% 0:3)))
     time2 <- ifelse(event == 0L, NA, ifelse(event %in% c(1L,2L), time, time2))
     time <- ifelse(event == 2L, NA, time)
     status <- event
   } else if (type == "interval2") {
     stopifnot("`time2` must be a numeric vector." = is.numeric(time2),
               "`time` and `time2` must be two vectors of the same length" = length(time2) == length(time))
     time <- ifelse(is.finite(time), time, NA)
     time2 <- ifelse(is.finite(time2), time2, NA)
     unknown <- (is.na(time) & is.na(time2))
     if(any(unknown)){
       stop("Invalid data; time is not resolved")
     }
   }
  status <- ifelse(is.na(time), 2L, ifelse(is.na(time2), 0L, ifelse(time == time2, 1L, 3L)))
  return(list(time = time, time2 = time2, status = status))
}

#' Nonparametric maximum likelihood estimation for arbitrary truncation
#'
#' The syntax is reminiscent of the \link[survival]{Surv} function, with
#' additional vectors for left-truncation and right-truncation.
#'
#' @note Contrary to the Kaplan-Meier estimator, the mass is placed in the interval
#' [\code{max(time), Inf}) so the resulting distribution function is not deficient.
#'
#' @param time excess time of the event of follow-up time, depending on the value of event
#' @param event status indicator, normally 0=alive, 1=dead. Other choices are \code{TRUE}/\code{FALSE} (\code{TRUE} for death).
#' For interval censored data, the status indicator is 0=right censored, 1=event at time, 2=left censored, 3=interval censored.
#' Although unusual, the event indicator can be omitted, in which case all subjects are assumed to have an event.
#' @param time2 ending excess time of the interval for interval censored data only.
#' @param type character string specifying the type of censoring. Possible values are "\code{right}", "\code{left}", "\code{interval}", "\code{interval2}".
#' @param ltrunc lower truncation limit, default to \code{NULL}
#' @param rtrunc upper truncation limit, default to \code{NULL}
#' @seealso \code{\link[survival]{Surv}}
#' @return a list with components
#' \itemize{
#' \item{\code{xval}: }{unique ordered values of observed failure times}
#' \item{\code{prob}: }{estimated probability of failure}
#' \item{\code{convergence}: }{logical; \code{TRUE} if the EM algorithm iterated until convergence}
#' \item{\code{niter}: }{logical; number of iterations for the EM algorithm}
#' \item{\code{cdf}: }{weighted empirical distribution function}
#' }
#' @export
npsurv <- function(time,
                   time2 = NULL,
                   event = NULL,
                   type = c("right","left","interval","interval2"),
                   ltrunc = NULL,
                   rtrunc = NULL){
  stopifnot("`time` must be a numeric vector." = is.numeric(time))
  survout <- .check_surv(time = time,
              time2 = time2,
              event = event,
              type = type)
  time <- survout$time
  time2 <- survout$time2
  status <- survout$status
  n <- length(time)
    # unique survival time
  if(sum(status == 1L) < 2){
    warning("There are not enough uncensored observations to estimate the distribution.")
  }
  unex <- sort(unique(time[status == 1L]))
  J <- length(unex)
  if(isTRUE(any(status != 1L))){
    cens <- TRUE
  } else{
    cens <- FALSE
  }
  # Find in which interval the observation lies
  cens_lb <- pmin(J-1, findInterval(
                 x = time + (status %in% c(0L,3L))*1e-10,
                 vec = unex,
                 left.open = TRUE,
                 all.inside = FALSE))
  # get the first observation
  cens_lb[is.na(cens_lb)] <- 0L
  # values before first observation
  if(!cens){
    cens_ub <- cens_lb
  } else{
    cens_ub <- pmax(0, findInterval(x = time2 - (status %in% c(2,3))*1e-10,
                          vec = unex,
                          left.open = TRUE,
                          all.inside = FALSE))
    cens_ub[is.na(cens_ub)] <- J - 1L
  }
  if(is.null(ltrunc) & is.null(rtrunc)){
    trunc <- FALSE
    trunc_lb <- 0L
    trunc_ub <- J-1L
  } else{
    trunc <- TRUE
  if(!is.null(ltrunc)){
    trunc_lb <- rep(x = findInterval(x = ltrunc,
                             vec = unex,
                             left.open = TRUE,
                             all.inside = FALSE), length.out = n)
  } else{
    trunc_lb <- rep(0L, n)
  }
    if(!is.null(rtrunc)){
      trunc_ub <- rep(pmin(J-1, findInterval(x = rtrunc,
                               vec = unex,
                               left.open = TRUE,
                               all.inside = FALSE)), length.out = n)
    } else{
      trunc_ub <- rep(J-1, n)
    }
  }
  survfit_res <- .turnbull_em(x = unex,
               censLow = as.integer(cens_lb),
               censUpp = as.integer(cens_ub),
               truncLow = as.integer(trunc_lb),
               truncUpp = as.integer(trunc_ub),
               cens = cens,
               trunc = trunc,
               tol = 1e-12
  )
  if(!survfit_res$conv){
    warning(paste("EM algorithm did not converge after", survfit_res$niter, "iterations."))
  }
  structure(
    list(cdf = .wecdf(x = unex, w = survfit_res$p),
       xval = unex,
       prob = as.numeric(survfit_res$p),
       niter = survfit_res$neval,
       convergence = survfit_res$conv),
    class = "npcdf"
  )
}
