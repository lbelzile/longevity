#' Piece-wise generalized Pareto distribution
#'
#' Density, distribution, quantile functions and random number generation from the
#' mixture model of Northrop and Coleman (2014), which consists of \code{m}
#' different generalized Pareto distributions over non-overlapping intervals
#' with \code{m} shape parameters and one scale parameter; the other scale parametrs are
#' constrained so that the resulting distribution is continuous over the domain
#' and reduces to a generalized Pareto distribution if all of the shape parameters are equal.
#'
#' @name gppiece
#' @references Northrop & Coleman (2014). Improved threshold diagnostic plots for extreme value
#' analyses, \emph{Extremes}, \bold{17}(2), 289--303.
#' @param x,q vector of quantiles
#' @param scale positive value for the first scale parameter
#' @param shape vector of \code{m} shape parameters
#' @param thresh vector of \code{m} thresholds
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{\Pr[X \leq x]} otherwise, \eqn{\Pr[X > x]}.
#' @param log,log.p logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#'
NULL

#' @rdname gppiece
#' @export
#' @keywords internal
dgppiece <- function(x,
                     scale,
                     shape,
                     thresh,
                     log = FALSE
                     ){
  #parametrization is (sigma_1, xi_1, ..., xi_m)

  stopifnot("Thresholds should be a numeric vector" = is.numeric(thresh),
            "Threshold should be sorted in increasing order" = !is.unsorted(thresh),
            "Threshold must be positive" = isTRUE(all(thresh >= 0)),
            "There should be as many shape parameters as there are thresholds." = length(shape) == length(thresh),
            "There should be only a single scale parameter" = length(scale) == 1L
            )
  m <- length(shape)
  w <- as.numeric(diff(thresh))
  scale <- scale + c(0, cumsum(shape[-m]*w))
  stopifnot("At least one scale parameter is negative" = isTRUE(all(scale > 0)))
  logp <- c(0, cumsum(ifelse(shape[-m] == 0,
                             -w/scale[-m],
                             -1/shape[-m]*log(pmax(0,1+shape[-m]*w/scale[-m])))
                      )
            )
  Ind <- as.integer(cut(x, c(thresh, Inf), include.lowest = TRUE))
  ldens <-  logp[Ind] - log(scale[Ind]) -
      ifelse(shape[Ind] == 0,
             (x-thresh[Ind])/scale[Ind],
      (1+1/shape[Ind])*log(pmax(0,1+shape[Ind]*(x-thresh[Ind])/scale[Ind]))
      )

   if(log){
     return(ldens)
   } else{
     return(exp(ldens))
   }
}

#' @rdname gppiece
#' @export
#' @keywords internal
pgppiece <- function(q,
                     scale,
                     shape,
                     thresh,
                     lower.tail = TRUE,
                     log.p = FALSE
){
  stopifnot("Thresholds should be a numeric vector" = is.numeric(thresh),
            "Threshold should be sorted in increasing order" = !is.unsorted(thresh),
            "Threshold must be positive" = isTRUE(all(thresh >= 0)),
            "There should be as many shape parameters as there are thresholds." = length(shape) == length(thresh),
            "There should be only a single scale parameter" = length(scale) == 1L
  )
  m <- length(shape)
  w <- as.numeric(diff(thresh))
  scale <- scale + c(0, cumsum(shape[-m]*w))
  stopifnot("At least one scale parameter is negative" = isTRUE(all(scale > 0)))
  logp <- c(0, cumsum(ifelse(shape[-m] == 0,
                             -w/scale[-m],
                             -1/shape[-m]*log(pmax(0,1+shape[-m]*w/scale[-m])))
                      )
            )
  prob_interv <- -diff(c(exp(logp),0))
  Ind <- as.integer(cut(q, c(thresh, Inf), include.lowest = TRUE))
  cum_prob <- 1-exp(logp)
  F_upper <- c(vapply(1:(m-1), function(i){pgpd(q = w[i], scale = scale[i], shape = shape[i])}, numeric(1)), 1)
  ret_p <- cum_prob[Ind] +
    prob_interv[Ind]*(1 - ifelse(shape[Ind] == 0,
              exp((thresh[Ind] - q) / scale[Ind]),
              pmax(0, (1 + shape[Ind] * (q - thresh[Ind]) / scale[Ind]))^(-1/shape[Ind])
    )) / F_upper[Ind]
  if(shape[m] < 0){
    # check that any point beyond the support gets cumulative value of 1
    outbound <- which(q[Ind == m] > thresh[m] - scale[m]/shape[m])
    if(length(outbound) > 0){
      ret_p[outbound] <- 1
    }
  }
  outzero <- which(q <= min(thresh))
  if(length(outzero) > 0){
    ret_p[outzero] <- 0
  }
  if(!lower.tail){
    ret_p <- 1-ret_p
  }
  if(log.p){
    return(log(ret_p))
  } else{
    return(ret_p)
  }
}

#' @rdname gppiece
#' @export
#' @keywords internal
qgppiece <- function(p,
                     scale,
                     shape,
                     thresh,
                     lower.tail = TRUE,
                     log.p = FALSE
){
  stopifnot("Thresholds should be a numeric vector" = is.numeric(thresh),
            "Threshold should be sorted in increasing order" = !is.unsorted(thresh),
            "Threshold must be positive" = isTRUE(all(thresh >= 0)),
            "There should be as many shape parameters as there are thresholds." = length(shape) == length(thresh),
            "There should be only a single scale parameter" = length(scale) == 1L
  )
  m <- length(shape)
  w <- as.numeric(diff(thresh))
  scale <- scale + c(0, cumsum(shape[-m]*w))
  stopifnot("At least one scale parameter is negative" = isTRUE(all(scale > 0)))
  logp <- c(0, cumsum(ifelse(shape[-m] == 0,
                             -w/scale[-m],
                             -1/shape[-m]*log(pmax(0,1+shape[-m]*w/scale[-m])))
                )
            )
  cum_prob_interv <- c(1-exp(logp), 1)
  F_upper <- c(vapply(1:(m-1), function(i){pgpd(q = w[i], scale = scale[i], shape = shape[i])}, numeric(1)), 1)
  if(log.p){
    p <- exp(p)
  }
  stopifnot("Some probabilities are smaller than zero or larger than one." = isTRUE(all(c(p >= 0, p <= 1))))
  if(!lower.tail){
    p <- 1 - p
  }
  Ind <- as.integer(cut(p, cum_prob_interv, include.lowest = TRUE))
  # Because of the division and numerical overflow, 1 gets mapped to *nearly* one
  ps <-  ifelse(p!=1, (p - (1-exp(logp)[Ind]))/(-diff(c(exp(logp),0))[Ind])*F_upper[Ind], 1)
  thresh[Ind] +
            ifelse(shape[Ind] == 0,
                  - scale[Ind] * log(1-ps),
                  scale[Ind] * ((1-ps)^(-shape[Ind]) - 1)/shape[Ind]
  )

}

#' @rdname gppiece
#' @export
#' @keywords internal
rgppiece <- function(n,
                     scale,
                     shape,
                     thresh
){
# Compute the probability p_j, then p_j - p_{j+1} is the
# probability of falling in the interval (v_j, v_{j+1})
# 1) simulate a multinomial vector
# 2) sample a truncated generalized Pareto on (v_j, v_{j+1})
# 3) add back the lowest threshold
  stopifnot("Thresholds should be a numeric vector" = is.numeric(thresh),
            "Threshold should be sorted in increasing order" = !is.unsorted(thresh),
            "Threshold must be positive" = isTRUE(all(thresh >= 0)),
            "There should be as many shape parameters as there are thresholds." = length(shape) == length(thresh),
            "There should be only a single scale parameter" = length(scale) == 1L
  )
  m <- length(shape)
  w <- as.numeric(diff(thresh))
  scale <- scale + c(0, cumsum(shape[-m]*w))
  stopifnot("At least one scale parameter is negative" = isTRUE(all(scale > 0)))
  logp <- c(0, cumsum(ifelse(shape[-m] == 0,
                             -w/scale[-m],
                             -1/shape[-m]*log(pmax(0,1+shape[-m]*w/scale[-m])))
                      )
            )
  prob_interv <- -diff(c(exp(logp),0))
  #sample.int(m, prob = prob_interv, replace = TRUE, size = n)
  n_mixt <- rmultinom(n = 1, size = n, prob = prob_interv)
  samp <- vector(mode = "numeric", length = n)
  pos <- 1L
  for(j in 1:(m-1)){
    if(n_mixt[j] > 0){
      samp[pos:(pos+n_mixt[j]-1)] <- thresh[j] +
        qgpd(runif(n_mixt[j])*pgpd(q = w[j], scale = scale[j], shape = shape[j]),
             scale = scale[j], shape = shape[j])
    }
    pos <- pos + n_mixt[j]
  }
  if(n_mixt[m] > 0){
    samp[pos:n] <- thresh[m] + qgpd(runif(n_mixt[m]), scale = scale[m], shape = shape[m])
  }
  return(samp[sample.int(n, n)])
}
