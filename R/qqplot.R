# For the RSOS paper, I did the following:
# Fit a Kaplan-Meier (left-truncated right-censored data)
# and transform these plotting positions to exponential
# For the uncertainty, bootstrap scheme:
# - Sample new observations from the model
# (doubly truncated data, people outside the window
# are right-censored)
# - Refit the parametric model and get MLE
# - Refit the product-limit estimator, extract timings
# - Use the output of the latter and compute a weighted empirical distribution function (ECDF)
# - Evaluate the ECDF on a regular grid of days
# - Transform these days to the exponential (with bootstrap param)
# Compute simultaneous confidence intervals using the envelope method

#' Goodness-of-fit plots for parametric models
#'
#' Because of censoring and truncation, the plotting
#' positions must be adjusted accordingly.
#' For right-censored data, the methodology is described
#' in Waller & Turnbull (1992). Only non-censored observations are
#' displayed, which can create distortion.
#'
#' For truncated data, we first estimate the distribution function
#' nonparametrically, \eqn{F_n}. The uniform plotting positions of the data
#' \deqn{v_i = [F_n(y_i) - F_n(a_i)]/[F_n(b_i) - F_n(a_i)]}.
#' For probability-probability plots, the empirical quantiles are transformed
#' using the same transformation, with \eqn{F_n} replaced by the postulated or estimated
#' distribution function \eqn{F_0}.
#' For quantile-quantile plots, the plotting positions \eqn{v_i} are mapped back
#' to the data scale viz. \deqn{F_0^{-1}\{F_0(a_i) + v_i[F_0(b_i) - F_0(a_i)]\}}
#' When data are truncated, the plotting positions need not be in the same order as the data
#'
#' @param object a parametric model of class \code{elife_par}
#' @param plot.type string, one of \code{base} for base R or \code{ggplot}
#' @param which.plot vector of string indicating the plots, among \code{pp} for probability-probability plot, \code{qq} for quantile-quantile plot, \code{erp} for empirically rescaled plot (only useful for censored data), \code{e} for exponential quantile-quantile plot
#' @param detrended logical; if \code{TRUE}, Tukey's mean difference plot. The pair \eqn{(x,y)} is mapped to \code{((x+y)/2,y-x)} are detrended (y-x, x)
plot.elife_par <- function(object,
                           plot.type = c("base","ggplot"),
                           which.plot = c("pp","qq"),
                           detrended = FALSE){
  if(detrended){
    warning("Option currently unimplemented.")
    detrended <- FALSE
  }
  stopifnot("Object should be of class `elife_par`" = inherits(object, what = "elife_par"))
     if(is.null(object$dat)){
          stop("`object` created using a call to `fit_elife` should include the data (`export=TRUE`).")
     }
  which.plot <- match.arg(which.plot, choices = c("pp","qq","exp"), several.ok = TRUE)
  # if("erp" %in% which.plot && is.null(rcens)){
  #   stop("Empirically rescaled plot is only useful for censored data.")
  # }
     thresh <- object$thresh
     exc <- object$dat > thresh[1]
     dat <- object$dat[exc] - thresh[1]
     if(length(dat) > 2000){
       warning("Nonparametric estimation of the function very expensive for the given sample size.")
     }
     if(!is.null(object$ltrunc)){
       ltrunc <- pmax(0, object$ltrunc[exc] - thresh[1])
     } else{
       ltrunc <- NULL
     }
     if(!is.null(object$rtrunc)){
       rtrunc <- object$rtrunc[exc] - thresh[1]
     } else{
       rtrunc <- NULL
     }
     if(!is.null(object$rcens)){
       rcens <- object$rcens[exc]
     } else{
       rcens <- NULL
     }
    thresh <- object$thresh - object$thresh[1]
     # Fit a nonparametric survival function (Turnbull, 1976)
     np <- np.ecdf(dat = dat,
                         thresh = 0,
                         rcens = rcens,
                         ltrunc = ltrunc,
                         rtrunc = rtrunc,
                         tol = 1e-12,
                         vcov = FALSE)
     # Create a weighted empirical CDF
     ecdffun <- .wecdf(x = np$x, w = np$par)
     if(is.null(rcens)){# no censoring
        xpos <- ecdffun(dat)
        if(!is.null(ltrunc) && is.null(rtrunc)){
          xpos <- (xpos - ecdffun(ltrunc))/(1-ecdffun(ltrunc))
        } else if(!is.null(rtrunc) && is.null(ltrunc)){
          xpos <- xpos/ecdffun(rtrunc)
        } else if(!is.null(rtrunc) && !is.null(ltrunc)){
          xpos <- (xpos - ecdffun(ltrunc))/(ecdffun(rtrunc) - ecdffun(ltrunc))
        }
     } else{ #right censoring
          xpos <- ecdffun(dat[!rcens])
         if(!is.null(ltrunc) && is.null(rtrunc)){
          xpos <- (xpos - ecdffun(ltrunc[!rcens]))/(1-ecdffun(ltrunc[!rcens]))
        } else if(!is.null(rtrunc) && is.null(ltrunc)){
          xpos <- xpos/ecdffun(rtrunc[!rcens])
        } else if(!is.null(rtrunc) && !is.null(ltrunc)){
          xpos <- (xpos - ecdffun(ltrunc[!rcens]))/(ecdffun(rtrunc[!rcens]) - ecdffun(ltrunc[!rcens]))
        }
     }
     scale <- as.numeric(object$par[1])
     shape <- as.numeric(object$par[-1])
     pmod <- function(q, scale, shape, family){
       switch(object$family,
                    exp = pexp(q = q, rate = 1/scale),
                    gp = pgpd(q = q, loc = 0, scale = scale, shape = shape),
                    gomp = pgomp(q = q, scale = scale, shape = shape),
                    weibull = pweibull(q = q, shape = shape, scale = scale),
                    extgp = pextgp(q = q, scale = scale, shape1 = shape[1], shape2 = shape[2]),
                    gppiece = pgppiece(q = q, scale = scale, shape = shape, thresh = thresh)
                    )
     }
     qmod <- function(p, scale, shape, family){
       switch(object$family,
              exp = qexp(p = p, rate = 1/scale),
              gp = qgpd(p = p, loc = 0, scale = scale, shape = shape),
              gomp = qgomp(p = p, scale = scale, shape = shape),
              weibull = qweibull(p = p, shape = shape, scale = scale),
              extgp = qextgp(p = p, scale = scale, shape1 = shape[1], shape2 = shape[2]),
              gppiece = qgppiece(p = p, scale = scale, shape = shape, thresh = thresh)
       )
     }
     if(!is.null(rcens)){
      ypos <- pmod(q = dat[!rcens], scale = scale, shape = shape, family = family)
      if(!is.null(ltrunc) && is.null(rtrunc)){
        ypos <- (ypos - pmod(q = ltrunc[!rcens], scale = scale, shape = shape, family = family))/(1-pmod(q = ltrunc[!rcens], scale = scale, shape = shape, family = family))
      } else if(!is.null(rtrunc) && is.null(ltrunc)){
        ypos <- ypos/pmod(q = rtrunc[!rcens], scale = scale, shape = shape, family = family)
      } else if(!is.null(rtrunc) && !is.null(ltrunc)){
        ypos <- (ypos - pmod(q = ltrunc[!rcens], scale = scale, shape = shape, family = family))/(pmod(q = rtrunc[!rcens], scale = scale, shape = shape, family = family) - pmod(q = ltrunc[!rcens], scale = scale, shape = shape, family = family))
      }
     } else{
       ypos <- pmod(q = dat, scale = scale, shape = shape, family = family)
       if(!is.null(ltrunc) && is.null(rtrunc)){
         ypos <- (ypos - pmod(q = ltrunc, scale = scale, shape = shape, family = family))/(1-pmod(q = ltrunc, scale = scale, shape = shape, family = family))
       } else if(!is.null(rtrunc) && is.null(ltrunc)){
         ypos <- ypos/pmod(q = rtrunc, scale = scale, shape = shape, family = family)
       } else if(!is.null(rtrunc) && !is.null(ltrunc)){
         ypos <- (ypos - pmod(q = ltrunc, scale = scale, shape = shape, family = family))/(pmod(q = rtrunc, scale = scale, shape = shape, family = family) - pmod(q = ltrunc, scale = scale, shape = shape, family = family))
       }
     }
     # if(any(c("qq","erp") %in% which.plot)){
     if(any(c("qq") %in% which.plot)){
       txpos <- xpos
       if(!is.null(ltrunc) && !is.null(rtrunc)){
         txpos <- xpos*pmod(q = switch(is.null(rcens), rtrunc, rtrunc[!rcens]), scale = scale, shape = shape, family = object$family) +
           (1-xpos)*pmod(q = switch(is.null(rcens), ltrunc, ltrunc[!rcens]), scale = scale, shape = shape, family = object$family)
       } else if(is.null(ltrunc) && !is.null(rtrunc)){
         txpos <- xpos*pmod(q = switch(is.null(rcens), rtrunc, rtrunc[!rcens]), scale = scale, shape = shape, family = object$family)
       } else if(is.null(rtrunc) && !is.null(ltrunc)){
         txpos <- (1-xpos)*pmod(q = switch(is.null(rcens), ltrunc, ltrunc[!rcens]), scale = scale, shape = shape, family = object$family)
       }
     }
     if(plot.type == "base"){
        if("pp" %in% which.plot){
      plot(y = ypos,
           x = xpos,
           bty = "l",
           pch = 20,
           xlab = "theoretical quantiles",
           ylab = "empirical quantiles")
        }
       if("exp" %in% which.plot){
         # exponential q-q plot
         plot(y = -log(1-ypos),
              x = -log(1-xpos),
              bty = "l",
              pch = 20,
              xlab = "theoretical quantiles",
              ylab = "empirical quantiles")
       }
       if("qq" %in% which.plot){
          plot(y = switch(is.null(rcens), dat, dat[!rcens]),
               x = qmod(p = txpos, scale = scale, shape = shape, family = object$family),
               bty = "l",
               pch = 20,
               xaxs = "i",
               yaxs = "i",
               xlab = "theoretical quantiles",
               ylab = "empirical quantiles",
               panel.first = {abline(a = 0, b = 1)})
       }
       # if("erp" %in% which.plot){
       #   np2 <- np.ecdf(dat = dat[!rcens],
       #                 thresh = 0,
       #                 ltrunc = ltrunc[!rcens],
       #                 rtrunc = rtrunc[!rcens],
       #                 tol = 1e-12,
       #                 vcov = FALSE)
       #   # Create a weighted empirical CDF
       #   ecdffun2 <- .wecdf(x = np2$x, w = np2$par)
       #   plot(y = ecdffun2(dat[!rcens]),
       #        x = ecdffun2(qmod(p = txpos, scale = scale, shape = shape, family = object$family)),
       #        bty = "l",
       #        pch = 20,
       #        xaxs = "i",
       #        yaxs = "i",
       #        xlab = "theoretical quantiles",
       #        ylab = "empirical quantiles",
       #        panel.first = {abline(a = 0, b = 1)})
       # }
       return(NULL)
     }

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
        rval <- approxfun(vals, 1-cumsum(wmatch),
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

