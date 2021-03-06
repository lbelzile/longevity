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
#' @param which.plot vector of string indicating the plots, among \code{pp} for probability-probability plot, \code{qq} for quantile-quantile plot, \code{erp} for empirically rescaled plot (only for censored data), \code{exp} for exponential quantile-quantile plot or \code{tmd} for Tukey's mean difference plot, which is a variant of the Q-Q plot in which we map the pair \eqn{(x,y)} is mapped to \code{((x+y)/2,y-x)} are detrended
#' @param confint logical; if \code{TRUE}, creates uncertainty diagnostic via a parametric bootstrap
#' @param plot logical; if \code{TRUE}, creates a plot. Useful for returning \code{ggplot} objects without printing the graphs
#' @examples
#' samp <- samp_elife(
#'  n = 200,
#'  scale = 2,
#'  shape = 0.3,
#'  family = "gomp",
#'  lower = 0, upper = runif(200, 0, 10),
#'  type = "ltrc")
#' fitted <- fit_elife(
#'  dat = samp$dat,
#'  thresh = 0,
#'  ltrunc = 0,
#'  rcens = samp$rcens,
#'  type = "ltrc",
#'  family = "exp",
#'  export = TRUE)
#' plot(fitted, plot.type = "ggplot")
#' # Left- and right-truncated data
#' samp <- samp_elife(
#'  n = 200,
#'  scale = 2,
#'  shape = 0.3,
#'  family = "gp",
#'  lower = ltrunc <- runif(200),
#'  upper = rtrunc <- ltrunc + runif(200, 0, 10),
#'  type = "ltrt")
#' fitted <- fit_elife(
#'  dat = samp,
#'  thresh = 0,
#'  ltrunc = ltrunc,
#'  rtrunc = rtrunc,
#'  type = "ltrt",
#'  family = "gp",
#'  export = TRUE)
#' plot(fitted, plot.type = "ggplot")
plot.elife_par <- function(object,
                           plot.type = c("base","ggplot"),
                           which.plot = c("pp","qq"),
                           confint = FALSE,
                           plot = TRUE){
  #TODO incomplete; add uncertainty
  stopifnot("Object should be of class `elife_par`" = inherits(object, what = "elife_par"))
     if(is.null(object$dat)){
          stop("`object` created using a call to `fit_elife` should include the data (`export=TRUE`).")
     }
   plot.type <- match.arg(plot.type)
   which.plot <- match.arg(which.plot, choices = c("pp","qq","erp","exp","tmd"), several.ok = TRUE)
     thresh <- object$thresh
     n <- length(object$dat)
     exc <- object$dat > thresh[1]
     dat <- object$dat[exc] - thresh[1]
     nexc <- length(dat)
     if(!is.null(object$ltrunc)){
       if(length(object$ltrunc) == 1L){
        ltrunc <- rep(pmax(0, object$ltrunc - thresh[1]), length.out = nexc)
       } else if(length(object$ltrunc) == n){
        ltrunc <- pmax(0, object$ltrunc[exc] - thresh[1])
       } else{
         stop("Invalid argument: `ltrunc` must be a scalar or a vector of the same length as `dat`")
       }
     } else{
       ltrunc <- NULL
     }
     if(!is.null(object$rtrunc)){
       if(length(object$rtrunc) == 1L){
         rtrunc <- rep(object$rtrunc - thresh[1], length.out = nexc)
       } else if(length(object$rtrunc) == n){
         rtrunc <- object$rtrunc[exc] - thresh[1]
       } else{
        stop("Invalid argument: `rtrunc` must be a scalar or a vector of the same length as `dat`")
       }
     } else{
       rtrunc <- NULL
     }
     if(!is.null(object$rcens)){
       stopifnot("`rcens` must be a vector of the same length as `dat`" = length(object$rcens) == n)
       rcens <- object$rcens[exc]
     } else{
       rcens <- NULL
     }
    thresh <- object$thresh - object$thresh[1]
     # Fit a nonparametric survival function (Turnbull, 1976)
     if(is.null(rcens)){
       np <- npsurv(time = dat,
                  type = "right",
                  ltrunc = ltrunc,
                  rtrunc = rtrunc)
     } else{
       np <- npsurv(time = dat,
                    event = 1-object$rcens,
                    type = "right",
                    ltrunc = ltrunc,
                    rtrunc = rtrunc)
     }
     # Create a weighted empirical CDF
     ecdffun <- np$cdf
     if(is.null(rcens)){# no censoring
        xpos <- length(dat) * ecdffun(dat) / (length(dat) + 1L)

        if(!is.null(ltrunc) && is.null(rtrunc)){
          xpos <- (xpos - ecdffun(ltrunc))/(1-ecdffun(ltrunc))
        } else if(!is.null(rtrunc) && is.null(ltrunc)){
          xpos <- xpos/ecdffun(rtrunc)
        } else if(!is.null(rtrunc) && !is.null(ltrunc)){
          xpos <- (xpos - ecdffun(ltrunc))/(ecdffun(rtrunc) - ecdffun(ltrunc))
        }
     } else{ #right censoring
          xpos <- length(dat[!rcens]) * ecdffun(dat[!rcens])  / (length(dat[!rcens]) + 1L)
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
     if(any(c("qq","tmd","erp") %in% which.plot)){
       txpos <- xpos
       ind_rcens <- 1L + !is.null(rcens)  # FALSE = 1, TRUE = 2
       if(!is.null(ltrunc) && !is.null(rtrunc)){
         txpos <- xpos*pmod(q = switch(ind_rcens, rtrunc, rtrunc[!rcens]), scale = scale, shape = shape, family = object$family) +
           (1-xpos)*pmod(q = switch(ind_rcens, ltrunc, ltrunc[!rcens]), scale = scale, shape = shape, family = object$family)
       } else if(is.null(ltrunc) && !is.null(rtrunc)){
         txpos <- xpos*pmod(q = switch(ind_rcens, rtrunc, rtrunc[!rcens]), scale = scale, shape = shape, family = object$family)
       } else if(is.null(rtrunc) && !is.null(ltrunc)){
         txpos <- xpos + (1-xpos)*pmod(q = switch(ind_rcens, ltrunc, ltrunc[!rcens]), scale = scale, shape = shape, family = object$family)
       }
     }
     if("erp" %in% which.plot){
       if(object$type != "ltrc"){
         warning("`erp` is only useful for censored data. Use `which.plot = pp` for specifying the equivalent plot.")
         if("pp" %in% which.plot){
           which.plot <- which.plot[which.plot == "erp"]
         } else {
           which.plot[which.plot == "erp"] <- "pp"
         }
     } else{
       np2 <- np_elife(dat = dat[!rcens],
                     thresh = 0,
                     ltrunc = ltrunc[!rcens],
                     rtrunc = rtrunc[!rcens],
                     tol = 1e-12,
                     vcov = FALSE)
       # Create a weighted empirical CDF
       ecdffun2 <- np2$ecdf
       }
     }
     if(plot.type == "ggplot"){
       if(requireNamespace("ggplot2", quietly = TRUE)){
         library(ggplot2)
       } else{
         warning("`ggplot2` package is not installed. Switching to base R plots.")
         plot.type <- "base"
       }
     }
     if(plot.type == "base"){
       for(pl in which.plot){
          if(pl == "pp"){
      plot(y = ypos,
           x = xpos,
           bty = "l",
           pch = 20,
           xlab = "theoretical quantiles",
           ylab = "empirical quantiles",
           panel.first = {abline(a = 0, b = 1, col = "gray")})
        } else if(pl == "exp"){
         # exponential q-q plot
         plot(y = -log(1-ypos),
              x = -log(1-xpos),
              bty = "l",
              pch = 20,
              xlab = "theoretical quantiles",
              ylab = "empirical quantiles",
              panel.first = {abline(a = 0, b = 1, col = "gray")})
       } else if(pl == "qq"){
          plot(y = switch(ind_rcens, dat, dat[!rcens]),
               x = qmod(p = txpos, scale = scale, shape = shape, family = object$family),
               bty = "l",
               pch = 20,
               xlab = "theoretical quantiles",
               ylab = "empirical quantiles",
               panel.first = {abline(a = 0, b = 1, col = "gray")})
       } else if(pl == "tmd"){
         yp <- switch(ind_rcens, dat, dat[!rcens])
         xp <- qmod(p = txpos, scale = scale, shape = shape, family = object$family)
         plot(y = yp - xp,
              x = (xp+yp)/2,
              bty = "l",
              pch = 20,
              xlab = "average quantile",
              ylab = "quantile difference",
              panel.first = {abline(h = 0, col = "gray")})
       } else if(pl == "erp" && !is.null(rcens)){
       plot(y = ecdffun2(dat[!rcens]),
           x = ecdffun2(qmod(p = txpos, scale = scale, shape = shape, family = object$family)),
           bty = "l",
           pch = 20,
           xlab = "theoretical quantiles",
           ylab = "empirical quantiles",
           panel.first = {abline(a = 0, b = 1, col = "gray")})
       }
       }
       return(invisible(NULL))
     } else if(plot.type == "ggplot"){
       pl_list <- list()
       for(pl in which.plot){
         if(pl == "pp"){
           pl_list[["pp"]] <-
             ggplot(data = data.frame(y = ypos, x = xpos),
                  mapping = aes(x = x, y = y)) +
             geom_abline(intercept = 0, slope = 1, col = "gray") +
             geom_point() +
             labs(x = "theoretical quantiles",
                  y = "empirical quantiles")
         } else if(pl == "exp"){
           pl_list[["exp"]] <-
             ggplot(data = data.frame(y = -log(1-ypos),
                                      x = -log(1-xpos)),
                    mapping = aes(x = x, y = y)) +
             geom_abline(intercept = 0, slope = 1, col = "gray") +
             geom_point() +
             labs(x = "theoretical quantiles",
                  y = "empirical quantiles")
         } else if(pl == "qq"){
           # if(confint && object$type != "ltrc"){
           #     # TODO fix this
           #  confint_qq <- uq1_qqplot_elife(B = 1999L,
           #                                 dat = dat,
           #                                 par = object$par,
           #                                 lower = ltrunc,
           #                                 upper = rtrunc,
           #                                 rcens = rcens,
           #                                 type = object$type,
           #                                 family = object$family)
           # }
           # if(!confint){
             pl_list[["qq"]] <-
             ggplot(data = data.frame(y = switch(ind_rcens, dat, dat[!rcens]),
                                      x = qmod(p = txpos, scale = scale, shape = shape, family = object$family)),
                  mapping = aes(x = x, y = y)) +
             geom_abline(intercept = 0, slope = 1, col = "gray") +
             geom_point() +
             labs(x = "theoretical quantiles",
                  y = "empirical quantiles")
             # } else if(confint){
             #   pl_list[["qq"]] <-
             #     ggplot(data = data.frame(y = switch(ind_rcens, dat, dat[!rcens]),
             #                              x = qmod(p = txpos, scale = scale, shape = shape, family = object$family)),
             #            mapping = aes(x = x, y = y)) +
             #     geom_abline(intercept = 0, slope = 1, col = "gray") +
             #     labs(x = "theoretical quantiles",
             #          y = "empirical quantiles") +
             #     geom_line(data = data.frame(x = confint_qq$point[1,][order(object$dat)], y = sort(object$dat)), mapping = aes(x=x,y=y), col = "grey") +
             #     geom_line(data = data.frame(x = confint_qq$point[2,], y = object$dat),  mapping = aes(x=x,y=y), col = "grey") +
             #     geom_line(data = data.frame(x = confint_qq$overall[1,], y = object$dat),  mapping = aes(x=x,y=y), lty = 2) +
             #     geom_line(data = data.frame(x = confint_qq$overall[2,], y = object$dat),  mapping = aes(x=x,y=y), lty = 2) +
             #     geom_point()
             # }
         } else if(pl == "tmd"){
           pl_list[["tmd"]] <-
             ggplot(data = data.frame(yp = switch(ind_rcens, dat, dat[!rcens]),
                                      xp = qmod(p = txpos, scale = scale, shape = shape, family = object$family)),
                                     mapping = aes(x = (xp + yp) / 2, y = yp - xp)) +
             geom_hline(yintercept = 0, col = "gray") +
             geom_point() +
             labs(x = "average quantile",
                  y = "quantile difference")

         } else if(pl == "erp" && !is.null(rcens)){
           pl_list[["erp"]] <-
             ggplot(data = data.frame(y = ecdffun2(dat[!rcens]),
                                      x = ecdffun2(qmod(p = txpos, scale = scale, shape = shape, family = object$family))),
                   mapping = aes(x = x, y = y)) +
             geom_abline(intercept = 0, slope = 1, col = "gray") +
             geom_point() +
             labs(x = "theoretical quantiles",
                  y = "empirical quantiles")
         }
       }
       if(plot){
         lapply(pl_list, print)
       }
       return(invisible(pl_list))
     }
}

#' Uncertainty quantification for quantile-quantile plots
#' @export
#' @param B number of bootstrap samples
#' @param dat vector of data
#' @param par parameter estimates of the model
#' @param lower lower bounds (truncation or lowest possible value)
#' @param upper upper bound for right-censoring or right-truncation
#' @param level level of the confidence intervals
#' @inheritParams nll_elife
#' @keywords internal
uq1_qqplot_elife <-
  function(B = 9999L,
           dat,
           par,
           lower,
           upper,
           rcens = NULL,
           level = 0.95,
           type = c("none","ltrt","ltrc"),
           family = c("exp","gp","gomp","weibull","extgp")
  ){
    cens <- !is.null(rcens)
    type <- match.arg(type)
    family <- match.arg(family)
    if(missing(lower)){
      ltrunc <- 0
    }
    if(missing(upper)){
      rtrunc <- Inf
    }
    if(!missing(lower) && !missing(upper)){
      if(length(upper) != 1 && length(upper) != 1){
      stopifnot( "`upper` and `lower` should be vectors of the same length." = length(lower) == length(upper),
                 "`Length of data `dat` does not match vector of lower and upper bounds." = length(dat) == length(upper))
      }
    }
    stopifnot("Invalid `rcens` argument." = (is.null(rcens) && type != "ltrc") || (is.logical(rcens) && type == "ltrc"),
              "`lower` should be smaller than `upper`." = isTRUE(all(lower <= upper)),
              "Number of bootstrap samples must be larger than what is prescribed by the level." = B >= 1/(1-level) - 1L)
    # parametric bootstrap samples
    # - simulate new data with the same sampling scheme
    # - estimate parameters of the distribution
    # - compute quantiles corresponding to plotting positions
    n <- ifelse(cens, length(dat[!rcens]), length(dat))
    ppos <- matrix(NA, nrow = B, ncol = n)
    if(cens){
      ddat <- dat[!rcens]
      ltrunc <- lower[!rcens]
      rtrunc <- upper[!rcens]
    }
    if(type == "none"){
      for(b in seq_len(B)){
        dat_boot <- samp_elife(n = n,
                               scale = par[1],
                               shape = par[-1],
                               family = family,
                               type = type)
        fit_boot <- fit_elife(
          dat = dat_boot,
          thresh = 0,
          family = family,
          type = type)
        ppos[b,] <- qelife(p = ppoints(n),
                           scale = fit_boot$par[1],
                           shape = fit_boot$par[-1],
                           family = family)

      }
    } else if(type == "ltrc"){
      for(b in seq_len(B)){
      boot_dat <- samp_elife(n = n,
                             scale = par[1],
                             shape = par[-1],
                             lower = lower,
                             upper = upper,
                             family = family,
                             type = type)
       fit_boot <- fit_elife(
        dat = boot_dat$dat,
        thresh = 0,
        ltrunc = lower,
        rcens = boot_dat$rcens,
        family = family,
        type = type)
       np_boot <- npsurv(time = boot_dat$dat,
                event = !boot_dat$rcens,
                type = "right",
                ltrunc = lower)
       ppos[b,] <- qelife(p = n/(n+1)*np_boot$cdf(dat),
                          scale = fit_boot$par[1],
                          shape = fit_boot$par[-1],
                          family = family)
      }
      } else if(type == "ltrt"){
        for(b in seq_len(B)){
          boot_dat <- samp_elife(n = n,
                                 scale = par[1],
                                 shape = par[-1],
                                 lower = lower,
                                 upper = upper,
                                 family = family,
                                 type = type)
          fit_boot <- fit_elife(
            dat = boot_dat,
            thresh = 0,
            ltrunc = lower,
            rtrunc = upper,
            family = family,
            type = type)
          np_boot <- npsurv(time = boot_dat,
                            type = "interval",
                            event = rep(1L, n),
                            ltrunc = lower,
                            rtrunc = upper)$cdf
          xpos <- n / (n + 1) * (np_boot(dat) - np_boot(ltrunc))/(np_boot(rtrunc) - np_boot(ltrunc))
          ppos[b,] <- qelife(p = pelife(q = ltrunc,
                             scale = fit_boot$par[1],
                             shape = fit_boot$par[-1],
                             family = family)*(1-xpos) +
            xpos * pelife(q = rtrunc,
                         scale = fit_boot$par[1],
                         shape = fit_boot$par[-1],
                         family = family),
            scale = fit_boot$par[1],
            shape = fit_boot$par[-1],
            family = family)
        }
      }
    return(ppos)
    # return(boot::envelope(mat = ppos, level = level))
  }

#' Approximate uncertainty for diagnostic plots
#'
#' Approximate the uncertainty by resampling parameters
#' estimates using a normal approximation to the maximum
#' likelihood estimators.
#'
#' @param B integer; number of simulations
#' @param object an object of class \code{elife_par}
#' @param logscale logical; if \code{TRUE} (default), compute the approximation using \eqn{\log(\sigma)} rather than \eqn{\sigma}
#' @keywords internal
# uq2_qqplot_elife <- function(B = 1e4,
#                              object,
#                              logscale = TRUE){
#   stopifnot("Cannot compute approximation to sampling distribution of maximum likelihood estimators: missing `vcov` or `par` values" = !is.null(object$par) && !is.null(object$vcov))
#   par <- object$par
#   vcov <- object$vcov
#
# }
