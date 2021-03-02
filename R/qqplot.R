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
#' plot(fitted)
#' # Left- and right-truncated data
#' samp <- samp_elife(
#'  n = 200,
#'  scale = 2,
#'  shape = 0.3,
#'  family = "gp",
#'  lower = ltrunc <- runif(200),
#'  upper = rtrunc <- runif(200, 0, 10),
#'  type = "ltrc")
#' fitted <- fit_elife(
#'  dat = samp$dat,
#'  thresh = 0,
#'  ltrunc = ltrunc,
#'  rtrunc = rtrunc,
#'  type = "ltrt",
#'  family = "gp",
#'  export = TRUE)
#' plot(fitted)
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
     if(n > 2000){
       warning("Nonparametric estimation of the function very expensive for the given sample size.")
     }

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
       np2 <- np.ecdf(dat = dat[!rcens],
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
           pl_list[["qq"]] <-
             ggplot(data = data.frame(y = switch(ind_rcens, dat, dat[!rcens]),
                                      x = qmod(p = txpos, scale = scale, shape = shape, family = object$family)),
                  mapping = aes(x = x, y = y)) +
             geom_abline(intercept = 0, slope = 1, col = "gray") +
             geom_point() +
             labs(x = "theoretical quantiles",
                  y = "empirical quantiles")
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

#' @keywords internal
uq1_qqplot_elife <-
  function(dat,
           lower,
           upper,
           rcens = NULL,
           type = c("none","ltrt","ltrc"),
           family = c("exp","gp","gomp","weibull","extgp")
  ){
    # parametric bootstrap samples
    # - simulate new data with the same sampling scheme
    # - estimate parameters of the distribution
    # - compute quantiles corresponding to plotting positions
  }

#' Approximate the uncertainty
#' @keywords internal
uq2_qqplot_elife <-
  function(dat,
           lower,
           upper,
           rcens = NULL,
           type = c("none","ltrt","ltrc"),
           family = c("exp","gp","gomp","weibull","extgp")
  ){
    # parametric bootstrap samples
    # - simulate new data with the same sampling scheme
    # - estimate parameters of the distribution
    # - compute quantiles corresponding to plotting positions
  }
