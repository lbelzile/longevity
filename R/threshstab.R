#' Threshold stability plots
#'
#' The generalized Pareto and exponential distribution
#' are threshold stable. This property, which is used
#' for extrapolation purposes, can also be used to diagnose
#' goodness-of-fit: we expect the parameters \eqn{\xi} and \eqn{\tilde{\sigma} = \sigma + \xi u}
#' to be constant over a range of thresholds. The threshold stability
#' plot consists in plotting maximum likelihood estimates with pointwise confidence interval.
#' This function handles interval truncation and right-censoring.
#'
#' @details The shape estimates are constrained
#' @inheritParams nll_elife
#' @param method string giving the type of pointwise confidence interval, either Wald (\code{wald}) or profile likelihood (\code{lrt})
#' @param level probability level for the pointwise confidence intervals
#' @param plot.type string; either \code{base} for base R plots or \code{ggplot} for \code{ggplot2} plots
#' @param which.plot string; which parameters to plot;
#' @param plot logical; should a plot be returned alongside with the estimates? Default to \code{TRUE}
#' @seealso \code{\link[mev]{tstab.gpd}}, \code{\link[ismev]{gpd.fitrange}}, \code{\link[evd]{tcplot}}
#' @export
#' @return an invisible list with pointwise estimates and confidence intervals for the scale and shape parameters
tstab <- function(dat,
                  thresh,
                  ltrunc = NULL,
                  rtrunc = NULL,
                  rcens = NULL,
                  type = c("none","ltrt","ltrc"),
                  family = c('gp','exp'),
                  method = c("wald","profile"),
                  level = 0.95,
                  plot = TRUE,
                  plot.type = c("base","ggplot"),
                  which.plot = c("scale","shape")
                  ){
  family <- match.arg(family)
  method <- match.arg(method)
  if(plot){
    plot.type <- match.arg(plot.type)
    which.plot <- match.arg(which.plot, choices = c("scale","shape"), several.ok = TRUE)
  }
  wald_confint <- function(par, std.error, level = 0.95){
    c(par[1], par[1] +  qnorm(0.5+level/2)*std.error[1]*c(-1,1))
  }
  stopifnot("Provide multiple thresholds" = length(thresh) > 1)
  scale_par_mat <- matrix(NA, nrow = length(thresh), ncol = 3)
  if(family == "gp"){
    shape_par_mat <- matrix(NA, nrow = length(thresh), ncol = 3)
  }
  for(i in 1:length(thresh)){
    opt_mle <- fit_elife(dat = dat,
                           thresh = thresh[i],
                           ltrunc = ltrunc,
                           rtrunc = rtrunc,
                           rcens = rcens,
                           type = type,
                           family = family)
    if(method == "profile"){
      if(family == "gp"){
        if("shape" %in% which.plot){
          shape_i <- try(prof_gp_shape_confint(mle = opt_mle,
                                        dat = dat,
                                        thresh = thresh[i],
                                        ltrunc = ltrunc,
                                        rtrunc = rtrunc,
                                        rcens = rcens,
                                        type = type,
                                        level = level))
          if(!is.character(shape_i)){
            shape_par_mat[i,] <- shape_i
          }
        }
        if("scale" %in% which.plot){
         scalet_i <- try(prof_gp_scalet_confint(mle = opt_mle,
                                               dat = dat,
                                               thresh = thresh[i],
                                               ltrunc = ltrunc,
                                               rtrunc = rtrunc,
                                               rcens = rcens,
                                               type = type,
                                               level = level))
         if(!is.character(scalet_i)){
           scale_par_mat[i,] <- scalet_i
         }
        }
      } else if(family == "exp"){
        scale_i <- try(prof_exp_scale_confint(mle = opt_mle,
                                           dat = dat,
                                           thresh = thresh[i],
                                           ltrunc = ltrunc,
                                           rtrunc = rtrunc,
                                           rcens = rcens,
                                           type = type,
                                           level = level))
        if(!is.character(scale_i)){
          scale_par_mat[i,] <- scale_i
        }
      }
    } else if(method == "wald"){
      if(family == "gp"){
        sigma_t_mle <- opt_mle$par[1] - opt_mle$par[2]*thresh[i]
        dV <- matrix(c(1, -thresh[i]), ncol = 1)
        if(!is.null(opt_mle$vcov)){
        v <- t(dV) %*% opt_mle$vcov %*% dV
          shape_par_mat[i,] <- wald_confint(par = opt_mle$par[2], opt_mle$std.error[2], level = level)
          scale_par_mat[i,] <- wald_confint(par = sigma_t_mle, std.error = sqrt(v), level = level)
        }
      } else if(family == "exp"){
        if(!is.null(opt_mle$vcov)){
        scale_par_mat[i,] <- wald_confint(par = opt_mle$par, opt_mle$std.error)
        }
      }
    }

  }
  colnames(scale_par_mat) <- c("estimate","lower","upper")
  if(family == "gp"){
    colnames(shape_par_mat) <- c("estimate","lower","upper")
    if(!"scale" %in% which.plot){
      scale_par_mat <- NULL
    }
    if(! "shape" %in% which.plot){
      shape_par_mat <- NULL
    }
    res <- structure(list(family = family,
                          scale = scale_par_mat,
                          shape = shape_par_mat,
                          thresh = thresh),
                     class = "elife_tstab")

  } else if(family == "exp"){
    res <- structure(list(family = family,
                          scale = scale_par_mat,
                          thresh = thresh),
                     class = "elife_tstab")
  }
  if(plot){
    plot(res, plot.type = plot.type, which.plot = which.plot)
  }
  res$nexc <- as.integer(sapply(thresh, function(u){sum(dat > u, na.rm = TRUE)}))
  invisible(res)
}

#' @export
plot.elife_tstab <- function(object,
                             plot.type = c("base","ggplot"),
                             which.plot = c("scale","shape"),
                             plot = TRUE,
                             ...){
  plot.type <- match.arg(plot.type)
  which.plot <- match.arg(which.plot, choices = c("scale","shape"), several.ok = TRUE)
  if(plot.type == "ggplot" && requireNamespace("ggplot2", quietly = TRUE)){
    library(ggplot2)
    ggplot_thstab <- function(x, thresh, ylab){
      stopifnot(all.equal(colnames(x), c("estimate","lower","upper")))
      g <- ggplot(data = as.data.frame(cbind(thresh = thresh,
                                             x)),
                  aes(x = thresh, y = estimate)) +
        geom_pointrange(aes(ymin=lower, ymax=upper), size = 0.5, shape = 20) +
        labs(x = "threshold", y = ylab, main = "threshold stability plot")
      return(g)
    }
    graphs <- list()
    if(object$family == "gp"){
      if("scale" %in% which.plot){
        g1 <- ggplot_thstab(x = object$scale, thresh = object$thresh, ylab = "modified scale")
        if(length(which.plot) == 1L && plot){
            print(g1)
        }
        graphs$g1 <- g1
      }
      if("shape" %in% which.plot){
        g2 <- ggplot_thstab(x = object$shape, thresh = object$thresh, ylab = "shape")
        if(length(which.plot) == 1L && plot){
          print(g2)
        }
        graphs$g2 <- g2
      }
      if(length(which.plot) == 2L && plot){
        if(requireNamespace("patchwork", quietly = TRUE)){
          library(patchwork)
          g1 + g2
        } else{
          print(g1)
          print(g2)
        }
      }
    } else if(object$family == "exp"){
      g1 <- ggplot_thstab(x = object$scale, thresh = object$thresh, ylab = "scale")
      if(plot){
        print(g1)
      }
    graphs$g1 <- g1
    }

  } else{ # Default to base plot
    base_tstab_plot <- function(x, thresh, ylab = "scale"){
      plot(x = thresh, x[,1],
           type= "p",
           pch = 19,
           bty = "l",
           ylab = ylab,
           xlab = "threshold",
           ylim = range(x), ...)
      for(i in 1:length(thresh)){
        arrows(x0 = thresh[i],
               y0 = x[i,2],
               y1 = x[i,3],
               length = 0.02,
               angle = 90,
               code = 3,
               lwd = 2)
      }
    }
    if(object$family == "exp"){
      base_tstab_plot(x = object$scale, thresh = object$thresh)
    } else{ #generalized Pareto
      if("scale" %in% which.plot){
        base_tstab_plot(x = object$scale, thresh = object$thresh, ylab = "modified scale")
      }
      if("shape" %in% which.plot){
        base_tstab_plot(x = object$shape, thresh = object$thresh, ylab = "shape")
      }
    }
  }
}

#' Profile log likelihood for the shape parameter of the generalized Pareto distribution
#'
#' This internal function is used to produce threshold stability plots.
#'
#' @inheritParams nll_elife
#' @param mle an object of class \code{elife_par}
#' @param level level of the confidence interval
#' @keywords internal
#' @return a vector of length three with the maximum likelihood of the shape and profile-based confidence interval
prof_gp_shape_confint <-
  function(mle = NULL,
           dat,
           thresh,
           ltrunc,
           rtrunc,
           rcens,
           type = c("none","ltrt","ltrc"),
           level = 0.95){
    type <- match.arg(type)
    stopifnot("Provide a single value for the level" = length(level) == 1L,
              "Level should be a probability" = level < 1 && level > 0,
              "Provide a single threshold" = length(thresh) == 1L)
    if(!is.null(mle)){
      stopifnot("`mle` should be an object of class `elife_par` as returned by `fit_elife`" =  inherits(mle, "elife_par"))
    } else{
      mle <- fit_elife(dat = dat,
                         thresh = thresh,
                         ltrunc = ltrunc,
                         rtrunc = rtrunc,
                         rcens = rcens,
                         type = type,
                         family = "gp")
    }

  xis <- seq(-0.99, 2, length = 100L)
  mdat <- max(dat)
  dev <- sapply(xis, function(xi){
    opt <- optimize(f = function(lambda){
      nll_elife(par = c(lambda, xi),
                dat = dat,
                thresh = thresh,
                type = type,
                rcens = rcens,
                ltrunc = ltrunc,
                rtrunc = rtrunc,
                family = "gp")
      },
      interval = c(ifelse(xi < 0, mdat*abs(xi), 1e-8), 10*mdat), tol = 1e-10)
    c(-2*opt$objective, opt$minimum)
  })
  prof <- list(psi = xis, pll = -2*mle$loglik+dev[1,], maxpll = 0, mle = mle$par,
               psi.max = mle$par['shape'], std.error = sqrt(mle$vcov[2,2]))
  conf_interv(prof, level = level, print = FALSE)

}

#' Profile log likelihood for the transformed scale parameter of the generalized Pareto distribution
#'
#' This internal function is used to produce threshold stability plots.
#'
#' @inheritParams nll_elife
#' @param mle an object of class \code{elife_par}
#' @param level level of the confidence interval
#' @keywords internal
#' @return confidence interval
prof_gp_scalet_confint <-
  function(mle = NULL,
           dat,
           thresh,
           ltrunc,
           rtrunc,
           rcens,
           type = c("none","ltrt","ltrc"),
           level = 0.95){
    type <- match.arg(type)
    stopifnot("Provide a single value for the level" = length(level) == 1L,
              "Level should be a probability" = level < 1 && level > 0,
              "Provide two thresholds" = length(thresh) == 1L)
    if(!is.null(mle)){
      stopifnot("`mle` should be an object of class `elife_par` as returned by `fit_elife`" =  inherits(mle, "elife_par"))
    } else{
      mle <- fit_elife(dat = dat,
                         thresh = thresh,
                         ltrunc = ltrunc,
                         rtrunc = rtrunc,
                         rcens = rcens,
                         type = type,
                         family = "gp")
    }
    mdat <- max(dat)
    sigma_t_mle <- mle$par[1] - mle$par[2]*thresh
    dV <- matrix(c(1, -thresh), ncol = 1)
    sigma_t_se <- sqrt(as.numeric(t(dV) %*% mle$vcov %*% dV))
    psi <- sigma_t_mle + seq(-5*sigma_t_se, 10*sigma_t_se, length.out = 101)
    psi <- psi[psi>0]
    # Optimize is twice as fast
    # dev <- matrix(0, ncol = 2, nrow = length(psi))
    # i_start <- which.min(abs(psi - sigma_t_mle))
    # for(i in c(i_start:nrow(dev), (i_start-1):1)){
    #   scalet <- psi[i]
    #   opt <- Rsolnp::solnp(pars = ifelse(i >= i_start,
    #                                      ifelse(i == i_start, mle$par[2]+0.01, dev[i-1,2]),
    #                                      dev[i+1,2]),
    #                        fun = function(xi){
    #                          nll_elife(par = c(scalet + xi * thresh, xi),
    #                                    dat = dat,
    #                                    thresh = thresh,
    #                                    type = type,
    #                                    rcens = rcens,
    #                                    ltrunc = ltrunc,
    #                                    rtrunc = rtrunc,
    #                                    family = "gp")
    #                        },
    #                        ineqfun = function(xi){
    #                          xi
    #                        },
    #                        ineqLB = pmax(-1, ifelse(thresh == 0, -scalet/thresh, -scalet/mdat)),
    #                        ineqUB = 6,
    #                        control = list(trace = 0))
    #   dev[i,] <- c(-2*opt$values[length(opt$values)], opt$pars)
    # }
    dev <- t(sapply(psi, function(scalet){
      opt <- optimize(f = function(xi){
        nll_elife(par = c(scalet + xi * thresh, xi),
                  dat = dat,
                  thresh = thresh,
                  type = type,
                  rcens = rcens,
                  ltrunc = ltrunc,
                  rtrunc = rtrunc,
                  family = "gp")
      },
      interval = c(pmax(-1, ifelse(thresh == 0, -scalet/thresh, -scalet/mdat)), 3),
      tol = 1e-10)
      c(-2*opt$objective, opt$minimum)
    }))
    prof <- list(psi = psi,
                 pll = -2*mle$loglik+dev[,1],
                 maxpll = 0,
                 mle = mle$par,
                 psi.max = sigma_t_mle,
                 std.error = sigma_t_se)
    conf_interv(prof, level = level, print = FALSE)
}


#' Profile log likelihood for the scale parameter of the exponential distribution
#'
#' This internal function is used to produce threshold stability plots.
#'
#' @inheritParams nll_elife
#' @param mle an object of class \code{elife_par}
#' @param level level of the confidence interval
#' @keywords internal
#' @return a vector of length three with the maximum likelihood of the scale and profile-based confidence interval
#' @keywords internal
prof_exp_scale_confint <- function(mle = NULL,
                                   dat,
                                   thresh,
                                   ltrunc,
                                   rtrunc,
                                   rcens,
                                   type = c("none","ltrt","ltrc"),
                                   level = 0.95){
  type <- match.arg(type)
  stopifnot("Provide a single value for the level" = length(level) == 1L,
            "Level should be a probability" = level < 1 && level > 0,
            "Provide a single threshold" = length(thresh) == 1L)
  if(!is.null(mle)){
    stopifnot("`mle` should be an object of class `elife_par` as returned by `fit_elife`" =  inherits(mle, "elife_par"))
  } else{
    mle <- fit_elife(dat = dat,
                       thresh = thresh,
                       ltrunc = ltrunc,
                       rtrunc = rtrunc,
                       rcens = rcens,
                       type = type,
                       family = "exp")
  }
  psi <- mle$par + seq(pmax(-mle$par + 1e-4, -4*mle$std.error), 4*mle$std.error, length.out = 201)
  pll <- sapply(psi, function(scale){
    nll_elife(par = scale,
              dat = dat,
              thresh = thresh,
              type = type,
              rcens = rcens,
              ltrunc = ltrunc,
              rtrunc = rtrunc,
              family = "exp")
  })
  conf_interv(list(psi = psi,
                   pll = -2*(mle$loglik + pll),
                   maxpll = 0,
                   psi.max = mle$par,
                   std.error = mle$std.error,
                   mle = mle$par), level = level)
}



#' Confidence intervals for profile likelihoods
#'
#' This code is adapted from the \code{mev} package.
#' @param object a list containing information about the profile likelihood in the same format as the \code{hoa} package
#' @param level probability level of the confidence interval
#' @param prob vector of length 2 containing the bounds, by default double-sided
#' @param print logical indicating whether the intervals are printed to the console
#' @param ... additional arguments passed to the function
#' @return a table with confidence intervals.
#' @keywords internal
conf_interv <- function(object,
                        level = 0.95,
                        prob = c((1 - level)/2,  1 - (1 - level)/2),
                        print = FALSE,
                        ...){
  if (!isTRUE(all.equal(diff(prob), level, check.attributes = FALSE))) {
    warning("Incompatible arguments: `level` does not match `prob`.")
  }
  args <- list(...)
  if ("warn" %in% names(args) && is.logical(args$warn)) {
    warn <- args$warn
  }  else {
    warn <- TRUE
  }
  if (length(prob) != 2) {
    stop("`prob` must be a vector of size 2")
    prob <- sort(prob)
  }
  qulev <- qnorm(1 - prob)
  conf <- rep(0,3)
  if (is.null(object$pll) && is.null(object$r)) {
    stop("Object should contain arguments `pll` or `r` in order to compute confidence intervals.")
  }
  if (is.null(object$r)) {
    object$r <- sign(object$psi.max - object$psi) *
      sqrt(2 * (object$maxpll - object$pll))
  } else {
    object$r[is.infinite(object$r)] <- NA
  }
  if (is.null(object$normal)) {
    object$normal <- c(object$psi.max, object$std.error)
  }
    fit.r <- stats::smooth.spline(x = na.omit(cbind(object$r,
                                                    object$psi)), cv = FALSE)
    pr <- predict(fit.r, c(0, qulev))$y
    pr[1] <- object$normal[1]
    conf <- pr
  if (warn) {
    if (!any(object$r > qnorm(prob[1]))) {
      warning("Extrapolating the lower confidence interval for the profile likelihood ratio test")
    }
    if (!any(object$r < qnorm(prob[2]))) {
      warning("Extrapolating the upper confidence interval for the profile likelihood ratio test")
    }
  }
    names(conf) <- c("Estimate", "Lower CI", "Upper CI")
    conf[2] <- ifelse(conf[2] > conf[1], NA, conf[2])
    conf[3] <- ifelse(conf[3] < conf[1], NA, conf[3])
    if (print) {
      cat("Point estimate for the parameter of interest psi:\n")
      cat("Maximum likelihood          :", round(object$psi.max,
                                                 3), "\n")
      cat("\n")
      cat("Confidence intervals, levels :", prob, "\n")
      cat("Wald intervals               :", round(object$psi.max +
                                                    sort(qulev) * object$std.error, 3), "\n")
      cat("Profile likelihood           :", round(conf[2:3], 3), "\n")
    }
    return(invisible(conf))
}
