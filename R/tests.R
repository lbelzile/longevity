#' Likelihood ratio test for covariates
#'
#' This function fits separate models for each distinct
#' value of covariate and computes a likelihood ratio test
#' to test whether there are significante differences between
#' groups.
#'
#' @export
#' @inheritParams nll_elife
#' @param covariate vector of factors, logical or integer whose distinct values are
#' @return a list with elements
#' \itemize{
#' \item{\code{stat}: }{likelihood ratio statistic}
#' \item{\code{df}: }{degrees of freedom}
#' \item{\code{pval}: }{the p-value obtained from the asymptotic chi-square approximation.}
#' }
#' @examples
#' with(uk110,
#' test_elife(dat = ndays,
#' thresh = 40178L, covariate = gender,
#' ltrunc = ltrunc, rtrunc = rtrunc,
#' family = "exp", type = "ltrt"))
test_elife <- function(time,
                       time2 = NULL,
                       event = NULL,
                       covariate,
                       thresh = 0,
                       ltrunc = NULL,
                       rtrunc = NULL,
                       type = c("right", "left","interval","interval2"),
                       family = c("exp", "gp", "weibull", "gomp", "gompmake", "extgp"),
                       weights = rep(1, length(time))) {
  family <- match.arg(family)
  type <- match.arg(type)
  stopifnot("Covariate must be provided" = !missing(covariate),
            "Object `covariate` should be of the same length as `dat`" = length(covariate) == length(time),
            "Provide a single threshold" = !missing(thresh) && length(thresh) == 1L)
  npar <- switch(family,
               "exp" = 1L,
               "gp" = 2L,
               "gomp" = 2L,
               "gompmake" = 3L,
               "extgp" = 3L,
               "weibull" = 2L)
  survout <- .check_surv(time = time,
                         time2 = time2,
                         event = event,
                         type = type)
  time <- survout$time
  time2 <- survout$time2
  status <- survout$status
  # Transform to factor
  covariate <- as.factor(covariate)
  nobs_cov <- table(covariate)
  m <- length(nobs_cov)
  stopifnot("There should be more than one group in `covariate`." = m > 1,
            "There are too few observations (less than 5 times the number of parameters) for some modalities of `covariate`." = min(nobs_cov) >= 5*npar)
  # Fit the pooled model
  labels <- names(nobs_cov)
  fit_null <- try(fit_elife(time = time,
                            time2 = time2,
                            event = event,
                            status = status,
                            thresh = thresh,
                            ltrunc = ltrunc,
                            rtrunc = rtrunc,
                            type = type,
                            family = family,
                            weights = weights))
  loglik0 <- ifelse(is.character(fit_null), NA, fit_null$loglik)
  fit_alternative <- list()
  loglik1 <- rep(0, m)
  n_levels <- rep(0L, m)
  for(i in 1:m){
    fit_alternative[[i]] <- try(fit_elife(time = time[covariate == labels[i]],
                                          time2 = time2[covariate == labels[i]],
                                          status = status[covariate == labels[i]],
                                thresh = thresh,
                                ltrunc = ltrunc[covariate == labels[i]],
                                rtrunc = rtrunc[covariate == labels[i]],
                                type = type,
                                family = family,
                                weights = weights[covariate == labels[i]]))
    loglik1[i] <- ifelse(is.character(fit_alternative[[i]]), NA, fit_alternative[[i]]$loglik)
    n_levels[i] <- fit_alternative[[i]]$nexc
  }
  lrt_stat <- 2*as.numeric((sum(loglik1)-loglik0))
  names(n_levels) <- labels
  p_value <- pchisq(q = lrt_stat, df = (m - 1) * npar, lower.tail = FALSE)
    invisible(structure(list(stat = lrt_stat,
              df = (m - 1) * npar,
              pval = p_value,
              nobs_covariate = n_levels,
              thresh = thresh,
              family = family),
    class = "elife_par_test"))
}
#' @export
print.elife_par_test <-   function(x,
             digits = min(3, getOption("digits")),
             na.print = "", ...){
      cat("Model:", switch(x$family,
                           exp = "exponential",
                           gomp = "Gompertz",
                           gompmake = "Gompertz-Makeham",
                           weibull = "Weibull",
                           extgp = "extended generalized Pareto",
                           gp = "generalized Pareto",
                           gppiece = "piecewise generalized Pareto"),
          "distribution.", "\n")
      cat("Threshold:", round(x$thresh, digits), "\n")
      cat("Number of exceedances per covariate level:\n")
      print.default(x$nobs_covariate)
      cat("\nLikelihood ratio statistic:", format(x$stat, digits = digits))
      cat(paste0("\nNull distribution: chi-square (", x$df, ")\n"))
      cat("Asymptotic p-value: ", format(x$pval, digits = digits),"\n")
      invisible(x)
}

#' @export
anova.elife_par <- function(object,
                            object2,
                            ...,
                            test = c("Chisq","bootstrap")){
  if (any(missing(object), missing(object2))){
    stop("Two models must be specified.")
  }
  test <- match.arg(test)
  model1 <- deparse(substitute(object))
  model2 <- deparse(substitute(object2))
  models <- c(model1, model2)
  narg <- 2L
  for (i in 1:narg) {
    if (!inherits(get(models[i], envir = parent.frame()), "elife_par")){
      stop("Invalid input: use only with objects of class 'elife_par'.")
    }
  }

  npar <- rep(0, length(models))
  dev <- rep(0, length(models))
  thresh <- nobs <- rep(0, length(models))
  family <- rep("", length(models))
  conv <- rep(FALSE, 2)
  for (i in 1:narg) {
    elifemod <- get(models[i], envir = parent.frame())
    dev[i] <- 2*elifemod$loglik
    npar[i] <- length(elifemod$par)
    thresh[i] <- elifemod$thresh[1]
    nobs[i] <- elifemod$nexc
    conv[i] <- elifemod$convergence
    family[i] <- elifemod$family
  }
  if(npar[1] < npar[2]){
    dev <- dev[2:1]
    npar <- npar[2:1]
    family <- family[2:1]
  }
  if(!isTRUE(all.equal(thresh[1], thresh[2]))){
    stop("Invalid arguments: the thresholds should be the same.")
  }
  if(!isTRUE(all(conv))){
    stop("At least one of the optimization failed to converge.")
  }
  if(!isTRUE(all.equal(nobs[1], nobs[2]))){
    stop("Invalid arguments: the observations should be the same.")
  }
  if(family[1] == family[2]){
    stop("Models should be of different families.")
  }
  # Cases considered
  nmods <- rbind(
    c("exp","weibull","regular"),
    c("exp","gp","regular"),
    c("exp","gppiece","regular"),
    c("gomp","extgp","regular"),
    c("gp","gppiece","regular"),
    c("gomp","gompmake","regular"),
    c("exp","gompmake","boundary"),
    c("exp","gomp","boundary"),
    c("exp","extgp","boundary"),
    c("gp","extgp","boundary")
  )
  match_family <- which(apply(nmods[,1:2], 1, function(fam){isTRUE(all(family %in% fam))}))
  stopifnot("Invalid input: models are not nested" = length(match_family) == 1L)

  df <- -diff(npar)
  dvdiff <- -diff(dev)
  if(dvdiff < 0 && dvdiff > -1e-4){
    # Numerical tolerance for zero
    dvdiff <- 0
  }
  if(dvdiff < 0){
    stop("The alternative model has a lower likelihood value than the null model, indicating convergence problems.")
  }
  if(test == "Chisq"){
    if(nmods[match_family,3] == "regular"){ #regular model
     pval <- pchisq(dvdiff, df = df, lower.tail = FALSE)
    } else if(nmods[match_family,3] == "boundary"){
      pval <- 0.5*pchisq(dvdiff, df = df, lower.tail = FALSE) + 0.5*pchisq(dvdiff, df = df - 1L, lower.tail = FALSE)
    }
  }
  table <- data.frame(npar, dev, c(NA, df), c(NA, dvdiff), c(NA,pval))
  dimnames(table) <- list(models, c("npar", "Deviance", "Df",
                                    "Chisq", "Pr(>Chisq)"))
  rownames(table) <- family
  structure(table, heading = c("Analysis of Deviance Table\n"),
            class = c("anova", "data.frame"))
}

#' Score test of Northrop and Coleman
#'
#' This function computes the score test
#' with the piecewise generalized Pareto distribution
#' under the null hypothesis that the generalized Pareto
#' with a single shape parameter is an adequate simplification.
#' The test statistic is calculated using the observed information
#' matrix; both hessian and score vector are obtained through numerical differentiation.
#'
#' The reference distribution is chi-square
#' @inheritParams fit_elife
#' @export
#' @param thresh a vector of thresholds
#' @return a data frame with the following variables:
#' \itemize{
#' \item{\code{thresh}: }{threshold for the generalized Pareto}
#' \item{\code{nexc}: }{number of exceedances}
#' \item{\code{score}: }{score statistic}
#' \item{\code{df}: }{degrees of freedom}
#' \item{\code{pval}: }{the p-value obtained from the asymptotic chi-square approximation.}
#' }
nc_score_test <- function(time,
                          time2 = NULL,
                          event = NULL,
                          thresh = 0,
                          ltrunc = NULL,
                          rtrunc = NULL,
                          type = c("right", "left", "interval", "interval2"),
                          weights = rep(1, length(time))){
  stopifnot("Threshold is missing" = !missing(thresh))
  nt <- length(thresh)
  thresh <- sort(unique(thresh))
  stopifnot("Threshold should be at least length two" = nt >= 2L)
  res <- as.data.frame(matrix(NA, ncol = 5, nrow = nt - 1L))
  colnames(res) <- c("thresh","nexc","score","df","pval")
  for(i in 1:(nt-1L)){
    fit0 <- fit_elife(time = time,
                      time2 = time2,
                      event = event,
                      thresh = thresh[i],
                      ltrunc = ltrunc,
                      rtrunc = rtrunc,
                      type = type,
                      family = "gp",
                      weights = weights,
                      export = FALSE)
    score0 <- try(numDeriv::grad(func = function(x){
                     nll_elife(par = x,
                               time = time,
                               time2 = time2,
                               event = event,
                               thresh = thresh[i:nt],
                               type = type,
                               ltrunc = ltrunc,
                               rtrunc = rtrunc,
                               family = "gppiece",
                               weights = weights
                     )},
                     x = c(fit0$par['scale'], rep(fit0$par['shape'], nt-i+1L))
                   ))
    hess0 <- try(numDeriv::hessian(func = function(x){
                    nll_elife(par = x,
                              time = time,
                              time2 = time2,
                              event = event,
                              thresh = thresh[i:nt],
                              type = type,
                              ltrunc = ltrunc,
                              rtrunc = rtrunc,
                              family = "gppiece",
                              weights = weights
                    )},
                    x = c(fit0$par['scale'],rep(fit0$par['shape'], nt-i+1L))
                  ))
    score_stat <- try(as.numeric(score0 %*% solve(hess0) %*% score0))
    if(!is.character(score_stat)){
      res$score[i] <- score_stat
      res$df[i] <- nt - i
      res$pval[i] <- pchisq(q = score_stat, df = nt - i, lower.tail = FALSE)
    }
    res$nexc[i] <- fit0$nexc
  }
  res$thresh <- thresh[-nt]
  return(res)
}


#' Goodness-of-fit diagnostics
#'
#' Warning: EXPERIMENTAL
#' Compute the Kolmogorov-Smirnov
#' test statistic and compare it with a simulated null
#' distribution obtained via a parametric bootstrap.
#'
#' @note The bootstrap scheme requires simulating new data,
#' fitting a parametric model and estimating the nonparametric
#' maximum likelihood estimate for each new sample.
#' This is computationally intensive in large samples.
#'
#' @inheritParams fit_elife
#' @param B number of bootstrap simulations
#' @return a list with elements
#' \itemize{
#' \item{\code{stat}: }{the value of the test statistic}
#' \item{\code{pval}: }{p-value obtained via simulation}
#' }
#' @export
ks_test <- function(time,
                    time2 = NULL,
                    event = NULL,
                    thresh = 0,
                    ltrunc = NULL,
                    rtrunc = NULL,
                    type = c("right", "left","interval","interval2"),
                    family = c("exp", "gp", "weibull", "gomp", "extgp","gppiece"),
                    B = 999L){
  family <- match.arg(family)
  type <- match.arg(type)
  n <- length(time)
  survout <- .check_surv(time = time,
                         time2 = time2,
                         event = event,
                         type = type)
  time <- survout$time
  time2 <- survout$time2
  status <- survout$status
  if(thresh[1] > 0){
    # Keep only exceedances, shift observations
    # We discard left truncated observations and interval censored
    # if we are unsure whether there is an exceedance
    ind <- ifelse(status == 2, FALSE, time > thresh[1])
    weights <- weights[ind] # in this order
    time <- time[ind] - thresh[1]
    time2 <- time2[ind] - thresh[1]
    status <- status[ind]
    if(!is.null(ltrunc)){ #both ltrc and ltrt
      stopifnot("`ltrunc` must be of the same length as `time`." = length(ltrunc) == n)
      ltrunc <- pmax(0, ltrunc[ind] - thresh[1])
    }
    if(!is.null(rtrunc)){
      stopifnot("`rtrunc` must be of the same length as `time`." = length(rtrunc) == n)
      rtrunc <- rtrunc[ind] - thresh[1]
    }
  }
  thresh <- 0
  # Fit parametric model
  F0 <- try(fit_elife(time = time,
                      time2 = time2,
                      status = status,
                      thresh = thresh,
                      ltrunc = ltrunc,
                      rtrunc = rtrunc,
                      type = type,
                      family = family
                  ))
  if(is.character(F0) || !F0$convergence){
    stop("Could not estimate the parametric model.")
  }
  # Fit NPMLE of ECDF
  Fn <- npsurv(time = dat,
               type = "interval",
               event = status,
               ltrunc = ltrunc,
               rtrunc = rtrunc)
 # Compute test statistic
  if(family == "gompmake"){
    scale <- c(F0$par[1], F0$par[3])
    shape <- F0$par[2]
  } else{
    scale <- F0$par[1]
    shape <- F0$par[-1]
  }
 ks <- max(abs(Fn$cdf(Fn$x) - pelife(q = Fn$x, scale = scale, shape = shape, family = family)))
 stat <- rep(NA, B + 1L)
 stat[B + 1] <- ks
 if(is.null(ltrunc) && is.null(rtrunc) && isTRUE(all(status == 1L))){
   type2 <- "none"
 } else if(isTRUE(all(status == 1L)) && (!is.null(rtrunc) || !is.null(ltrunc))){
   type2 <- "ltrt"
 } else if(is.null(rtrunc) && isTRUE(all(status == c(0L,1L)))){
   type2 <- "ltrc"
 }
 for(b in 1:B){
   bootconv <- FALSE
   while(!bootconv){
 if(type2 == "ltrt"){
   boottime <- samp_elife(n = length(time),
                  scale = scale,
                  shape = shape,
                  lower = ltrunc,
                  upper = rtrunc,
                  family = family,
                  type2 = type2)
   bootstatus <- rep(1L, length(time))
 } else if(type2 == "none"){
  boottime <- relife(n = length(time),
                 scale = scale,
                 shape = shape,
                 family = family)
  bootstatus <- rep(1L, length(time))
 } else if(type2 == "ltrc"){
   boot_sim <- samp_elife(n = length(time),
                          scale = scale,
                          shape = shape,
                          lower = ltrunc,
                          upper = rtrunc,
                          family = family,
                          type2 = type2)
   boottime <- boot_sim$dat
   boottime2 <- ifelse(boot_sim$rcens, NA, boot_sim$dat)
   bootstatus <- as.integer(!boot_sim$rcens)
 }
   # bootstrap loop
   F0_b <- try(fit_elife(time = boottime,
                         time2 = boottime2,
                         status = bootstatus,
                         thresh = thresh,
                         ltrunc = ltrunc,
                         rtrunc = rtrunc,
                         type = "right",
                         family = family))
   # Fit NPMLE of ECDF
   if(type2 == "ltrc"){
   Fn_b <- try(npsurv(time = boottime,
                      event = bootstatus,
                      ltrunc = ltrunc,
                      rtrunc = rtrunc,
                      type = "interval"))
   } else{
     Fn_b <- try(npsurv(time = boottime,
                        ltrunc = ltrunc,
                        rtrunc = rtrunc))
   }
   if(family == "gompmake"){
     scale_boot <- c(F0_b$par[1], F0_b$par[3])
     shape_boot <- F0_b$par[2]
   } else{
     scale_boot <- F0_b$par[1]
     shape_boot <- F0_b$par[-1]
   }
   # Compute test statistic
   ks_boot <- try(max(abs(Fn_b$cdf(Fn_b$x) - pelife(q = Fn_b$x, scale = scale_boot, shape = shape_boot, family = family))))
   if(is.numeric(ks_boot)){
     stat[b] <- ks_boot
     bootconv <- TRUE
   }
  }
 }
  list(stat = ks,
       # null = stat,
       pval = mean(stat >= ks, na.rm = TRUE))
}
