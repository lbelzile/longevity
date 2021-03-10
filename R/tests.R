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
test_elife <- function(dat,
                       covariate,
                       thresh,
                       ltrunc = NULL,
                       rtrunc = NULL,
                       rcens = NULL,
                       type = c("none", "ltrc", "ltrt"),
                       family = c("exp", "gp", "weibull", "gomp", "extgp"),
                       weights = rep(1, length(dat))) {
  family <- match.arg(family)
  type <- match.arg(type)
  stopifnot("Covariate must be provided" = !missing(covariate),
            "Object `covariate` should be of the same length as `dat`" = length(covariate) == length(dat),
            "Provide a single threshold" = !missing(thresh) && length(thresh) == 1L)
  npar <- switch(family,
               "exp" = 1L,
               "gp" = 2L,
               "gomp" = 2L,
               "extgp" = 3L,
               "weibull" = 2L)
  # Transform to factor
  covariate <- as.factor(covariate)
  nobs_cov <- table(covariate)
  m <- length(nobs_cov)
  stopifnot("There should be more than one group in `covariate`." = m > 1,
            "There are too few observations (less than 5 times the number of parameters) for some modalities of `covariate`." = min(nobs_cov) >= 5*npar)
  # Fit the pooled model
  labels <- names(nobs_cov)
  if(is.null(weights)){
    weights <- rep(1, length(dat))
  }
  fit_null <- try(fit_elife(dat = dat,
                          thresh = thresh,
                          ltrunc = ltrunc,
                          rtrunc = rtrunc,
                          rcens = rcens,
                          type = type,
                          family = family,
                          weights = weights))
  loglik0 <- ifelse(is.character(fit_null), NA, fit_null$loglik)
  fit_alternative <- list()
  loglik1 <- rep(0, m)
  for(i in 1:m){
    fit_alternative[[i]] <- try(fit_elife(dat = dat[covariate == labels[i]],
                                thresh = thresh,
                                ltrunc = ltrunc[covariate == labels[i]],
                                rtrunc = rtrunc[covariate == labels[i]],
                                rcens = rcens[covariate == labels[i]],
                                type = type,
                                family = family,
                                weights = weights[covariate == labels[i]]))
    loglik1[i] <- ifelse(is.character(fit_alternative[[i]]), NA, fit_alternative[[i]]$loglik)
  }
  lrt_stat <- 2*as.numeric((sum(loglik1)-loglik0))
  p_value <- pchisq(q = lrt_stat, df = (m - 1) * npar, lower.tail = FALSE)
  return(list(stat = lrt_stat,
              df = (m - 1) * npar,
              pval = p_value))
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
  for (i in 1:narg) {
    elifemod <- get(models[i], envir = parent.frame())
    dev[i] <- 2*elifemod$loglik
    npar[i] <- length(elifemod$par)
    thresh[i] <- elifemod$thresh[1]
    nobs[i] <- elifemod$nexc
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
    c("exp","gomp","boundary"),
    c("exp","extgp","boundary"),
    c("gp","extgp","boundary")
  )
  match_family <- which(apply(nmods[,1:2], 1, function(fam){isTRUE(all(family %in% fam))}))
  stopifnot("Invalid input: models are not nested" = length(match_family) == 1L)

  df <- -diff(npar)
  dvdiff <- -diff(dev)
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
nc_score_test <- function(dat,
                          thresh,
                          ltrunc = NULL,
                          rtrunc = NULL,
                          rcens = NULL,
                          type = c("none", "ltrc", "ltrt"),
                          weights = rep(1, length(dat))){
  stopifnot("Threshold is missing" = !missing(thresh))
  nt <- length(thresh)
  thresh <- sort(unique(thresh))
  stopifnot("Threshold should be at least length two" = nt >= 2L)
  res <- as.data.frame(matrix(NA, ncol = 5, nrow = nt - 1L))
  colnames(res) <- c("thresh","nexc","score","df","pval")
  for(i in 1:(nt-1L)){
    fit0 <- fit_elife(dat = dat,
                      thresh = thresh[i],
                      ltrunc = ltrunc,
                      rtrunc = rtrunc,
                      rcens = rcens,
                      type = type,
                      family = "gp",
                      weights = weights,
                      export = FALSE)
    score0 <- try(numDeriv::grad(func = function(x){
                     nll_elife(par = x,
                               dat = dat,
                               thresh = thresh[i:nt],
                               type = type,
                               rcens = rcens,
                               ltrunc = ltrunc,
                               rtrunc = rtrunc,
                               family = "gppiece",
                               weights = weights
                     )},
                     x = c(fit0$par['scale'], rep(fit0$par['shape'], nt-i+1L))
                   ))
    hess0 <- try(numDeriv::hessian(func = function(x){
                    nll_elife(par = x,
                              dat = dat,
                              thresh = thresh[i:nt],
                              type = type,
                              ltrunc = ltrunc,
                              rtrunc = rtrunc,
                              rcens = rcens,
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
#' Compute the Kolmogorov-Smirnov or the Anderson-Darling
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
ks_test <- function(dat,
                     thresh,
                     ltrunc = NULL,
                     rtrunc = NULL,
                     rcens = NULL,
                     type = c("none", "ltrc", "ltrt"),
                     family = c("exp", "gp", "weibull", "gomp", "extgp","gppiece"),
                     B = 999L){
  family <- match.arg(family)
  type <- match.arg(type)
  ntot <- length(dat)
  wexc <- dat > thresh[1]
  dat <- dat[wexc] - thresh[1]
  n <- length(dat)
  if(!is.null(ltrunc)){
    ltrunc <- pmax(0, ltrunc[wexc] - thresh[1])
  }
  if(!is.null(rtrunc)){
    rtrunc <- rtrunc[wexc] - thresh[1]
  }
  if(type == "ltrc"){
    stopifnot("Right-censoring indicator is lacking." = !is.null(rcens) && length(rcens) == ntot)
    rcens <- rcens[wexc]
    status <- as.integer(!rcens)
  } else{
    rcens <- NULL
    status <- rep(1L, n)
  }
  thresh <- 0
  # Fit parametric model
  F0 <- try(fit_elife(dat = dat,
                  thresh = thresh,
                  rcens = rcens,
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
 ks <- max(abs(Fn$cdf(Fn$x) - pelife(q = Fn$x, scale = F0$par[1], shape = F0$par[-1], family = family)))
 stat <- rep(NA, B + 1L)
 stat[B + 1] <- ks
 for(b in 1:B){
   bootconv <- FALSE
   while(!bootconv){
 if(type == "ltrt"){
   bootsamp <- samp_elife(n = length(dat),
                  scale = F0$par[1],
                  shape = F0$par[-1],
                  lower = ltrunc,
                  upper = rtrunc,
                  family = family,
                  type = type)
   bootrcens <- NULL
 } else if(type == "none"){
  bootsamp <- relife(n = length(dat),
                 scale = F0$par[1],
                 shape = F0$par[-1],
                 family = family)
  bootrcens <- NULL
 } else if(type == "ltrc"){
   boot_sim <- samp_elife(n = length(dat),
                          scale = F0$par[1],
                          shape = F0$par[-1],
                          lower = ltrunc,
                          upper = rtrunc,
                          family = family,
                          type = type)
   bootsamp <- boot_sim$dat
   bootrcens <- boot_sim$rcens
 }
   # bootstrap loop
   F0_b <- try(fit_elife(dat = bootsamp,
                   thresh = thresh,
                   rcens = bootrcens,
                   ltrunc = ltrunc,
                   rtrunc = rtrunc,
                   type = type,
                   family = family))
   # Fit NPMLE of ECDF
   if(type == "ltrc"){
   Fn_b <- try(npsurv(time = bootsamp,
                      event = as.integer(!bootrcens),
                      ltrunc = ltrunc,
                      rtrunc = rtrunc,
                      type = "interval"))
   } else{
     Fn_b <- try(npsurv(time = bootsamp,
                        ltrunc = ltrunc,
                        rtrunc = rtrunc))
   }

   # Compute test statistic
   ks_boot <- try(max(abs(Fn_b$cdf(Fn_b$x) - pelife(q = Fn_b$x, scale = F0_b$par[1], shape = F0_b$par[-1], family = family))))
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
