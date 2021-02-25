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
#' \item{\code{lrt_stat}: }{likelihood ratio statistic}
#' \item{\code{df}: }{degrees of freedom}
#' \item{\code{p_value}: }{the p-value obtained from the asymptotic chi-square approximation.}
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
  return(list(lrt_stat = lrt_stat,
              df = (m - 1) * npar,
              p_value = p_value))
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

