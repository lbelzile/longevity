#' Sample observations from a doubly truncated generalized Pareto distribution
#'
#' @param n sample size
#' @param scale scale parameter
#' @param shape shape parameter(s)
#' @param family string; choice of parametric family, either exponential (\code{exp}), Weibull (\code{weibull}), generalized Pareto (\code{gp}), Gompertz (\code{gomp}) or extended generalized Pareto (\code{extgp}).
#' @param lower vector of lower bounds
#' @param upper vector of upper bounds
#' @return a vector of \code{n} observations
rdtrunc_elife <- function(n,
                          scale,
                          shape,
                          lower,
                          upper,
                          family = c("exp","gp","gomp","weibull","extgp")
                          ){
  family <- match.arg(family)
  if(length(lower) > 1){
    stopifnot("Sample size should match length of 'lower'" = length(lower) == n)
  }
  if(length(upper) > 1){
    stopifnot("Sample size should match length of 'upper'" = length(upper) == n)
  }
  if(family == "gomp"){
    stopifnot("Scale and shape parameters must be positive" = isTRUE(scale > 0 && shape[1] >= 0))
    if(isTRUE(all.equal(shape[1], 0, ignore.attributes = TRUE))){
      family == "gp"
      shape[1] <- 0
    } else{
      return(qgomp(pgomp(lower, scale = scale, shape = shape[1]) +
                     runif(n)*(pgomp(upper, scale = scale, shape = shape[1]) - pgomp(lower, scale = scale, shape = shape[1])),
                   scale = scale, shape = shape[1])
              )
    }
  }
  if(family == "exp"){
    stopifnot("Scale parameter must be positive" = isTRUE(scale > 0))
    shape <- 0
    family <- "gp"
  }
  if(family == "gp"){
    stopifnot("Scale and shape parameters must be in range of admissible values" = isTRUE(scale > 0 && shape[1] >= -1),
              "Upper bound must be lower than the maximum value of the support." = ifelse(shape[1] < 0, max(upper) < -scale/shape[1], TRUE)
    )
    return(qgpd(pgpd(lower, scale = scale, shape = shape[1]) +
         runif(n)*(pgpd(upper, scale = scale, shape = shape[1]) - pgpd(lower, scale = scale, shape = shape[1])),
       scale = scale, shape = shape[1])
  )
  } else if(family == "extgp"){
    stopifnot("Scale and shape parameters must be in range of admissible values" = isTRUE(scale > 0 && shape[1] >= 0),
              "Shape parameters must be in the range of admissible values" = 1 - shape[1] / shape[2] > 0,
              "Upper bound must be lower than the maximum value of the support." = ifelse(shape[2] < 0, max(upper) < scale / shape[1] * log(1 - shape[1] / shape[2]), TRUE)
    )
    return(qextgp(pextgp(lower, scale = scale, shape1 = shape[1], shape2 = shape[2]) +
                  runif(n)*(pextgp(upper, scale = scale, shape1 = shape[1], shape2 = shape[2]) - pextgp(lower, scale = scale, shape1 = shape[1], shape2 = shape[2])),
                scale = scale, shape1 = shape[1], shape2 = shape[2])
    )
  } else if(family == "weibull"){
    stopifnot("Scale and shape parameters must be in range of admissible values" = isTRUE(scale > 0 && shape[1] > 0))
    return(qweibull(pweibull(lower, scale = scale, shape = shape[1]) +
                  runif(n)*(pweibull(upper, scale = scale, shape = shape[1]) - pweibull(lower, scale = scale, shape = shape[1])),
                scale = scale, shape = shape[1])
    )
  }
}
