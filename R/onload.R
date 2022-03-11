
.s3_register <-
function (generic, class, method = NULL){
  stopifnot(is.character(generic), length(generic) == 1)
  stopifnot(is.character(class), length(class) == 1)
  pieces <- strsplit(generic, "::")[[1]]
  stopifnot(length(pieces) == 2)
  package <- pieces[[1]]
  generic <- pieces[[2]]
  caller <- parent.frame()
  get_method_env <- function() {
    top <- topenv(caller)
    if (isNamespace(top)) {
      asNamespace(environmentName(top))
    }
    else {
      caller
    }
  }
  get_method <- function(method, env) {
    if (is.null(method)) {
      get(paste0(generic, ".", class), envir = get_method_env())
    }
    else {
      method
    }
  }
  register <- function(...) {
    envir <- asNamespace(package)
    method_fn <- get_method(method)
    stopifnot(is.function(method_fn))
    if (exists(generic, envir)) {
      registerS3method(generic, class, method_fn, envir = envir)
    }
    else if (identical(Sys.getenv("NOT_CRAN"), "true")) {
      warning(sprintf("Can't find generic `%s` in package %s to register S3 method.",
                      generic, package))
    }
  }
  setHook(packageEvent(package, "onLoad"), register)
  if (isNamespaceLoaded(package)) {
    register()
  }
  invisible()
}

.onLoad <- function(...) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    .s3_register("ggplot2::autoplot", "elife_par")
    .s3_register("ggplot2::autoplot", "elife_northropcoleman")
    .s3_register("ggplot2::autoplot", "elife_hazard")
    .s3_register("ggplot2::autoplot", "elife_tstab")
  }
}


#' @importFrom graphics abline
#' @importFrom graphics arrows
#' @importFrom graphics hist
#' @importFrom graphics mtext
#' @importFrom graphics rug
#' @importFrom stats approxfun
#' @importFrom stats dexp
#' @importFrom stats dweibull
#' @importFrom stats family
#' @importFrom stats na.omit
#' @importFrom stats optim
#' @importFrom stats optimize
#' @importFrom stats pchisq
#' @importFrom stats pexp
#' @importFrom stats ppoints
#' @importFrom stats predict
#' @importFrom stats pweibull
#' @importFrom stats qchisq
#' @importFrom stats qexp
#' @importFrom stats qnorm
#' @importFrom stats qweibull
#' @importFrom stats rexp
#' @importFrom stats rmultinom
#' @importFrom stats runif
#' @importFrom stats rweibull
#' @importFrom stats smooth.spline
#' @importFrom stats weighted.mean
NULL
