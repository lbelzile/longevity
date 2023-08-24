# Check that plots work with both base R and with ggplot2

library(tinytest)
library(longevity)

set.seed(1234)
n <- 100L
x <- samp_elife(n = n,
           scale = 2,
           shape = 0.5,
           lower = low <- runif(n),
           upper = upp <- runif(n, min = 3, max = 10),
           type2 = "ltrt",
           family = "weibull")
fit <- fit_elife(
  time = x,
  ltrunc = low,
  rtrunc = upp,
  family = "weibull",
  export = TRUE)

# x <- relife(n = n,
#             scale = 2,
#             shape = 0.5,
#             family = "weibull")
# fit <- fit_elife(time = x, family = "weibull", export = TRUE)

print(fit)
plot(fit, which.plot = "pp", plot.type = "base")
plot(fit, which.plot = "qq", plot.type = "base")
# Warning because erp == pp w/o censoring
expect_warning(plot(fit, which.plot = "erp", plot.type = "base"))
plot(fit, which.plot = "exp", plot.type = "base")
plot(fit, which.plot = "tmd", plot.type = "base")

plot(fit, which.plot = "pp", plot.type = "ggplot")
plot(fit, which.plot = "qq", plot.type = "ggplot")
# plot(fit, which.plot = "erp", plot.type = "ggplot")
plot(fit, which.plot = "exp", plot.type = "ggplot")
plot(fit, which.plot = "tmd", plot.type = "ggplot")


coef(fit)
expect_true(isTRUE(all.equal(deviance(fit), as.numeric(logLik(fit))*-2)))
expect_true(length(x) == nobs(fit))


# print.elife_npar
#
# print.elife_par_test
# print.elife_northropcoleman
# print.elife_profile
# confint.elife_profile
#
# plot.elife_hazard
# plot.elife_northropcoleman
# plot.elife_tstab
# plot.elife_profile
# plot.elife_ecdf
