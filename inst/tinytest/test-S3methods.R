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
fit2 <- fit_elife(
  time = x,
  ltrunc = low,
  rtrunc = upp,
  family = "extweibull",
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

np_out <- longevity::np_elife(time = x,
                    ltrunc = low,
                    rtrunc = upp)
print(np_out)
summary(np_out)
plot(np_out)
plot(np_out$cdf, main = "")

anova(fit, fit2)
covar <- sample(c(0,1), size = length(x), replace = TRUE)
fit3 <- test_elife(
  time = x,
  covariate = covar,
  ltrunc = low,
  rtrunc = upp,
  family = "exp")
print(fit3)

#
# # The following returns an error message
# test <- hazard_elife(x = 1.5,
#                      time = x,
#                      ltrunc = low,
#                      rtrunc = upp,
#                      family = "gp")
# hazard_elife(x = seq(0.2, 3, by = 0.25),
#              time = x,
#              ltrunc = low,
#              rtrunc = upp,
#              family = "gp")
set.seed(1234)
x <- samp_elife(n = n,
                scale = 2,
                shape = -0.2,
                lower = low <- runif(n),
                upper = upp <- runif(n, min = 3, max = 20),
                type2 = "ltrt",
                family = "gp")

test <- nc_test(
  time = x,
  ltrunc = low,
  rtrunc = upp,
  thresh = quantile(x, seq(0, 0.5, by = 0.1)))
print(test)
plot(test)
plot(test, plot.type = "ggplot")

tstab_plot <- tstab(time = x,
      ltrunc = low,
      rtrunc = upp,
      thresh = quantile(x, seq(0, 0.5, by = 0.1)))
plot(tstab_plot, plot.type = "ggplot")

endpt <- prof_gp_endpt(time = x,
              ltrunc = low,
              rtrunc = upp,
              psi = seq(max(x) + 1e-4, 70, length.out = 51))
print(endpt)
plot(endpt)
plot(endpt, plot.type = "ggplot")
confint(endpt)

