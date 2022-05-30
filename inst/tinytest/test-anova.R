# Check nesting structure

library(longevity)
library(tinytest)
set.seed(123)
n <- rpois(n = 1, lambda = 150)
lower <- ifelse(runif(n) < 0.8,
                0,
                qgomp(runif(n, max = 0.9),
                      scale = 1,
                      shape = 1.5))
upper <- lower + rexp(n, rate = 1/10)
samp <- longevity::samp_elife(
  n = n,
  scale = 2,
  shape = 1.5,
  type2 = "ltrc",
  lower = lower,
  upper = upper,
  family = "gomp")
tinytest::expect_true(
  isTRUE(all(samp$dat - lower > 0)))
# Fit the model with all families
fit_exp <-
  fit_elife(time = samp$dat,
            type = "right",
            ltrunc = lower,
            event = samp$rcens == 1L,
            family = "exp")
fit_gomp <-
  fit_elife(time = samp$dat,
            type = "right",
            ltrunc = lower,
            event = samp$rcens == 1L,
            family = "gomp")
fit_gp <-
  fit_elife(time = samp$dat,
            type = "right",
            ltrunc = lower,
            event = samp$rcens == 1L,
            family = "gp")
fit_gompmake <-
  fit_elife(time = samp$dat,
            type = "right",
            ltrunc = lower,
            event = samp$rcens == 1L,
            family = "gompmake")
fit_extgp <-
  fit_elife(time = samp$dat,
            type = "right",
            ltrunc = lower,
            event = samp$rcens == 1L,
            family = "extgp")
fit_weibull <-
  fit_elife(time = samp$dat,
            type = "right",
            ltrunc = lower,
            event = samp$rcens == 1L,
            family = "weibull")
devs <-
c(deviance(fit_exp),
  deviance(fit_gp),
  deviance(fit_gomp),
  deviance(fit_gompmake),
  deviance(fit_weibull),
  deviance(fit_extgp))
length(devs) == length(unique(devs))

# Regular testing
tinytest::expect_equal(
  anova(fit_exp, fit_gp)[2,5],
  pchisq(deviance(fit_exp) - deviance(fit_gp),
         df = 1, lower.tail = FALSE))
tinytest::expect_equal(
  anova(fit_exp, fit_weibull)[2,5],
  pchisq(deviance(fit_exp) -
           deviance(fit_weibull),
       df = 1, lower.tail = FALSE))
tinytest::expect_equal(
  anova(fit_gomp, fit_extgp)[2,5],
  pchisq(deviance(fit_gomp) - deviance(fit_extgp),
         df = 1, lower.tail = FALSE))


# Nonregular, one parameter
tinytest::expect_equal(
  anova(fit_exp, fit_gomp)[2,5],
  0.5*pchisq(deviance(fit_exp) - deviance(fit_gomp),
         df = 1, lower.tail = FALSE)
  )
tinytest::expect_equal(
  anova(fit_gp, fit_extgp)[2,5],
  0.5*pchisq(deviance(fit_gp) - deviance(fit_extgp),
             df = 1, lower.tail = FALSE)
)
# Non-regular comparison
tinytest::expect_equal(
  anova(fit_exp, fit_gompmake)[2,5],
  0.5*pchisq(deviance(fit_exp) - deviance(fit_gompmake),
         df = 1, lower.tail = FALSE) +
    0.25*pchisq(deviance(fit_exp) - deviance(fit_gompmake),
               df = 2, lower.tail = FALSE))


# Invalid comparisons - non-nested models
tinytest::expect_error(
  anova(fit_weibull, fit_gp))
tinytest::expect_error(
  anova(fit_weibull, fit_extgp))
tinytest::expect_error(
  anova(fit_weibull, fit_gompmake))
tinytest::expect_error(
  anova(fit_weibull, fit_gomp))
tinytest::expect_error(
  anova(fit_gomp, fit_gp))
tinytest::expect_error(
  anova(fit_gompmake, fit_gp))
tinytest::expect_error(
  anova(fit_gompmake, fit_extgp))
