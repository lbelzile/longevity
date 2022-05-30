# # Optimization routines
# # Check that we get the same values for the MLE
# # as for "mev" package for GP and for special
# # cases of the exponential distribution where
# # the solution is available in closed-form
#
library(mev)
library(evd)
library(fitdistrplus)
set.seed(1234)
n <- 1e6L
lower <- runif(n, max = 10)
upper <- runif(n, min = 10, max = 20)
samp1 <- longevity::samp_elife(
  n = 1e6,
  scale = 10,
  shape = 0.1,
  lower = lower,
  upper = upper,
  family = "gp",
  type2 = "ltrt")
fit1 <- longevity::fit_elife(
  time = samp1,
  ltrunc = lower,
  rtrunc = upper,
  family = "gp")

samp2 <- evd::rgpd(1e3, loc = 2, scale = 2, shape = -0.1)
fit_2a <- longevity::fit_elife(
  time = samp2,
  ltrunc = rep(0, length(samp2)),
  rtrunc = rep(Inf, length(samp2)),
  thresh = 2,
  family = "gp")
fit_2b <- longevity::fit_elife(
  time = samp2,
  thresh = 2,
  family = "gp")
fit_2c <- mev::fit.gpd(
  xdat = samp2,
  threshold = 2)

fit_2a$par - fit_2c$par
fit_2a$par - fit_2b$par

samp3 <- rexp(n = 100, rate = 0.5)
fit_3a <- longevity::fit_elife(
  time = samp3,
  family = "exp")
fit_3a$par - mean(samp3)

samp4 <- longevity::samp_elife(
  n = 100,
  scale = 2,
  upper = 3,
  family = "exp",
  type2 = "ltrc")
fit4 <- longevity::fit_elife(
  time = samp4$dat,
  event = !samp4$rcens,
  type = "right",
  family = "exp")
fit4$par - sum(samp4$dat)/sum(!samp4$rcens)

# Check maximum likelihood without censoring or truncation
# Check that specifying truncation bounds that are 0/Inf
# does not impact the MLE
# Check that censored data with no censoring gives
# the same point estimates

