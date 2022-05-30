# Checking the nonparametric maximum
# likelihood estimator of Turnbull
# against the implementation.

# [1] Toy examples
#
# Two observations: A1: [1,3], A2: 4
# Probability of 0.5
time <- c(1,3)
time2 <- c(2,3)
test_simple1 <- np_elife(time,
         time2,
         event = c(3,1),
         type = "interval",
         vcov = FALSE)
test_simple2 <- npsurv(time,
       time2,
       event = c(3,1),
       type = "interval")
tinytest::expect_equal(test_simple1$prob, test_simple2$prob)

# Second example with right-censoring and a tie
rcens <- c(0, 1, 1)
time <- c(1, 1, 2)
# risk set at t=1 is 3, 1 failure
# risk set at t=2 is 1, 1 failure
test_simple3 <- npsurv(time,
                         event = rcens,
                         type = "right")
tinytest::expect_equal(c(1/3, 2/3),
                       test_simple3$prob)


# [2] Left truncated data
n <- 100L
u <- rnorm(n = n, sd = 1)
y <- u + rexp(n = n, rate = 1/10)
# Bell-Lynden estimator
utimes <- sort(unique(y))
ratio <- sapply(utimes, function(qj){
  1 - sum(y == qj)/
    sum(u <= qj & qj <=y)
})
density_bl <- -diff(c(1, cumprod(ratio)))
## This is unexpectedly long
npmle_ltrunc_long <-
  npsurv(time = y,
         ltrunc = u,
         tol = 1e-9)
tinytest::expect_equal(npmle_ltrunc_long$prob,
                       density_bl)

# [3] Left-truncated right censored data
# Create fake data with ties
set.seed(1234)
n <- 100L
ltrunc <- pmax(0, rnorm(n))
# Since about half of lower truncation bounds are null
# The following creates ties!
time <- rexp(10) + ltrunc
# 0 for right-censored, 1 for observed
rightcens <-
  sample(c(0L, 1L),
         size = n,
         replace = TRUE,
         prob = c(0.2, 0.8))
# But some tied times are censored!
time <- time + ifelse(rightcens == 0L, 1e-10, 0)
# The ECDF is supported on unique failure times
tun <- sort(unique(time[rightcens == 1L]))
# Product limit estimator
wgt <- sapply(tun, function(ti) {
  1 - sum(time == ti & rightcens == 1L) /
    sum(ltrunc <= ti & ti <= time)
  # & ifelse(rightcens == 1L, TRUE, ti < time))
})
# Test the function from the 'longevity' package
npmle_ltruncrcens_long <-
  npsurv(
    time = time,
    event = rightcens,
    type = "right",
    ltrunc = ltrunc
  )
npmle_ltruncrcens_long2 <-
  np_elife(
    time = time,
    event = rightcens,
    type = "right",
    ltrunc = ltrunc
  )
probs <- -diff(c(1, cumprod(wgt)))
# Expect only unique failure times
tinytest::expect_true(
  with(npmle_ltruncrcens_long,
     isTRUE(all(xval[,1] == xval[,2]))))
# Same right and left censoring
tinytest::expect_equal(npmle_ltruncrcens_long$xval[,1],
                       tun)
# Same probabilities
tinytest::expect_equal(npmle_ltruncrcens_long$prob,
                       probs)


# [4] right censored data (Kaplan-Meier)
n <- rpois(n = 1, lambda = 100)
time <- rpois(n, lambda = 100)
status <- ifelse(runif(n) < 0.5, 0, 1)
utime <- time[status == 1]
cstat <- .check_surv(time,
                     event = status,
                     type = 'right')
unex <- turnbull_intervals(
  time = cstat$time,
  time2 = cstat$time2,
  status = cstat$status)
cLimits <- .censTruncLimits(
  tsets = unex,
  n = length(cstat$time),
  rcens = cstat$time2,
  lcens = cstat$time,
  ltrunc = c(rep(-Inf, length(cstat$time))),
  rtrunc = c(rep(Inf, length(cstat$time))),
  trunc = FALSE)
# Fit model
test_rcens_long <-
  npsurv(time = time,
         event = status,
         type = 'right')

# Compare to Kaplan-Meier
km <- survival::survfit(survival::Surv(time, status) ~ 1)
probsKM <- -diff(c(1, km$surv))


if(max(time) == max(utime)){
  tinytest::expect_equal(test_rcens_long$xval[,1],
                         test_rcens_long$xval[,2])
  tinytest::expect_equal(km$time[probsKM > 0],
                         test_rcens_long$xval[,1])
  tinytest::expect_equal(test_rcens_long$prob,
                         probsKM[probsKM>0])
} else{
  # largest observation is right-censored
  nm <- nrow(test_rcens_long$xval)
  tinytest::expect_equal(test_rcens_long$xval[-nm, 1],
                         test_rcens_long$xval[-nm, 2])
  tinytest::expect_equal(km$time[probsKM > 0],
                         test_rcens_long$xval[-nm,1])
  tinytest::expect_equal(test_rcens_long$prob[-nm],
                         probsKM[probsKM>0])
}


# [5] Example of Frydman (1994)
time <- c(2,4,6)
time2 <- c(3,5,8)
ltrunc <- c(1,1,2.5)

F_int <- turnbull_intervals(
  time = time,
  time2 = time2,
  ltrunc = ltrunc,
  status = rep(3,3))
test <- npsurv(
  time = time,
  time2 = time2,
  event = rep(3,3),
  type = "interval",
  ltrunc = ltrunc)
Fref <- cbind(time, c(ltrunc[3], time2[-1]))
tinytest::expect_equal(F_int,
                       Fref,
                       check.attributes = FALSE)
tinytest::expect_equal(test$prob,
                       c(0.5,0.25,0.25))

