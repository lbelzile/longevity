# Test functions in distributions.R

# Check parameter constraints

# Create fake ordered data
ddata <- c(-1, 0, sort(rexp(100)))
endpt_gp_neg <- c(4, 4 + 1e-8)
endpt_extgp_neg <- c(4 * log(2), 4 * log(2) + 1e-8)
# Evaluate distribution functions
p_exp_v <- pexp(ddata, rate = 1 / 2)
p_gpd_v_pos <- longevity::pgpd(ddata,
                               scale = 2,
                               shape = 0.5)
p_gpd_v_neg <- longevity::pgpd(ddata,
                               scale = 2,
                               shape = -0.5)
p_gpd_v_exp <- longevity::pgpd(ddata,
                               scale = 2,
                               shape = 0)
p_gomp_v <- longevity::pgomp(ddata,
                             scale = 2,
                             shape = 0.5)
p_gomp_v_exp <- longevity::pgomp(ddata,
                                 scale = 2,
                                 shape = 0)
p_gompmake_v <- longevity::pgompmake(ddata,
                                     scale = 2,
                                     shape = 0.5,
                                     lambda = 0.5)
p_gompmake_v_gomp <- longevity::pgompmake(ddata,
                                          scale = 2,
                                          shape = 0.5,
                                          lambda = 0)
p_gompmake_v_com <- longevity::pgompmake(ddata,
                                         scale = 2,
                                         shape = 0,
                                         lambda = 0.5)
p_gompmake_v_exp <- longevity::pgompmake(ddata,
                                         scale = 2,
                                         shape = 0,
                                         lambda = 0)
p_extgp_v_neg <- longevity::pextgp(ddata,
                                   scale = 2,
                                   shape1 = 0.5,
                                   shape2 = -0.5)
p_extgp_v_pos <- longevity::pextgp(ddata,
                                   scale = 2,
                                   shape1 = 0.5,
                                   shape2 = 0.5)
p_extgp_v_gpd <- longevity::pextgp(ddata,
                                   scale = 2,
                                   shape1 = 0,
                                   shape2 = -0.5)
p_extgp_v_gomp <- longevity::pextgp(ddata,
                                    scale = 2,
                                    shape1 = 0.5,
                                    shape2 = 0)



p_exp_r <- pexp(ddata,
                rate = 1 / 2,
                lower.tail = FALSE,
                log.p = TRUE)
p_gpd_r_pos <- longevity::pgpd(
  ddata,
  scale = 2,
  shape = 0.5,
  lower.tail = FALSE,
  log.p = TRUE
)
p_gpd_r_neg <- longevity::pgpd(
  ddata,
  scale = 2,
  shape = -0.5,
  lower.tail = FALSE,
  log.p = TRUE
)
p_gpd_r_exp <- longevity::pgpd(
  ddata,
  scale = 2,
  shape = 0,
  lower.tail = FALSE,
  log.p = TRUE
)
p_gomp_r <- longevity::pgomp(
  ddata,
  scale = 2,
  shape = 0.5,
  lower.tail = FALSE,
  log.p = TRUE
)
p_gomp_r_exp <- longevity::pgomp(
  ddata,
  scale = 2,
  shape = 0,
  lower.tail = FALSE,
  log.p = TRUE
)
p_gompmake_r <- longevity::pgompmake(
  ddata,
  scale = 2,
  shape = 0.5,
  lambda = 0.5,
  lower.tail = FALSE,
  log.p = TRUE
)
p_gompmake_r_gomp <- longevity::pgompmake(
  ddata,
  scale = 2,
  shape = 0.5,
  lambda = 0,
  lower.tail = FALSE,
  log.p = TRUE
)
p_gompmake_r_com <- longevity::pgompmake(
  ddata,
  scale = 2,
  shape = 0,
  lambda = 0.5,
  lower.tail = FALSE,
  log.p = TRUE
)
p_gompmake_r_exp <- longevity::pgompmake(
  ddata,
  scale = 2,
  shape = 0,
  lambda = 0,
  lower.tail = FALSE,
  log.p = TRUE
)
p_extgp_r_neg <- longevity::pextgp(
  ddata,
  scale = 2,
  shape1 = 0.5,
  shape2 = -0.5,
  lower.tail = FALSE,
  log.p = TRUE
)
p_extgp_r_pos <- longevity::pextgp(
  ddata,
  scale = 2,
  shape1 = 0.5,
  shape2 = 0.5,
  lower.tail = FALSE,
  log.p = TRUE
)
p_extgp_r_gpd <- longevity::pextgp(
  ddata,
  scale = 2,
  shape1 = 0,
  shape2 = -0.5,
  lower.tail = FALSE,
  log.p = TRUE
)
p_extgp_r_gomp <- longevity::pextgp(
  ddata,
  scale = 2,
  shape1 = 0.5,
  shape2 = 0,
  lower.tail = FALSE,
  log.p = TRUE
)
# Check that distribution functions evaluate to 0 below origin
models <-
  c(
    "p_gpd_v_neg",
    "p_gpd_v_pos",
    "p_gpd_v_exp",
    "p_gomp_v",
    "p_gomp_v_exp",
    "p_gompmake_v",
    "p_gompmake_v_exp",
    "p_gompmake_v_gomp",
    "p_gompmake_v_com",
    "p_extgp_v_gomp",
    "p_extgp_v_gpd",
    "p_extgp_v_neg",
    "p_extgp_v_pos"
  )

zerocdf <- sapply(models,
                  function(model) {
                    all.equal(get(model)[1:2], c(0, 0))
                  })
expect_true(current = isTRUE(all(zerocdf)),
            info = "Distribution function is zero for negative values.")

# Distribution function is nondecreasing

nondecreasing <- sapply(models,
                        function(model) {
                          isTRUE(all(diff(get(model)) >= 0))
                        })
expect_true(current = isTRUE(all(nondecreasing)),
            info = "Distribution function is nondecreasing.")

# Distribution function maps to unit interval
unitinterval <- sapply(models,
                       function(model) {
                         isTRUE(all(get(model) >= 0 & get(model) <= 1))
                       })
expect_true(current = isTRUE(all(unitinterval)),
            info = "Distribution function maps to [0,1].")


# Check special cases
expect_equal(p_exp_v,
             p_gomp_v_exp,
             info = "CDF: gomp vs exp")
expect_equal(p_exp_v,
             p_gompmake_v_exp,
             info = "CDF: gomp vs exp")
expect_equal(p_exp_v,
             p_gpd_v_exp,
             info = "CDF: gp vs exp")
expect_equal(p_gompmake_v_gomp, p_gomp_v,
             info = "CDF: gompmake vs gomp")
expect_equal(p_extgp_v_gpd, p_gpd_v_neg,
             info = "CDF: extgp vs gp")
expect_equal(p_extgp_v_gomp, p_gomp_v,
             info = "CDF: extgp vs gomp")

# Check log survival matches

logsurv <-
  sapply(models,
         function(model) {
           all.equal(log(1 - get(model)),
                     get(sub(
                       x = model,
                       pattern = "_v",
                       replacement = "_r"
                     )))
         })
expect_true(current = isTRUE(all(logsurv)),
            info = "log survival functions match.")


# Check that density is zero outside of the domain

# Evaluate density functions
d_exp_v <- dexp(ddata, rate = 1 / 2)
d_gpd_v_pos <- longevity::dgpd(ddata,
                               scale = 2,
                               shape = 0.5)
d_gpd_v_neg <- longevity::dgpd(ddata,
                               scale = 2,
                               shape = -0.5)
d_gpd_v_exp <- longevity::dgpd(ddata,
                               scale = 2,
                               shape = 0)
d_gomp_v <- longevity::dgomp(ddata,
                             scale = 2,
                             shape = 0.5)
d_gomp_v_exp <- longevity::dgomp(ddata,
                                 scale = 2,
                                 shape = 0)
d_gompmake_v <- longevity::dgompmake(ddata,
                                     scale = 2,
                                     shape = 0.5,
                                     lambda = 0.5)
d_gompmake_v_gomp <- longevity::dgompmake(ddata,
                                          scale = 2,
                                          shape = 0.5,
                                          lambda = 0)
d_gompmake_v_com <- longevity::dgompmake(ddata,
                                         scale = 2,
                                         shape = 0,
                                         lambda = 0.5)
d_gompmake_v_exp <- longevity::dgompmake(ddata,
                                         scale = 2,
                                         shape = 0,
                                         lambda = 0)
d_extgp_v_neg <- longevity::dextgp(ddata,
                                   scale = 2,
                                   shape1 = 0.5,
                                   shape2 = -0.5)
d_extgp_v_pos <- longevity::dextgp(ddata,
                                   scale = 2,
                                   shape1 = 0.5,
                                   shape2 = 0.5)
d_extgp_v_gpd <- longevity::dextgp(ddata,
                                   scale = 2,
                                   shape1 = 0,
                                   shape2 = -0.5)
d_extgp_v_gomp <- longevity::dextgp(ddata,
                                    scale = 2,
                                    shape1 = 0.5,
                                    shape2 = 0)

models <-
  c(
    "d_gpd_v_neg",
    "d_gpd_v_pos",
    "d_gpd_v_exp",
    "d_gomp_v",
    "d_gomp_v_exp",
    "d_gompmake_v",
    "d_gompmake_v_exp",
    "d_gompmake_v_gomp",
    "d_gompmake_v_com",
    "d_extgp_v_gomp",
    "d_extgp_v_gpd",
    "d_extgp_v_neg",
    "d_extgp_v_pos"
  )

zerodens <- sapply(models,
                   function(model) {
                     all.equal(get(model)[1], 0)
                   })
expect_true(current = isTRUE(all(zerodens)),
            info = "Density is zero for negative values.")

posdens <- sapply(models,
                  function(model) {
                    isTRUE(all(get(model) >= 0))
                  })
expect_true(current = isTRUE(all(posdens)),
            info = "Density is non-negative.")

# Density
expect_equal(dgpd(
  endpt_gp_neg,
  loc = 0,
  scale = 2,
  shape = -0.5
), c(0, 0),
info = "GP density at endpoint is zero with negative shape")

expect_equal(dextgp(
  endpt_extgp_neg,
  scale = 2,
  shape1 = 0.5,
  shape2 = -0.5
),
c(0, 0),
info = "extgp density at endpoint is zero with negative shape")

# Check special cases
expect_equal(d_exp_v,
             d_gomp_v_exp,
             info = "density: gomp vs exp")
expect_equal(d_exp_v,
             d_gompmake_v_exp,
             info = "density: gomp vs exp")
expect_equal(d_exp_v,
             d_gpd_v_exp,
             info = "density: gp vs exp")
expect_equal(d_gompmake_v_gomp,
             d_gomp_v,
             info = "density: gompmake vs gomp")
expect_equal(d_extgp_v_gpd,
             d_gpd_v_neg,
             info = "density: extgp vs gp")
expect_equal(d_extgp_v_gomp,
             d_gomp_v,
             info = "density: extgp vs gomp")

# Check that CDF is inverse of quantile function
unif <- sort(runif(100))
scale_p <- rexp(1, 10)
shape_p_gp <- runif(1, min = -1, max = 1)
shape_p_gomp <- runif(1)
lambda_p <- rexp(1, rate = 1/10)

gp_inv <- pgpd(qgpd(unif, scale = scale_p, shape = shape_p_gp),
               scale = scale_p,
               shape = shape_p_gp)
expect_equal(gp_inv, unif,
             info = "gpd: inverse function")
gomp_inv <-
  pgomp(qgomp(unif, scale = scale_p, shape = shape_p_gomp),
        scale = scale_p,
        shape = shape_p_gomp)
expect_equal(gomp_inv, unif,
             info = "gomp: inverse function")

gompmake_inv <-
  pgompmake(
    qgompmake(
      unif,
      scale = scale_p,
      shape = shape_p_gomp,
      lambda = lambda_p
    ),
    scale = scale_p,
    shape = shape_p_gomp,
    lambda_p
  )
expect_equal(gompmake_inv, unif,
             info = "gompmake: inverse function")

extgp_inv <- pextgp(
  qextgp(
    unif,
    scale = scale_p,
    shape1 = shape_p_gomp,
    shape2 = shape_p_gp
  ),
  scale = scale_p,
  shape1 = shape_p_gomp,
  shape2 = shape_p_gp
)
expect_equal(extgp_inv, unif,
             info = "extgp: inverse function")

# Check log density

d_gpd_l_pos <- longevity::dgpd(ddata,
                               scale = 2,
                               shape = 0.5,
                               log = TRUE)
d_gpd_l_neg <- longevity::dgpd(ddata,
                               scale = 2,
                               shape = -0.5,
                               log = TRUE)
d_gpd_l_exp <- longevity::dgpd(ddata,
                               scale = 2,
                               shape = 0,
                               log = TRUE)
d_gomp_l <- longevity::dgomp(ddata,
                             scale = 2,
                             shape = 0.5,
                             log = TRUE)
d_gomp_l_exp <- longevity::dgomp(ddata,
                                 scale = 2,
                                 shape = 0,
                                 log = TRUE)
d_gompmake_l <- longevity::dgompmake(
  ddata,
  scale = 2,
  shape = 0.5,
  lambda = 0.5,
  log = TRUE
)
d_gompmake_l_gomp <- longevity::dgompmake(
  ddata,
  scale = 2,
  shape = 0.5,
  lambda = 0,
  log = TRUE
)
d_gompmake_l_com <- longevity::dgompmake(
  ddata,
  scale = 2,
  shape = 0,
  lambda = 0.5,
  log = TRUE
)
d_gompmake_l_exp <- longevity::dgompmake(
  ddata,
  scale = 2,
  shape = 0,
  lambda = 0,
  log = TRUE
)
d_extgp_l_neg <- longevity::dextgp(
  ddata,
  scale = 2,
  shape1 = 0.5,
  shape2 = -0.5,
  log = TRUE
)
d_extgp_l_pos <- longevity::dextgp(
  ddata,
  scale = 2,
  shape1 = 0.5,
  shape2 = 0.5,
  log = TRUE
)
d_extgp_l_gpd <- longevity::dextgp(
  ddata,
  scale = 2,
  shape1 = 0,
  shape2 = -0.5,
  log = TRUE
)
d_extgp_l_gomp <- longevity::dextgp(
  ddata,
  scale = 2,
  shape1 = 0.5,
  shape2 = 0,
  log = TRUE
)
logdens <-
  sapply(models,
         function(model) {
           all.equal(log(get(model)),
                     get(sub(
                       x = model,
                       pattern = "_v",
                       replacement = "_l"
                     )))
         })
expect_true(current = isTRUE(all(logdens)),
            info = "log density matches.")

# Density evaluated at NA/-Inf/Inf
nfnt <- c(NA, -Inf, Inf)
nfnt_ev <- c(NA, 0, 0)
expect_equal(nfnt_ev,
             dgpd(x = nfnt,
                  shape = runif(1, -1, 1)),
             info = "gpd dens at non-finite values")
expect_equal(nfnt_ev,
             dgompmake(
               x = nfnt,
               shape = runif(1),
               lambda = runif(1)
             ),
             info = "gompmake dens at non-finite values")
expect_equal(nfnt_ev,
             dgomp(x = nfnt,
                   shape = runif(1)),
             info = "gomp dens at non-finite values")
expect_equal(nfnt_ev,
             dextgp(
               x = nfnt,
               shape1 = runif(1),
               shape2 = runif(1, -1, 1)
             ),
             info = "extgp dens at non-finite values")


# Check distribution function at endpoints + NA
expect_equal(c(NA, 0, 1),
             pgpd(nfnt,
                  shape = runif(1, -1, 1)),
             info = "gpd CDF at non-finite values")
expect_equal(c(NA, 0, 1),
             pgompmake(nfnt,
                       shape = runif(1),
                       lambda = runif(1)),
             info = "gompmake CDF at non-finite values")
expect_equal(c(NA, 0, 1),
             pgomp(nfnt,
                   shape = runif(1)),
             info = "gomp CDF at non-finite values")
expect_equal(c(NA, 0, 1),
             pextgp(nfnt,
                    shape1 = runif(1),
                    shape2 = runif(1, -1, 1)),
             info = "extgp CDF at non-finite values")


# Check invalid parameters yield errors
expect_error(dextgp(runif(10), shape1 = -0.1), info = "Invalid parameter")
expect_error(dextgp(runif(10), scale = -1), info = "Invalid parameter")
expect_error(dgomp(runif(10), shape = -0.1), info = "Invalid parameter - gompertz")
expect_error(dgompmake(runif(10), lambda = -0.1), info = "Invalid parameter - gompmake")
expect_error(dgompmake(runif(10), shape = -0.1), info = "Invalid parameter - gompmake")


# Check that *elife gives the same thing
# as functions called
expect_equal(
  dgompmake(ddata,
            lambda =  0.1,
            shape = 0.2),
  delife(
    ddata,
    scale = c(1, 0.1),
    shape = 0.2,
    family = "gompmake"
  ),
  info = "delife (gompmake)"
)

expect_equal(
  dextgp(ddata,
         shape1 = 0.1,
         shape2 = 0.5),
  delife(ddata,
         shape = c(0.1, 0.5),
         family = "extgp"),
  info = "delife (extgp)"
)

expect_equal(
  dgomp(ddata,
        shape = 0.1,
        scale = 2),
  delife(
    ddata,
    shape = 0.1,
    scale = 2,
    family = "gomp"
  ),
  info = "delife (gomp)"
)


expect_equal(
  pgompmake(ddata,
            lambda =  0.1,
            shape = 0.2),
  pelife(
    ddata,
    scale = c(1, 0.1),
    shape = 0.2,
    family = "gompmake"
  ),
  info = "pelife (gompmake)"
)

expect_equal(
  pextgp(ddata,
         shape1 = 0.1,
         shape2 = 0.5),
  pelife(ddata,
         shape = c(0.1, 0.5),
         family = "extgp"),
  info = "pelife (extgp)"
)

expect_equal(
  pgomp(ddata,
        shape = 0.1,
        scale = 2),
  pelife(
    ddata,
    shape = 0.1,
    scale = 2,
    family = "gomp"
  ),
  info = "pelife (gomp)"
)
