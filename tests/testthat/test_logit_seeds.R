require(ggdmc); require(data.table); require(ggplot2); require(gridExtra)
require(testthat)

  rm(list = ls())
  setwd("/media/yslin/KIWI/Documents/ggdmc/")
  dat <- list(S = c(10, 23, 23, 26, 17,  5, 53, 55, 32, 46,  ## Success
                    10,  8, 10,  8, 23,  0,  3, 22, 15, 32,
                    3),
              N = c(39, 62, 81, 51, 39,  6, 74, 72, 51, 79,
                    13, 16, 30, 28, 45,  4, 12, 41, 30, 51,
                    7),
              X1 = c(0, 0, 0, 0, 0,  0, 0, 0, 0, 0,  ## seed variety; 0 = aegytpiao 75 1 = aegyptiao 73
                     0, 1, 1, 1, 1,  1, 1, 1, 1, 1,
                     1),
              X2 = c(0, 0, 0, 0, 0,  1, 1, 1, 1, 1,  ## root extract; 0 = bean; 1 = cucumber
                     1, 0, 0, 0, 0,  0, 1, 1, 1, 1,
                     1),
              Ns = 21)

  d <- data.table(S = dat$X1, E = dat$X2, R = "r1", Y = dat$S, N = dat$N)
  dplyr::tbl_df(d)
  dat <- data.frame(d)

  model <- BuildModel(
    p.map     = list(a0 = "1", a1 = "1", a2 = "1", a3 = "1", b = "1"),
    match.map = NULL,
    regressors= NULL,
    factors   = list(S = c("s1", "s2"), E = c("e1", "e2")),
    responses = "r1",
    constants = NULL,
    type      = "logit")

dmi <- BuildDMI(data.frame(d), model)
dplyr::tbl_df(dmi)

npar <- length(GetPNames(model))
start  <- BuildPrior(
    dists = c(rep("tnorm2", npar)),
    p1    = c(a0 = -.556, a1 = .096, a2 = 1.344, a3 = -.824, b = .278),
    p2    = rep(1e-6, npar),
    lower = rep(NA, npar),
    upper = rep(NA, npar))

  p.prior  <-BuildPrior(
    dists = c(rep("tnorm2", npar)),
    p1    = c(a0 = 0, a1 = 0, a2 = 0, a3 = 0, b = 0),
    p2    = rep(1e-6, npar),
    lower = rep(NA, npar),
    upper = rep(NA, npar))

p.vector <- c(a0 = -.556, a1 = .096, a2 = 1.344, a3 = -.824, b = .278)
  plot(p.prior, ps = p.vector)
print(p.prior)

## Sampling -----------
fit0 <- run(Startlogit(5e2, dmi, start, p.prior, thin = 2))

fit <- fit0
fit <- run(RestartSamples(5e2, fit, thin = 8))
rhat <- gelman(fit, verbose = TRUE)


  cat("Done ", path[1], "\n")
  setwd("/media/yslin/KIWI/Documents/ggdmc/")

  ## Analysis -----------
  p0 <- plot(fit)
  p1 <- plot(fit, start = 101)
  ## plot(fit, pll =F, den=T)

est <- summary(fit, recover = TRUE, ps = p.vector, verbose = TRUE, digits = 3)
##                    a0     a1    a2     a3     b
## True           -0.556  0.096 1.344 -0.824 0.278
## 2.5% Estimate  -0.939 -0.408 0.954 -2.032 0.388
## 50% Estimate   -0.538  0.149 1.499 -1.253 0.500
## 97.5% Estimate -0.162  0.737 2.042 -0.482 0.692
## Median-True     0.018  0.053 0.155 -0.429 0.222
es <- effectiveSize(fit)
##       a0       a1       a2       a3        b 
## 3823.417 3797.726 3821.115 3752.681 3219.640

## Inference for Bugs model at "/tmp/Rtmp6gOxmy/model.txt", 
#             mean    sd   2.5%    25%     50%     75%   97.5%  Rhat n.eff
# alpha0    -0.557 0.197 -0.964 -0.675  -0.556  -0.430  -0.174 1.005  4500
# alpha1     0.086 0.317 -0.560 -0.114   0.096   0.286   0.688 1.002  4500
# alpha2     1.348 0.276  0.834  1.172   1.344   1.514   1.916 1.001  4500
# alpha12   -0.824 0.445 -1.736 -1.104  -0.820  -0.538   0.048 1.002  1900
# sigma      0.286 0.146  0.044  0.184   0.278   0.376   0.600 1.024   330
# deviance 101.989 6.824 90.355 96.937 101.400 106.700 115.600 1.001  4500

#            mean    sd    2.5%  97.5% n_eff
# alpha0  -0.5643 0.207 -0.9630 -0.156   470
# alpha1   0.0917 0.335 -0.5626  0.727   435
# alpha2   1.3642 0.285  0.7954  1.957   500
# alpha12 -0.8303 0.450 -1.6577  0.088   552
# sigma    0.2982 0.137  0.0869  0.621   163
