cat("\n-------------------- Testing HDDM --------------------")
  rm(list = ls())

  model <- BuildModel(
    p.map     = list(a = "1", v = "F", z = "1", d = "1", sz = "1", sv = "1",
                     t0 = "1", st0 = "1"),
    match.map = list(M = list(s1 = "r1", s2 = "r2")),
    factors   = list(S = c("s1", "s2"), F = c("f1", "f2")),
    constants = c(st0 = 0, d = 0),
    responses = c("r1", "r2"),
    type      = "rd")
  npar <- length(GetPNames(model))

  ## Population distribution
  pop.mean  <- c(a=2,   v.f1=4,  v.f2=3,  z=0.5, sz=0.3, sv=1,  t0=0.3)
  pop.scale <- c(a=0.5, v.f1=.5, v.f2=.5, z=0.1, sz=0.1, sv=.3, t0=0.05)
  pop.prior <- BuildPrior(
    dists = rep("tnorm", npar),
    p1    = pop.mean,
    p2    = pop.scale,
    lower = c(0,-5, -5, 0, 0, 0, 0),
    upper = c(5, 7,  7, 1, 2, 1, 1))

  ## Simulate some data
  dat <- simulate(model, nsub = 4, nsim = 10, prior = pop.prior)
  dmi <- BuildDMI(dat, model)
  ps <- attr(dat, "parameters")

  p.prior <- BuildPrior(
    dists = rep("tnorm", npar),
    p1    = pop.mean,
    p2    = pop.scale*5,
    lower = c(0,-5, -5, 0, 0, 0, 0),
    upper = c(5, 7,  7, 1, 2, 1, 1))

  mu.prior <- ggdmc::BuildPrior(
    dists = rep("tnorm", npar),
    p1    = pop.mean,
    p2    = pop.scale*5,
    lower = c(0,-5, -5, 0, 0, 0, 0),
    upper = c(5, 7,  7, 1, 2, 1, 1)
  )
  sigma.prior <- BuildPrior(
    dists = rep("beta", npar),
    p1    = c(a=1, v.f1=1,v.f2 = 1, z=1, sz=1, sv=1, t0=1),
    p2    = rep(1, npar),
    upper = rep(2, npar))

  priors <- list(pprior=p.prior, location=mu.prior, scale=sigma.prior)

  ## Sampling ------------
  fit0 <- StartNewsamples(dmi, priors)
  fit  <- run(fit0)
  fit  <- run(fit, 1e2, add=TRUE)

  res <- hgelman(fit, verbose = TRUE)
  est0 <- summary(fit, recovery = TRUE, ps = ps, verbose = TRUE)
  est1 <- summary(fit, hyper = TRUE, recovery = TRUE, ps = pop.mean,  type = 1, verbose = TRUE)
  est2 <- summary(fit, hyper = TRUE, recovery = TRUE, ps = pop.scale, type = 2, verbose = TRUE)

  pdf(file = "HDDM.pdf")
  p1 <- plot(fit0, hyper = TRUE, start = 51)
  p2 <- plot(fit, pll = FALSE)
  p3 <- plot(fit)
  dev.off()

  ##                    a    sv    sz    t0  v.f1  v.f2    z
  ## True            2.00  1.00  0.30  0.30  4.00  3.00 0.50
  ## 2.5% Estimate   1.35  0.13  0.01  0.20  3.06  2.02 0.12
  ## 50% Estimate    1.81  0.86  0.22  0.28  3.68  2.58 0.51
  ## 97.5% Estimate  2.13  0.99  0.43  0.33  4.26  3.23 0.90
  ## Median-True    -0.19 -0.14 -0.08 -0.02 -0.32 -0.42 0.01
  ##
  ##                    a    sv   sz   t0 v.f1 v.f2    z
  ## True            0.50  0.30 0.10 0.05 0.50 0.50 0.10
  ## 2.5% Estimate   0.21  0.05 0.09 0.04 0.39 0.42 0.07
  ## 50% Estimate    0.38  0.23 0.24 0.06 0.71 0.75 0.13
  ## 97.5% Estimate  1.00  1.88 0.59 0.89 1.46 1.53 1.83
  ## Median-True    -0.12 -0.07 0.14 0.01 0.21 0.25 0.03

