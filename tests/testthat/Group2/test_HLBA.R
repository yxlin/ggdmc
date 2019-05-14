cat("\n-------------------- Testing HLBA --------------------")
  rm(list = ls())

  # 2x2 (Stimulus and factor F) specified in model.dmc call
  model <- BuildModel(
    p.map     = list(A = "1", B = "R", t0 = "1", mean_v = c("D", "M"),
                     sd_v = "M", st0 = "1"),
    match.map = list(M = list(s1 = 1, s2 = 2)),
    factors   = list(S = c("s1", "s2"), D = c("d1", "d2")),
    constants = c(sd_v.false = 1, st0 = 0),
    responses = c("r1", "r2"),
    type      = "norm")

  ## Population distribution, rate effect on F
  pop.mean <- c(A=.4, B.r1=.6, B.r2=.8, t0=.3,
                mean_v.d1.true  = 1.5,
                mean_v.d2.true  = 1.0,
                mean_v.d1.false = .15,
                mean_v.d2.false = .2,  sd_v.true = .25)
  pop.scale <-c(A=.1, B.r1=.1, B.r2=.1, t0=.05,
                mean_v.d1.true  =.2,
                mean_v.d2.true  =.2,
                mean_v.d1.false =.2,
                mean_v.d2.false =.2,  sd_v.true = .1)
  pop.prior <- BuildPrior(
    dists = rep("tnorm", 9),
    p1 = pop.mean,
    p2 = pop.scale,
    lower = c(0,0,0,   .1, NA,NA,NA,NA, 0),
    upper = c(NA,NA,NA, 1, NA,NA,NA,NA, NA))


  ## Simulate some data ----------
  dat <- simulate(model, nsub = 30, nsim = 30, prior = pop.prior)
  dmi <- BuildDMI(dat, model)
  ps <- attr(dat, "parameters")

  p.prior <- BuildPrior(
    dists = rep("tnorm", 9),
    p1   = pop.mean,
    p2   = pop.scale*5,
    lower=c(0,0,0,   .1, NA,NA,NA,NA, 0),
    upper=c(NA,NA,NA,NA, NA,NA,NA,NA, NA))
  mu.prior <- BuildPrior(
    dists = rep("tnorm",  9),
    p1    = pop.mean,
    p2    = c(1,   1,  1,  1,   2,  2,  2, 2,  1),
    lower = c(0,   0,  0, .01, NA, NA, NA, NA, 0),
    upper = c(NA, NA, NA,  NA, NA, NA, NA, NA, NA))
  sigma.prior <- BuildPrior(
    dists = rep("beta", length(p.prior)),
    p1    = c(A = 1, B.r1 = 1, B.r2 = 1, t0 = 1, mean_v.d1.true = 1,
              mean_v.d2.true = 1, mean_v.d1.false = 1, mean_v.d2.false = 1,
              sd_v.true = 1),
    p2    = rep(1, 9))

  ## Sampling -------------
  priors <- list(pprior=p.prior, location=mu.prior, scale=sigma.prior)

  fit0 <- StartNewsamples(dmi, priors, thin = 8)
  fit  <- run(fit0, thin = 8)
  fit  <- run(fit, 1e2, add=TRUE)

  pdf(file = "HLBA.pdf")
  plot(pop.prior, ps = pop.mean)
  plot(p.prior, ps = ps)
  plot(sigma.prior, ps = pop.scale)

  p0 <- plot(fit, hyper = TRUE)
  p1 <- plot(fit, hyper = TRUE, den = TRUE, pll=FALSE)
  dev.off()

  ## Analysis -----------
  res <- hgelman(fit, verbose = TRUE)
  est0 <- summary(fit, recovery = TRUE, ps = ps, verbose = TRUE)
  est1 <- summary(fit, hyper = TRUE, recovery = TRUE, ps = pop.mean,  type = 1, verbose = TRUE)
  est2 <- summary(fit, hyper = TRUE, recovery = TRUE, ps = pop.scale, type = 2, verbose = TRUE)

  for(i in 1:length(fit))
  {
      est <- summary(fit[[i]], recovery = TRUE, ps = ps[i,], verbose=TRUE)
  }

  #                   A  B.r1  B.r2 mean_v.d1.false mean_v.d1.true
  # True           0.40  0.60  0.80            0.00           1.50
  # 2.5% Estimate  0.17  0.18  0.69           -0.36           1.35
  # 50% Estimate   0.44  0.58  0.76           -0.10           1.53
  # 97.5% Estimate 0.57  0.75  0.84            0.14           1.72
  # Median-True    0.04 -0.02 -0.04           -0.10           0.03
  #                mean_v.d2.false mean_v.d2.true sd_v.true    t0
  # True                      0.00           1.00      0.25  0.30
  # 2.5% Estimate            -0.46           0.85      0.03  0.23
  # 50% Estimate             -0.08           1.02      0.24  0.29
  # 97.5% Estimate            0.25           1.23      0.35  2.38
  # Median-True              -0.08           0.02     -0.01 -0.01

  #                   A B.r1  B.r2 mean_v.d1.false mean_v.d1.true
  # True           0.10 0.10  0.10            0.20           0.20
  # 2.5% Estimate  0.06 0.11  0.02            0.05           0.12
  # 50% Estimate   0.13 0.19  0.05            0.15           0.21
  # 97.5% Estimate 0.49 0.64  0.13            0.49           0.44
  # Median-True    0.03 0.09 -0.05           -0.05           0.01
  #                mean_v.d2.false mean_v.d2.true sd_v.true    t0
  # True                      0.20           0.20      0.10  0.05
  # 2.5% Estimate             0.19           0.13      0.07  0.02
  # 50% Estimate              0.40           0.23      0.14  0.04
  # 97.5% Estimate            0.84           0.48      0.40  0.14
  # Median-True               0.20           0.03      0.04 -0.01