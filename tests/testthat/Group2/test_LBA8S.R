cat("\n-------------------- Testing LBA Multi-core 8 Subjects ----------------")

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

  # pdf(file = "HLBA.pdf")
  # plot(pop.prior, ps = pop.mean)

  ## Simulate some data ----------
  dat <- simulate(model, nsub = 8, nsim = 10, prior = pop.prior)
  dmi <- BuildDMI(dat, model)
  ps <- attr(dat, "parameters")

  p.prior <- BuildPrior(
    dists = rep("tnorm", 9),
    p1   = pop.mean,
    p2   = pop.scale*5,
    lower=c(0,0,0,   .1, NA,NA,NA,NA, 0),
    upper=c(NA,NA,NA,NA, NA,NA,NA,NA, NA))

  ## Sampling separately ----------
  fit0 <- StartNewsamples(dmi, p.prior, ncore=4)
  fit  <- run(fit0, 5e2, ncore=2)
  fit  <- run(fit, 1e2, add=TRUE, ncore=4)

  est0 <- summary(fit, recovery = TRUE, ps = ps, verbose =TRUE)
