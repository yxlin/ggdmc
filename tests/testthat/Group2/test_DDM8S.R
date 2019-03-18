cat("\n-------------------- Testing DDM Multi-core 8 Subjects ----------------")
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
  dat <- simulate(model, nsub = 8, nsim = 10, prior = pop.prior)
  dmi <- BuildDMI(dat, model)
  ps <- attr(dat, "parameters")

  p.prior <- BuildPrior(
    dists = rep("tnorm", npar),
    p1    = pop.mean,
    p2    = pop.scale*5,
    lower = c(0,-5, -5, 0, 0, 0, 0),
    upper = c(5, 7,  7, 1, 2, 2, 1))

  ## Sampling separately ----------
  fit0 <- StartNewsamples(dmi, p.prior, ncore=2)
  fit  <- run(fit0, 5e2, ncore=2)
  fit  <- run(fit, 1e2, add=TRUE, ncore=2)

  est0 <- summary(fit, recovery = TRUE, ps = ps, verbose =TRUE)

