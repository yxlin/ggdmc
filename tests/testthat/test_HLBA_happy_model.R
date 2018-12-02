rm(list = ls())
require(ggdmc)
setwd("/media/yslin/KIWI/Documents/ggdmc")
  
## (P)eople factor, p1: young+female, p2: young+male, p3: old-female, p4: old-male
## (S)timuli factor, s0: 0% happiness, s1: 30%, s2: 60%, s3: 90%, s4: 100%
## (R)esponse, r1: sad face, r2: happy face
## S has no influence3
  model <- BuildModel(
    p.map     = list(A = "1", B = "1", t0 = "1", mean_v = c("P", "M"), 
                     sd_v = "1", st0 = "1"),
    match.map = list(M = list(s0 = 1, s1 = 1, s2 = 2, s3 = 2, s4 = 2)),
    factors   = list(S = c("s0", "s1", "s2", "s3", "s4"),
                     P = c("p1", "p2", "p3", "p4")),
    constants = c(st0 = 0, sd_v = 1),
    responses = c("r1", "r2"),
    type      = "norm")

## Presume effect sizes at the population distribution, P factor
## affect drift rates. I presume (1) young faces (vs old faces) tend to
## elicit happy responses and (2) female faces (vs male faces) tend to
## elicit happy responses and (3) gender plays a more significant role
## than age; hence, p1 > p3 > p2 > p4
npar <- length(GetPNames(model))
pnames <- GetPNames(model)
pop.mean <- c(A = .4, B = .6, t0 = .1, 
              mean_v.p1.true  = 5.5, mean_v.p2.true  = 4.0,
              mean_v.p3.true  = 2.5, mean_v.p4.true  = 0.5,
              mean_v.p1.false = 1.1, mean_v.p2.false  = 1.2,
              mean_v.p3.false = 1.3, mean_v.p4.false  = 1.4)
pop.scale <- rtnorm(npar, 0, 1, 0, 10)
names(pop.scale) <- pnames

pop.prior <- BuildPrior(
    dists = rep("tnorm", npar),
    p1 = pop.mean,
    p2 = pop.scale,
    lower = c(c(0,0,   .1), rep(NA, 8)),
    upper = c(c(NA,NA, 1, rep(NA, 8))))
plot(pop.prior, ps = pop.mean)
  
## Simulate some data ----------
dat <- simulate(model, nsub = 30, nsim = 100, prior = pop.prior)
dmi <- BuildDMI(dat, model)
ps <- attr(dat, "parameters")
round(colMeans(ps), 2)
round(matrixStats::colSds(ps), 2)

mu.prior <- BuildPrior(
    dists = rep("tnorm", npar),
    p1 = pop.mean,
    p2 = pop.scale * 10,
    lower = c(c(0,0,   .1), rep(NA, 8)),
    upper = c(c(NA,NA, 1, rep(NA, 8))))
plot(mu.prior, ps = ps)

sigma.prior <- BuildPrior(
    dists = rep("beta", npar),
    p1    = rep(1, npar),
    p2    = rep(1, npar))
names(sigma.prior) <- pnames
plot(sigma.prior)
pp.prior <- list(mu.prior, sigma.prior)

## Sampling ------------
setwd("/media/yslin/KIWI/Documents/ggdmc_lesson/")
path <- c("data/test_HLBA_happy_model.RData")

  fit0 <- run(StartNewHypersamples(5e2, dmi, pop.prior, pp.prior))
  fit  <- fit0
  thin <- 1
  repeat {
    fit <- run(RestartHypersamples(5e2, fit), pm = .05, hpm = .05)
    save(fit0, fit, file = path[1])
    rhat <- hgelman(fit, verbose = TRUE)
    if (all(rhat < 1.1)) break
    thin <- thin * 2
  }
  cat("Done ", path[1], "\n")
  setwd("/media/yslin/KIWI/Documents/ggdmc/")
  
  p0 <- plot(fit, hyper = TRUE)
  p1 <- plot(fit, hyper = TRUE, start = 51)
  # p2 <- plot(fit, pll = FALSE)
  # p3 <- plot(fit)  
  est1 <- summary(fit, hyper = TRUE, recovery = TRUE, ps = pop.mean, type = 1, verbose = TRUE)
  est2 <- summary(fit, hyper = TRUE, recovery = TRUE, ps = pop.scale, type = 2, verbose = TRUE)



