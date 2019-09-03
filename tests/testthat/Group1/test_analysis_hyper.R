cat("\n-------------------- Testing analysis_hyper --------------------")

rm(list = ls())

cat("\nBuildModel report: \n")
model <- BuildModel(
  p.map     = list(a = "1", v="1", z="1", d="1", sz="1", sv="1", t0="1", st0="1"),
  match.map = list(M = list(s1 = "r1", s2 = "r2")),
  factors   = list(S = c("s1", "s2")),
  responses = c("r1","r2"),
  constants = c(st0 = 0, d = 0, sv = 0, sz = 0),
  type      = "rd")

npar <- length(GetPNames(model))
pop.mean  <- c(a=2,   v=4, z=0.5, t0=0.3)
pop.scale <- c(a=0.5, v=.5, z=0.1, t0=0.05)
pop.prior <- BuildPrior(
  dists = rep("tnorm", npar),
  p1    = pop.mean,
  p2    = pop.scale,
  lower = c(0,-5,  0, 0),
  upper = c(5, 7,  1, 1))

## Simulate some data
dat <- simulate(model, nsub = 8, nsim = 10, prior = pop.prior)
dmi <- BuildDMI(dat, model)
ps <- attr(dat, "parameters")

p.prior <- BuildPrior(
  dists = rep("tnorm", npar),
  p1    = pop.mean,
  p2    = pop.scale*5,
  lower = c(0,-5, 0, 0),
  upper = c(5, 7, 1, 1))
mu.prior <- BuildPrior(
  dists = rep("tnorm", npar),
  p1    = pop.mean,
  p2    = pop.scale*5,
  lower = c(0,-5,  0, 0),
  upper = c(5, 7,  1, 1))
sigma.prior <- BuildPrior(
  dists = rep("beta", npar),
  p1    = c(a=1, v=1, z=1, t0=1),
  p2    = rep(1, npar),
  upper = rep(1, npar))
priors <- list(pprior=p.prior, location=mu.prior, scale=sigma.prior)

## Fit hierarchical model ----
cat("Starting a new hierarchical model fit: \n")
fit0 <- StartNewsamples(dmi, priors, nmc=100)
fit  <- run(fit0)

cat("Testing six scenarios 'theta2mcmclist': \n")

tmp1 <- theta2mcmclist(fit[[1]])
tmp2 <- theta2mcmclist(fit[[2]], start = 10, end = 90)
tmp3 <- theta2mcmclist(fit[[3]], split = TRUE)
tmp4 <- theta2mcmclist(fit[[4]], subchain = TRUE)
tmp5 <- theta2mcmclist(fit[[5]], subchain = TRUE, nsubchain = 4)
tmp6 <- theta2mcmclist(fit[[6]], thin = 2)

#################################40
## effectiveSize example
#################################40
cat("Testing ten scenarios of 'effectiveSize': \n")
es1 <- effectiveSize_one(fit[[1]], 1, 100, 2, TRUE)
es2 <- effectiveSize_one(fit[[1]], 1, 100, 2, FALSE)
es3 <- effectiveSize_many(fit, 1, 100, TRUE)
es4 <- effectiveSize_many(fit, 1, 100, FALSE)
es5 <- effectiveSize_hyper(fit, 1, 100, 2, TRUE)

es6 <- effectiveSize(fit, TRUE, 1, 100, 2, TRUE)
es7 <- effectiveSize(fit, TRUE, 1, 100, 2, FALSE)
es8 <- effectiveSize(fit, FALSE, 1, 100, 2, TRUE)
es9 <- effectiveSize(fit, FALSE, 1, 100, 2, FALSE)
es10 <- effectiveSize(fit[[1]], FALSE, 1, 100, 2, TRUE)

cat("\nTesting 12 scenarios of 'summary': \n")
est1 <- summary(fit[[1]], FALSE)
est2 <- summary(fit[[1]], FALSE, 1, 100)
est3 <- ggdmc:::summary_one(fit[[1]], 1, 100, c(.025, .5, .975), verbose = TRUE)
est4 <- ggdmc:::summary_one(fit[[1]], 1, 100, c(.025, .5, .975), verbose = F)

est5 <- ggdmc:::summary_many(fit, 1, 100, c(.025, .5, .975), FALSE)
est6 <- ggdmc:::summary_many(fit, 1, 100, c(.025, .5, .975), TRUE)
est7 <- summary(fit)
est8 <- summary(fit, verbose = TRUE)
est9 <- summary(fit, verbose = FALSE)

hest1 <- ggdmc:::summary_hyper(fit, 1, 100, F, F, c(.025, .5, .975), 2, F)
hest2 <- ggdmc:::summary_hyper(fit, 1, 100, F, F, c(.05, .5, .9), 2, F)
hest3 <- summary(fit, TRUE)
