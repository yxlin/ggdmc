cat("\n-------------------- Testing convergence checks --------------------")
rm(list = ls())

## Set up model ----
model <- BuildModel(
  p.map     = list(a = "1", v="1", z="1", d="1", sz="1", sv="1", t0="1",
                   st0="1"),
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
    upper = c(5, 7,  1, 1)
    )

  sigma.prior <- BuildPrior(
    dists = rep("beta", npar),
    p1    = c(a=1, v=1, z=1, t0=1),
    p2    = rep(1, npar),
    upper = rep(1, npar))

## Note the names are important
  priors <- list(pprior=p.prior, location=mu.prior, scale=sigma.prior)

## Fit hierarchical model ----
fit0 <- StartNewsamples(dmi, priors)
fit  <- run(fit0)

isstuck(fit0[[1]])
isstuck(fit[[1]])
isstuck(fit, hyper = TRUE)

isflat(fit[[1]])
ismixed(fit[[1]])
iseffective(fit[[1]], minN = 500, nfun = "mean", FALSE)

ggdmc:::CheckConverged(fit[[1]])


