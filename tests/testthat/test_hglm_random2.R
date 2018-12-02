require(ggdmc); require(data.table); require(ggplot2); require(gridExtra)
require(testthat)
context("GLM")

test_that("GLM", {
    
  rm(list = ls())
  setwd("/media/yslin/KIWI/Documents/ggdmc/")
  model <- BuildModel(
    p.map      = list(a = "1", b = "1", tau = "1"),
    match.map  = NULL,
    regressors = c(8, 15, 22, 29, 36),
    factors    = list(S = c("x1")),
    responses  = "r1",
    constants  = NULL,
    type       = "glm")
  
  npar <- length(GetPNames(model))
  pop.mean <- c(a = 242.7, b = 6.185, tau = 6.086)
  pop.scale <- c(a = 1e-4, b = .04, tau = .04)
  ntrial <- 500
  pop.prior  <-BuildPrior(
    dists = rep("tnorm2", npar),
    p1    = pop.mean,
    p2    = pop.scale,
    lower = c(NA, 0, 0),
    upper = rep(NA, npar))
  
  plot(pop.prior, ps = pop.mean)
  
  dat <- simulate(model, nsub = 100, nsim = ntrial, prior = pop.prior)
  dmi <- BuildDMI(dat, model)
  ps <- attr(dat, "parameters")
  colMeans(ps)
  matrixStats::colSds(ps)

  p0 <- ggplot(dat, aes(x = X, y = RT)) +
    geom_point() +
    theme_bw(base_size = 14) +
    facet_wrap(.~s)
  plot(p0)

  ## alpha0    106.527  3.644  99.390 104.100 106.500 108.900 113.702 1.001  4400
  ## alpha.c   242.681  2.766 237.100 240.900 242.700 244.500 248.100 1.001  5900
  ## beta.c      6.189  0.108   5.974   6.118   6.189   6.259   6.404 1.001  6700
  ## sigma       6.103  0.472   5.273   5.768   6.071   6.397   7.117 1.001 12000
  ## alpha.tau   0.005  0.001   0.003   0.004   0.005   0.006   0.008 1.001 12000
  ## beta.tau    4.153  1.565   1.952   3.078   3.879   4.922   8.026 1.001  4400
  ## deviance  967.271 14.901 940.800 956.800 966.300 976.600 999.700 1.001  4100

  pstart <- BuildPrior(
    dists = c("tnorm", "tnorm", "tnorm"),
    p1    = c(a = 242,  b = 6.19, tau = .027),
    p2    = c(a = 14, b = .49, tau = .01),
    lower = c(NA, NA, 0),
    upper = rep(NA, npar))
  lstart <- BuildPrior(
    dists = c("tnorm", "tnorm", "tnorm"),
    p1    = c(a = 200,  b = 5, tau = .01),
    p2    = c(a = 50, b = 1, tau = .01),
    lower = c(NA, NA, 0),
    upper = rep(NA, npar))
  sstart <- BuildPrior(
    dists = c("tnorm", "tnorm", "tnorm"),
    p1    = c(a = 10,  b = .5, tau = .01),
    p2    = c(a = 5, b = .1, tau = .01),
    lower = c(NA, NA, 0),
    upper = rep(NA, npar))
  plot(lstart, ps = c(a = 242,  b = 6.19, tau = .027))
  plot(sstart, ps = c(a = 14, b = .49, tau = .01))
  start <- list(pstart, lstart, sstart)

  ## Hierarchical model must set p.prior appropriately; NA indicates getting 
  ## values from upper level; with value indicates setting values as constant.
  p.prior  <- BuildPrior(
    dists = c("tnorm2", "tnorm2", "gamma"),
    p1    = c(a = NA,  b = NA, tau = .001),
    p2    = c(rep(NA, 2), 100),
    lower = c(NA, 0, 0),
    upper = rep(NA, npar))
  mu.prior  <- BuildPrior(
    dists = c("tnorm2", "tnorm2", NA),
    p1    = c(a = 242.7, b = 6.185, tau = NA),
    p2    = c(a = 1e-3, b = 1e-3, tau = NA),
    lower = c(0, 0, NA),
    upper = rep(NA, npar))
  tau.prior <- BuildPrior(
    dists = c("gamma", "gamma", NA),
    p1    = c(a = .01, b = .01, tau = NA),
    p2    = c(a = 100,  b = 100, tau = NA),
    lower = c(0, 0, NA),
    upper = rep(NA, npar))
  print(p.prior)
  print(mu.prior)
  print(tau.prior)
  prior <- list(p.prior, mu.prior, tau.prior)

  fit0 <- ggdmc:::init_newhier_start(500, dmi, start, prior, .001, 4, 9)
  fit <- run(fit0)
  fit <- run(RestartHypersamples(500, fit))
  hyper <- attr(fit, "hyper")
  
  print(lstart)
  print(fit0$"1"$p.prior)
  print(hyper$pp.prior[[1]])
  print(hyper$pp.prior[[2]])
  hyper$p.names
  hyper$n.pars  
  fit0[[1]]$n.pars
  
  
  hyper$phi[[1]][,,500]
  hyper$phi[[2]][,,1]
  
  str(hyper$h_summed_log_prior)
  hyper$h_summed_log_prior[1,]
  hyper$h_log_likelihoods[1,]
  fit[[1]]$theta[,,1]
  fit[[1]]$summed_log_prior[1,]
  fit[[1]]$log_likelihoods[1,]
  

  plot(fit, hyper = TRUE, start = 1)  
  plot(fit, hyper = TRUE, start = 101)  
  plot(fit, hyper = TRUE, start = 201)  

  est1 <- summary(fit, start= 201, recovery = TRUE, ps = pop.mean[c("a", "b")], 
                  hyper = TRUE, type = 1, verbose = TRUE)
  est2 <- summary(fit, start = 201, recovery = TRUE, ps = ps)
  #                     a     b
  # True           242.70  6.18
  # 2.5% Estimate  216.41  3.52
  # 50% Estimate   235.34  5.63
  # 97.5% Estimate 255.58  6.87
  # Median-True     -7.36 -0.55
  #           a    b  tau
  # Mean 236.26 6.82 6.75
  # True 236.26 6.82 6.79
  # Diff  -0.01 0.00 0.03
  # Sd   104.48 3.96 3.57
  # True 104.49 3.96 3.64
  # Diff   0.01 0.00 0.07
})
