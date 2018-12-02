require(ggdmc); require(data.table); require(ggplot2); require(gridExtra)
require(testthat)
context("GLM")

test_that("GLM", {
  rm(list = ls())
  setwd("/media/yslin/KIWI/Documents/ggdmc/")
  model <- BuildModel(
    p.map     = list(a = "1", b = "1", tau = "1"),
    match.map = NULL,
    regressors= c(8, 15, 22, 29, 36),
    factors   = list(S = c("x1")),
    responses = "r1",
    constants = NULL,
    type      = "glm")
  
  npar <- length(GetPNames(model))
  1/6.086^2
  1/100^2
  1/5^2
  pop.mean <- c(a = 242.7, b = 6.185, tau = .027)
  pop.scale <- c(b = 1e-4, b = .04, tau = .04)
  ## 500
  #       alpha beta sigma
  # Mean 247.49 7.16  7.37
  # True 247.55 7.16  7.47
  # Diff   0.06 0.00  0.10
  # Sd    97.18 4.25  4.18
  # True  97.14 4.25  4.28
  # Diff  -0.03 0.00  0.09
  
  #       alpha  beta sigma
  # Mean 249.35  7.55  6.38
  # True 250.08  7.52  6.42
  # Diff   0.73 -0.03  0.04
  # Sd   110.66  4.14  3.16
  # True 110.69  4.15  3.15
  # Diff   0.03  0.01 -0.01
  
  # Summary each participant separately
  # alpha  beta sigma
  # Mean 249.34  7.55  6.38
  # True 254.62  7.24  6.27
  # Diff   5.28 -0.30 -0.11
  # Sd   110.66  4.14  3.16
  # True  93.27  4.73  3.73
  # Diff -17.39  0.58  0.56
  
  #       alpha beta sigma
  # Mean 241.51 6.19  7.23
  # True 241.77 6.19  7.29
  # Diff   0.26 0.00  0.06
  # Sd   100.96 4.60  3.96
  # True 101.28 4.60  3.96
  # Diff   0.32 0.00  0.00
  
  ntrial <- 100
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
  round(colMeans(ps), 2)
  round(matrixStats::colSds(ps), 2)
  p0 <- ggplot(dat, aes(x = X, y = RT, group = s, colour = s)) +
    geom_point() + geom_line() +
    theme_bw(base_size = 14) +
    facet_wrap(.~s)
  ## plot(p0)

  1/150^2
  1/10^2
  p.prior  <-BuildPrior(
    dists = rep("tnorm2", npar),
    p1    = c(a = 200, b = 0, tau = 2),
    p2    = c(a = 1e-5, b = .01, tau = .01),
    lower = c(NA, NA, 0),
    upper = rep(NA, npar))
  # mu.prior  <-BuildPrior(
  #   dists = rep("tnorm", npar),
  #   p1    = c(alpha = 200, beta = 0, sigma = 2),
  #   p2    = c(alpha = 50,  beta = 5, sigma = 5),
  #   lower = c(NA, NA, 0),
  #   upper = rep(NA, npar))
  # sigma.prior  <-BuildPrior(
  #   dists = rep("tnorm", npar),
  #   p1    = c(alpha = 100, beta = 5, sigma = 5),
  #   p2    = c(alpha = 50,  beta = 5, sigma = 5),
  #   lower = c(NA, NA, 0),
  #   upper = rep(NA, npar))
  # pp.prior <- list(mu.prior, sigma.prior)
  plot(p.prior, ps = ps)
  # plot(mu.prior, ps = pop.mean)
  # plot(sigma.prior, ps = pop.scale)

  ## shape, scale
  # dprior(ps[1,], sigma.prior)  
  # dgamma(ps[1,], .001, scale = .001, log = TRUE)

  ## Sampling -----------
  setwd("/media/yslin/KIWI/Documents/ggdmc_lesson/")
  path <- c("data/Lesson4/ggdmc_4_0_hglm_fixed.rda")
  # load(path)
  fit0 <- run(StartManynewsamples(5e2, dmi, p.prior))
  fit  <- fit0
  thin <- 1
  
  repeat {
    fit <- run(RestartManysamples(5e2, fit, thin = thin), pm0 = .05, ncore = 4)
    save(fit0, fit, file = path[1])
    rhat <- gelman(fit, verbose = TRUE)
    if (all(sapply(rhat, function(x) x$mpsrf) < 1.2)) break
    thin <- thin * 2
  }
  cat("Done ", path[1], "\n")
  setwd("/media/yslin/KIWI/Documents/ggdmc/")

  p0 <- plot(fit)
  p0 <- plot(fit)
  est1 <- summary(fit, recover = TRUE, ps = ps, verbose = TRUE)
  colMeans(ps)
  matrixStats::colSds(ps)


})


