require(ggdmc); require(data.table); require(ggplot2); require(gridExtra)
require(testthat)

  rm(list = ls())
  setwd("/media/yslin/KIWI/Documents/ggdmc/")
  model <- BuildModel(
    p.map      = list(alpha = "1", beta = "1", epsilon = "1"),
    match.map  = NULL,
    regressors = c(8, 15, 22, 29, 36),
    factors    = list(S = c("x1")),
    responses  = "r1",
    constants  = NULL,
    type       = "glm")
  
  npar <- length(GetPNames(model))
  pop.mean <- c(alpha = 242.7, beta = 6.185, epsilon = 6.086)
  pop.scale <- c(alpha = 100, beta = 5, epsilon = 5)
  
  # pop.mean <- c(alpha = 22.7, beta = 6.185, epsilon = 6.086)
  # pop.scale <- c(alpha = 10, beta = 5, epsilon = 5)
  ntrial <- 1000
  pop.prior  <-BuildPrior(
    dists = rep("tnorm", npar),
    p1    = pop.mean,
    p2    = pop.scale,
    lower = c(NA, 0, 0),
    upper = rep(NA, npar))
  plot(pop.prior, ps = pop.mean)
  dat <- simulate(model, nsub = 500, nsim = ntrial, prior = pop.prior)
  dmi <- BuildDMI(dat, model)
  ps  <- attr(dat, "parameters")
  colMeans(ps)
  matrixStats::colSds(ps)

  # p0 <- ggplot(dat, aes(x = X, y = RT)) +
  #   geom_point() +
  #   theme_bw(base_size = 14) +
  #   facet_wrap(.~s)
  # plot(p0)

  p.prior  <-BuildPrior(
    dists = rep("tnorm", npar),
    p1    = pop.mean, 
    p2    = pop.scale,
    lower = c(NA, 0, 0),
    upper = rep(NA, npar))
  mu.prior  <-BuildPrior(
    dists = rep("tnorm", npar),
    p1    = pop.mean,
    p2    = pop.scale,
    lower = c(NA, 0, 0),
    upper = rep(NA, npar))
  sigma.prior  <-BuildPrior(
    dists = rep("tnorm", npar),
    p1    = c(alpha = 100, beta = 5, epsilon = 5),
    p2    = c(alpha = 50,  beta = 5, epsilon = 5),
    lower = c(0, 0, 0),
    upper = rep(NA, npar))
  pp.prior <- list(mu.prior, sigma.prior)

  ## Sampling random-effect-----------
  setwd("/media/yslin/KIWI/Documents/ggdmc_lesson/")
  path <- c("data/Lesson4/ggdmc_4_0_hglm_random.rda")
  # load(path)
  
  fit0 <- run(StartNewHypersamples(5e2, dmi, p.prior, pp.prior))
  fit  <- fit0
  save(fit0, fit, file = path[1])
  thin <- 1
  repeat {
    fit <- run(RestartHypersamples(5e2, fit, thin = thin))
    save(fit0, fit, file = path[1])
    rhat <- hgelman(fit, verbose = TRUE)
    if (all(rhat < 1.2)) break
    thin <- thin * 2
  }
  cat("Done ", path[1], "\n")
  setwd("/media/yslin/KIWI/Documents/ggdmc/")

  p0 <- plot(fit, hyper = TRUE)
  # p0 <- plot(fit0, hyper = TRUE, start = 301)
  est1 <- summary(fit, hyper = TRUE, recover = TRUE, ps = pop.mean,  type = 1, verbose = TRUE)
  est2 <- summary(fit, hyper = TRUE, recover = TRUE, ps = pop.scale, type = 2, verbose = TRUE)
  # alpha  beta epsilon
  # True           242.70  6.18    6.09
  # 2.5% Estimate  242.47  5.03    5.14
  # 50% Estimate   251.83  5.78    5.88
  # 97.5% Estimate 259.84  6.37    6.50
  # Median-True      9.13 -0.40   -0.21
  
  # alpha beta epsilon
  # True           100.00 5.00    5.00
  # 2.5% Estimate   95.47 4.60    4.57
  # 50% Estimate   101.04 5.04    4.99
  # 97.5% Estimate 107.51 5.59    5.56
  # Median-True      1.04 0.04   -0.01
  