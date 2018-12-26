require(ggdmc); require(data.table); require(ggplot2); require(gridExtra)
require(testthat)
context("hglm")

test_that("hglm", {
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
  pop.mean <- c(a = 242.7, b = 6.189, tau = .03)
  pop.scale <- c(a = .005, b = 3.879, tau = .04)
  ntrial <- 100
  pop.prior  <-BuildPrior(
    dists = rep("tnorm2", npar),
    p1    = pop.mean,
    p2    = pop.scale,
    lower = c(NA, 0, 0),
    upper = rep(NA, npar))
  
  # png("~/myblog/images/BUGS/pop_prior.png", 800, 600)
  plot(pop.prior, ps = pop.mean)
  # dev.off()
  dat <- simulate(model, nsub = 30, nsim = ntrial, prior = pop.prior)
  dmi <- BuildDMI(dat, model)
  ps <- attr(dat, "parameters")
  # p0 <- ggplot(dat, aes(x = X, y = RT, group = s, color = s)) +
  #   geom_point() + geom_line() + 
  #   ylab("Y") +
  #   theme_bw(base_size = 14) 
  #   # facet_wrap(.~s)
  # plot(p0)
  
  pstart <- BuildPrior(
    dists = c("tnorm", "tnorm", "tnorm"),
    p1    = c(a = 242,  b = 6.19, tau = .027),
    p2    = c(a = 14, b = .49, tau = 10),
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
  start <- list(pstart, lstart, sstart)
  
  p.prior  <- BuildPrior(
    dists = rep("tnorm2", npar),
    p1    = c(a = NA, b = NA, tau = NA), 
    p2    = rep(NA, 3),
    lower = c(NA, 0, 0),
    upper = rep(NA, npar))
  mu.prior  <- BuildPrior(
    dists = rep("tnorm2", npar),
    p1    = c(a = 0, b = 0, tau = 0),
    p2    = c(a = 1e-4, b = 1e-4, tau = 1e-4),
    lower = c(NA, 0, 0),
    upper = rep(NA, npar))
  sigma.prior <- BuildPrior(
    dists = rep("gamma", npar),
    p1    = c(a = .001, b = .001, tau = .001),
    p2    = c(a = 1000,  b = 1000, tau = 1000),
    lower = c(0, 0, 0),
    upper = rep(NA, npar))
  prior <- list(p.prior, mu.prior, sigma.prior)

  # DT <- data.table(d)
  # DT[, .N, .(s)]
  
  fit0 <- StartNewhiersamples(500, dmi, start, prior, thin = 8)
  fit <- run(fit0, hpm = 0)
  thin <- 2
  # repeat {
    fit <- run(RestartHypersamples(5e2, fit, thin = thin))
    save(fit0, fit, file = path[1])
    rhat <- hgelman(fit, verbose = TRUE)
    if (all(rhat < 1.2)) break
    thin <- thin * 2
  # }
  cat("Done ", path[1], "\n")
  setwd("/media/yslin/KIWI/Documents/ggdmc/")

  p0 <- plot(fit, hyper = TRUE)
  p0 <- plot(fit, hyper = TRUE, start = 201)
  p1 <- plot(fit, hyper = TRUE, pll = FALSE, den = TRUE)
  

  est1 <- summary(fit, hyper = TRUE, recover = TRUE,  
                  ps = pop.mean,  type = 1, verbose = TRUE, digits = 3)
  est2 <- summary(fit, hyper = TRUE, recover = TRUE, 
                  ps = pop.scale, type = 2, verbose = TRUE, digits = 3)
  est3 <- summary(fit, recover = TRUE, ps = ps, verbose = TRUE)
  
  hes <- effectiveSize(fit, hyper = TRUE)
  es <- effectiveSize(fit)
  round(apply(data.frame(es), 1, mean))
  round(apply(data.frame(es), 1, sd))
  round(apply(data.frame(es), 1, max))
  round(apply(data.frame(es), 1, min))
  
  
  ## est4 <- summary(fit)

})


