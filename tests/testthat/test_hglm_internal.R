require(ggdmc); require(data.table); require(ggplot2); require(gridExtra)
require(testthat)
context("GLM")

test_that("GLM", {
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
  
  ntrial <- 100
  pop.prior  <-BuildPrior(
    dists = rep("tnorm", npar),
    p1    = pop.mean,
    p2    = pop.scale,
    lower = c(NA, 0, 0),
    upper = rep(NA, npar))
  plot(pop.prior, ps = pop.mean)
  dat <- simulate(model, nsub = 30, nsim = ntrial, prior = pop.prior)
  dmi <- BuildDMI(dat, model)
  ps <- attr(dat, "parameters")
  colMeans(ps)

  p0 <- ggplot(dat, aes(x = X, y = RT)) +
    geom_point() +
    theme_bw(base_size = 14) +
    facet_wrap(.~s)
  plot(p0)

  npar <- length(GetPNames(model))
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
  plot(p.prior, ps = ps)
  plot(mu.prior, ps = pop.mean)
  plot(sigma.prior, ps = pop.scale)
  ## shape, scale
  # dprior(ps[1,], sigma.prior)  
  # dgamma(ps[1,], .001, scale = .001, log = TRUE)

  ## Sampling random-effect-----------
  setwd("/media/yslin/KIWI/Documents/ggdmc_lesson/")
  path <- c("data/Lesson4/ggdmc_4_0_hglm_tmp.rda")
  # load(path)
  
  fit0 <- run(StartNewHypersamples(5e2, dmi, p.prior, pp.prior))
  fit  <- fit0
  thin <- 2
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
  # p0 <- plot(fit[[17]], den = TRUE, pll = FALSE)
  est1 <- summary(fit, hyper = TRUE, recover = TRUE, ps = pop.mean,  type = 1, verbose = TRUE)
  est2 <- summary(fit, hyper = TRUE, recover = TRUE, ps = pop.scale, type = 2, verbose = TRUE)
  #                 alpha  beta epsilon
  # True            22.70  6.18    6.09
  # 2.5% Estimate   19.11  4.28    2.59
  # 50% Estimate    22.29  7.34    6.40
  # 97.5% Estimate  25.44  9.12    8.28
  # Median-True     -0.41  1.15    0.32
  # 2.5% Estimate  192.79  0.23    3.22
  # 50% Estimate   227.83  3.36    6.87
  # 97.5% Estimate 261.23  6.86    8.49
  # Median-True    -14.87 -2.82    0.78
  #                 alpha   beta epsilon
  # True            10.00   5.00    5.00
  # 2.5% Estimate    7.19   3.44    3.48
  # 50% Estimate     9.18   4.64    4.93
  # 97.5% Estimate  12.26   7.11    7.82
  # Median-True     -0.82  -0.36   -0.07
  # 2.5% Estimate   83.32   5.04    3.15
  # 50% Estimate   103.52   7.04    4.34
  # 97.5% Estimate 134.60   9.86    7.23
  # Median-True      3.52   2.04   -0.66
  
  fit.lm <- lm(RT ~ X, data = dat)
  summary(fit.lm)
  coef(fit.lm)
  sd(fit.lm$residuals)

  ## Internal ------------
  (pnames   <- names(attr(model, "p.vector")))
  (allpar   <- attr(model, "all.par"))
  (parnames <- attr(model, "par.names"))
  (type     <- attr(model, "type"))
  (n1       <- attr(model, "n1.order"))
  (resp     <- attr(model, "responses"))
  (cell     <- ggdmc:::check_cell(1, model))
  (isr1     <- ggdmc:::check_rd(type, model))
  (dim1 <- dimnames(model)[[1]])
  (dim2 <- dimnames(model)[[2]])
  (dim3 <- dimnames(model)[[3]])
  (res <- ggdmc:::p_df(p.vector, cell, pnames, allpar, parnames, model, type, 
                      dim1, dim2, dim3, isr1, n1, TRUE))
  TableParameters(p.vector, 1, model, FALSE)
  
  ise <- attr(dmi, "cell.empty")
  cellidx  <- ggdmc:::cellIdx2Mat(dmi)

  ## cell.empty internal ----------  
  res <- ggdmc:::check_BuildDMI(dat, model)
  subject_models <- res$issm
  modeli <- res$model
  (fnams <- names(attr(modeli, "factors")))
  tmp1 <- dat[, c(fnams, "R")]
  head(dat)
  head(tmp1)
  
  cells <- apply(dat[, c(fnams, "R")], 1, paste, collapse = ".")
  ## ncell == nobservation
  
  sapply(dat[, c(fnams, "R")], levels)
  cell.index <- vector("list", dim(model)[1])
  names(cell.index) <- row.names(model)
  ## scan trial-by-trial (every observation)
  for ( j in names(cell.index) ) cell.index[[j]] <- cells %in% j
  
  ## density_glm
  X <- dat$X
  Y <- dat$RT 
  nsim  <- attr(dmi, "n.pda")
  bw    <- attr(dmi, "bw")
  debug <- attr(dmi, "debug")
  posdrift <- attr(model, "posdrift")
  n1idx    <- attr(model, "n1.order")
  mc       <- attr(model, "match.cell")

  den1 <- ggdmc:::density_glm(p.vector, pnames, allpar, parnames, model, type,
                              dim1, dim2, dim3, ise, cellidx, X, Y)
  
  den2 <- dnorm(Y, mu, p.vector[3], FALSE)
  den3 <- ggdmc:::sumloglike_glm(p.vector, pnames, allpar, parnames, model, type,
                                 dim1, dim2, dim3, n1idx, ise, cellidx, dmi$RT, mc,
                                 isr1, X, posdrift, nsim, bw, 1, 0, debug)
  
  all.equal(den1[,1], den2)
  sum(log(den1)) == sum(log(den2))
  all.equal( sum(log(den1)), den3)


})


