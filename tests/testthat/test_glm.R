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

  # as.vector(level_array)
  # length(level_array)
  # str(col.par.levels)
  # col.fac
  1/sqrt(.01)
  p.vector <- c(a = 242.7, b = 6.185, tau = .01)
  ntrial <- 1000
  dat <- simulate(model, nsim = ntrial, ps = p.vector)

  # x <- rep(c(8, 15, 22, 29, 36), each = ntrial/5)
  # mu <- p.vector[1] + p.vector[2] * x
  # y  <- rnorm(ntrial, mu, p.vector[3])
  # dat <- data.frame(S = rep("x1", ntrial), R = rep("r1", ntrial), X = x, RT = y)
  dplyr::tbl_df(dat)
  dmi <- BuildDMI(dat, model)

  npar <- length(GetPNames(model))
  p.prior  <-BuildPrior(
    dists = c("tnorm2", "tnorm2", "gamma"),
    p1    = c(a = 200, b = 0, tau = .1),
    p2    = c(a = 1e-6, b = 1e-6, tau = .1),
    lower = c(NA, NA, NA),
    upper = rep(NA, npar))
  plot(p.prior, ps = p.vector)
  plot(p.prior)
  print(p.prior)
  ## Sampling -----------
  setwd("/media/yslin/KIWI/Documents/ggdmc_lesson/")
  path <- c("data/Lesson3/ggdmc_3.0_glm.rda")
  
  fit0 <- run(StartNewsamples(5e2, dmi, p.prior), pm0 = .05)
  fit  <- fit0
  thin <- 1
  repeat {
    fit <- run(RestartSamples(5e2, fit, thin = thin), pm = .05)
    save(fit0, fit, file = path[1])
    rhat <- gelman(fit, verbose = TRUE)
    if (all(rhat$mpsrf < 1.1)) break
    thin <- thin * 2
  }
  cat("Done ", path[1], "\n")
  setwd("/media/yslin/KIWI/Documents/ggdmc/")
  
  ## Analysis -----------
  p0 <- plot(fit)
  p1 <- plot(fit, start = 101)
  p1 <- plot(fit)
  
  est <- summary(fit, start = 101, recover = TRUE, ps = p.vector, verbose = TRUE)
  
  #                 alpha beta epsilon
  # True           242.70 6.18    6.09
  # 2.5% Estimate  241.66 6.17    5.76
  # 50% Estimate   242.44 6.20    5.98
  # 97.5% Estimate 243.19 6.24    6.22
  # Median-True     -0.26 0.02   -0.10
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

  mu <- p.vector[1] + p.vector[2] * X
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


