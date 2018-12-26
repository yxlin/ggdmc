require(ggdmc); require(data.table); require(ggplot2); require(gridExtra)
require(testthat)
context("logistic")

test_that("logistic", {
  rm(list = ls())
  setwd("/media/yslin/KIWI/Documents/ggdmc/")
  model <- BuildModel(
    p.map     = list(a0 = "1", a1 = "1", a2 = "1", a3 = "1", b = "1"),
    match.map = NULL,
    regressors= NULL,
    factors   = list(S = c("s1", "s2"), E = c("e1", "e2")),
    responses = "r1",
    constants = NULL,
    type      = "logit")

  p.vector <- c(a0 = -.55, a1 = .08, a2 = -.81, a3 = 1.35, b = .27)
  ntrial <- 500
  dat <- simulate(model, nsim = ntrial, ps = p.vector)
  dplyr::tbl_df(dat)
  dmi <- BuildDMI(dat, model)

  # datnames <- colnames(dat)
  # Xnames   <- datnames[!datnames %in% c("R", "N", "Y")]
  # Ynames   <- c("Y", "N")
  # X_ <- data.matrix( dat[,Xnames])
  # Y <- data.matrix(dat[,Ynames])
  # X <- cbind( rep(1, nrow(dat)), X_, X_[,1]*X_[,2])
  # linpred <- X %*% p.vector[1:4]
  #
  # prob <- Y[,1] / Y[,2]
  # yi <- qlogis(prob, 0, 1, TRUE, FALSE)
  # res0 <- dnorm(yi, linpred, p.vector[5])
  #
  # # prob <- plogis( linpred, 0, 1, 1, 0 );
  # # res0 <- dbinom(Y[,1], Y[,2], prob, FALSE)
  #
  # res1 <- ggdmc:::density_logit(p.vector, X, Y)
  # res2 <- ggdmc:::density_logit(p.vector, X, Y)
  # all(res1==res2)
  # sum(is.na(res1))
  # sum(is.na(res2))
  # all.equal(res1, res2)
  # all.equal(res0, res2[,1])
  # cbind(res0, res1, res2)

  npar <- length(GetPNames(model))
  start  <- BuildPrior(
    dists = c(rep("tnorm2", npar)),
    p1    = c(a0 = -.55, a1 = 0.08, a2 = -.81, a3 = 1.35, b = .267),
    p2    = rep(1e-6, npar),
    lower = rep(NA, npar),
    upper = rep(NA, npar))

  p.prior  <-BuildPrior(
    dists = c(rep("tnorm2", npar)),
    p1    = c(a0 = 0, a1 = 0, a2 = 0, a3 = 0, b = 0),
    p2    = rep(1e-6, npar),
    lower = rep(NA, npar),
    upper = rep(NA, npar))
  plot(p.prior, ps = p.vector)
  print(p.prior)
  ## Sampling -----------
  setwd("/media/yslin/KIWI/Documents/ggdmc_lesson/")
  path <- c("data/Lesson3/ggdmc_3.8_logit.rda")

  fit0 <- Startlogit(5e2, dmi, start, p.prior, thin = 2)
  fit <- run(fit0)

  thin <- 2
  repeat {
    fit <- run(RestartSamples(5e2, fit, thin = thin))
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
  ## plot(fit, pll =F, den=T)

  est <- summary(fit, start = 101, recover = TRUE, ps = p.vector, verbose = TRUE)
  #                   b0   b1    b2   b3       sd
  # True           -0.55 0.08 -0.81 1.35     0.27
  # 2.5% Estimate  -0.62 0.00 -0.97 1.17 -1742.63
  # 50% Estimate   -0.53 0.12 -0.83 1.35     0.97
  # 97.5% Estimate -0.44 0.25 -0.69 1.54  1660.99
  # Median-True     0.02 0.04 -0.02 0.00     0.70

  ## Gelman's simulation method
  thetas <- matrix(aperm(fit$theta, c(3,1,2)), ncol = fit$n.par)
  colnames(thetas) <- fit$pnames
  use <- sample(1:n.sims, n.sims, replace = FALSE)
  theta_rnd <- thetas[use,]

  n.sims <- 1000  ## 1e3 simulation study; each study has 1e3 trials
  X.tilde <- cbind (1, dat$X)
  ntrial <- nrow (X.tilde)
  y.tilde <- array (NA, c(n.sims, ntrial))
  for (s in 1:n.sims) {
    mu <- X.tilde %*% theta_rnd[s, 1:2]
    sd <- 1/sqrt(theta_rnd[s, 3])
    y.tilde[s,] <- rnorm(ntrial, mu, sd)
  }



  p1color <- "#FFFFFF"
  p2color <- "#A9A9A9"
  bin1 <- seq(min(dat$RT) - .1, max(dat$RT) + .1, .3)
  p1 <- hist(dat$RT, breaks = "fd", plot = FALSE)
  plot(p1, col = p1color, main = "")

  for(s in 1:n.sims) {
    p2 <- hist(y.tilde[s,], breaks = "fd", plot = FALSE)
    lines(p2$mids, p2$counts, col = p2color)
  }

  p3 <- hist(dat$RT, breaks = "fd", plot = FALSE)
  lines(p3$mids, p3$counts, col = "black", lwd = 3)

})


