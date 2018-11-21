context("LBA 1e4")
require(ggdmc); require(testthat); require(ggplot2); require(data.table)

test_that("LBA 1e4", {
  rm(list = ls())
  setwd("/media/yslin/KIWI/Documents/ggdmc")
  model <- BuildModel(
    p.map     = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1",
                     st0 = "1"),
    match.map = list(M = list(s1 = 1, s2 = 2)),
    factors   = list(S = c("s1", "s2")),
    constants = c(st0 = 0, sd_v = 1),
    responses = c("r1", "r2"),
    type      = "norm")

  mapinfo <- ggdmc:::check_BuildModel(p.map, responses, factors, match.map,
                                      constants, type)
  
  mapinfo <- ggdmc:::check_BuildModel(
    p_map     = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1",
                     st0 = "1"),
    responses = c("r1", "r2"),
    factors   = list(S = c("s1", "s2")),
    match_map = list(M = list(s1 = 1, s2 = 2)),
    constants = c(st0 = 0, sd_v = 1),
    type      = "norm")
  
  
  p.vector <- c(A = .75, B = 1.25, t0 = .15, mean_v.true = 2.5, mean_v.false = 1.5)
  ntrial <- 1e4
  dat <- simulate(model, nsim = ntrial, ps = p.vector)
  dmi <- BuildDMI(dat, model)
  p.prior <- BuildPrior(
    dists = c("tnorm", "tnorm", "beta", "tnorm", "tnorm"),
    p1    = c(A = 1, B = 1, t0 = 1, mean_v.true = 1, mean_v.false = 1),
    p2    = c(1,  1,  1, 1, 1),
    lower = c(rep(0, 3),  rep(NA, 2)),
    upper = c(rep(NA, 2), 1, rep(NA, 2)))
  plot(p.prior, ps = p.vector)

  ## Sampling ---------
  setwd("/media/yslin/KIWI/Documents/ggdmc_lesson")
  path <- c("data/Lesson3/ggdmc_3_2_LBA1S_1e4_tmp.rda")
  fit0 <- run(StartNewsamples(5e2, dmi, p.prior))
  fit  <- fit0
  thin <- 1
  repeat {
    fit <- run(RestartSamples(5e2, fit, thin = thin))
    save(fit0, fit, file = path[1])
    rhat <- gelman(fit, verbose = TRUE)
    if (all(rhat$mpsrf < 1.1)) break
    thin <- thin * 2
  }
  rhat <- gelman(fit, verbose = TRUE, start = 101)
  cat("Done ", path[1], "\n")
  p0 <- plot(fit0)
  p1 <- plot(fit0, start = 101)
  p2 <- plot(fit0, start = 201)
  gridExtra::grid.arrange(p0, p1, p2, ncol = 3)
  d <- data.table(fit$data)
  d[, .N, .(S)]
  ## Analysis -----------
  est <- summary(fit, start = 201, recovery = TRUE, ps = p.vector, verbose = TRUE)
  ##                   A    B   t0 mean_v.true mean_v.false
  ## True           0.75 0.25 0.20        2.50         1.50
  ## 2.5% Estimate  0.60 0.17 0.17        2.24         1.12
  ## 50% Estimate   0.75 0.25 0.20        2.59         1.48
  ## 97.5% Estimate 0.90 0.35 0.22        2.95         1.84
  ## Median-True    0.00 0.00 0.00        0.09        -0.02
  
  ##                    A     B mean_v.false mean_v.true   t0
  ## True            0.75  1.25         1.50        2.50 0.15
  ## 2.5% Estimate   0.56  1.08         1.37        2.40 0.12
  ## 50% Estimate    0.72  1.23         1.46        2.47 0.15
  ## 97.5% Estimate  0.84  1.41         1.55        2.53 0.18
  ## Median-True    -0.03 -0.02        -0.04       -0.03 0.00
  
  ## problem?
  #                    A    B mean_v.false mean_v.true    t0
  # True            0.75 1.25         1.50        2.50  0.15
  # 2.5% Estimate   0.10 1.07         1.40       -0.48  0.08
  # 50% Estimate    0.58 1.35         1.48        1.13  0.13
  # 97.5% Estimate  0.85 1.70         1.57        2.71  0.19
  # Median-True    -0.17 0.10        -0.02       -1.37 -0.02
  #                   A     B mean_v.false mean_v.true   t0
  # True           0.75  1.25         1.50        2.50 0.15
  # 2.5% Estimate  0.67  0.89         1.38       -0.54 0.15
  # 50% Estimate   0.93  1.04         1.45        1.03 0.20
  # 97.5% Estimate 1.07  1.28         1.53        2.66 0.23
  # Median-True    0.18 -0.21        -0.05       -1.47 0.05
})



