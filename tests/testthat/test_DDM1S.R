require(ggdmc); require(testthat); require(ggplot2); require(data.table)
context("DDM 1e4")

test_that("DDM 1e4", {
  rm(list = ls())
  setwd("/media/yslin/KIWI/Documents/ggdmc")
  model <- BuildModel(
    p.map     = list(a = "1", v = "1", z = "1", d = "1", sz = "1", sv = "1",
                     t0 = "1", st0 = "1"),
    match.map = list(M = list(s1 = "r1", s2 = "r2")),
    factors   = list(S = c("s1", "s2")),
    responses = c("r1", "r2"),
    constants = c(st0 = 0, d = 0),
    type      = "rd")
  
  p.vector <- c(a = 1, v = 1.2, z = .38, sz = .25, sv = .2, t0 = .15)
  dat <- simulate(model, nsim = 1e4, ps = p.vector)
  dmi <- BuildDMI(dat, model)
  p.prior  <- BuildPrior(
    dists = c(rep("tnorm", 2), "beta", "beta", "tnorm", "beta"),
    p1    = c(a = 1, v = 0, z = 1, sz = 1, sv = 1, t0 = 1),
    p2    = c(a = 1, v = 2, z = 1, sz = 1, sv = 1, t0 = 1),
    lower = c(0, -5, NA, NA, 0, NA),
    upper = c(5,  5, NA, NA, 5, NA))
  plot(p.prior, ps = p.vector)
  
  ## Sampling ---------
  setwd("/media/yslin/KIWI/Documents/ggdmc_lesson")
  path <- c("data/Lesson3/ggdmc_3_7_DDM1S_1e4_tmp.rda")
  ## 6 mins
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
  cat("Done ", path[1], "\n")
  setwd("/media/yslin/KIWI/Documents/ggdmc")
  
  p0 <- plot(fit)
  p1 <- plot(fit, start = 101)
  p2 <- plot(fit, start = 201)
  gridExtra::grid.arrange(p0, p1, p2, ncol = 3)
  d <- data.table(fit$data)
  d[, .N, .(S)]
  ## Analysis -----------
  est <- summary(fit, recovery = TRUE, start = 201, ps = p.vector, verbose = TRUE)
  #                   a     v    z    sz    sv   t0
  # True           1.00  1.20 0.38  0.25  0.20 0.15
  # 2.5% Estimate  1.00  1.16 0.38  0.05  0.01 0.15
  # 50% Estimate   1.00  1.19 0.38  0.18  0.11 0.15
  # 97.5% Estimate 1.01  1.22 0.38  0.24  0.30 0.15
  # Median-True    0.00 -0.01 0.00 -0.07 -0.09 0.00
  #                   a   sv    sz   t0    v    z
  # True           1.00 0.20  0.25 0.15 1.20 0.38
  # 2.5% Estimate  1.00 0.12  0.12 0.15 1.16 0.38
  # 50% Estimate   1.00 0.36  0.20 0.15 1.20 0.38
  # 97.5% Estimate 1.01 0.52  0.25 0.15 1.25 0.39
  # Median-True    0.00 0.16 -0.05 0.00 0.00 0.00
  
  #                   a    sv    sz   t0     v    z
  # True           1.00  0.20  0.25 0.15  1.20 0.38
  # 2.5% Estimate  0.99  0.02  0.11 0.15  1.15 0.38
  # 50% Estimate   1.00  0.18  0.20 0.15  1.18 0.38
  # 97.5% Estimate 1.01  0.40  0.25 0.15  1.21 0.38
  # Median-True    0.00 -0.02 -0.05 0.00 -0.02 0.00
  #                   a    sv    sz   t0     v    z
  # True           1.00  0.20  0.25 0.15  1.20 0.38
  # 2.5% Estimate  0.99  0.01  0.10 0.15  1.14 0.37
  # 50% Estimate   1.00  0.12  0.20 0.15  1.17 0.38
  # 97.5% Estimate 1.01  0.33  0.26 0.15  1.21 0.38
  # Median-True    0.00 -0.08 -0.05 0.00 -0.03 0.00
  #                   a   sv   sz   t0    v    z
  # True           1.00 0.20 0.25 0.15 1.20 0.38
  # 2.5% Estimate  0.99 0.03 0.24 0.15 1.17 0.38
  # 50% Estimate   1.00 0.26 0.30 0.15 1.20 0.38
  # 97.5% Estimate 1.01 0.48 0.33 0.15 1.24 0.38
  # Median-True    0.00 0.06 0.05 0.00 0.00 0.00
  
})

