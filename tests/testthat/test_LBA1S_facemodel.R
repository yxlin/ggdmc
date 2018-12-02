require(ggdmc); require(testthat); require(ggplot2); require(data.table)

  rm(list = ls())
  setwd("/media/yslin/KIWI/Documents/ggdmc")
  ## (P)eople factor, p1: young+female, p2: young+male, p3: old-female, p4: old-male
  ## (S)timuli factor, s0: 0% happiness, s1: 10%, s2: 20%, s3: 30%, s4: 40%,
  ## s5: 50%, s6: 60%, s7: 70%, s8: 80%, s9: 90%, st: 100%
  ## (R)esponse, r1: sad face, r2: happy face
  model <- BuildModel(
    p.map     = list(A = "1", B = "1", t0 = "1", mean_v = c("S", "M"), 
                     sd_v = "1", st0 = "1"),
    match.map = list(M = list(s0 = 1, s1 = 1, s2 = 1, s3 = 1, s4 = 1, s5 = 1,
                              s6 = 2, s7 = 2, s8 = 2, s9 = 2, st = 2)),
    factors   = list(S = c("s0", "s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8",
                           "s9", "st")),
    constants = c(st0 = 0, sd_v = 1),
    responses = c("r1", "r2"),
    type      = "norm")
  pnames <- GetPNames(model)

  p.vector <- c(A = .75, B = 1.25, t0 = .15, seq(1.5, 4.5, length.out = 11),
                round(runif(11, 0, 2), 2))
  names(p.vector) <- pnames
  ntrial <- 600
  dat <- simulate(model, nsim = ntrial, ps = p.vector)
  dmi <- BuildDMI(dat, model)
dplyr::tbl_df(dmi)

  p.prior <- BuildPrior(
    dists = c("tnorm", "tnorm", "beta", rep("tnorm", 22)),
    p1    = rep(1, 25),
    p2    = rep(5, 25),
    lower = c(rep(0, 3),  rep(NA, 22)),
    upper = c(rep(NA, 2), 1, rep(NA, 22)))
  names(p.prior) <- pnames
  plot(p.prior, ps = p.vector)

  ## Sampling ---------
  ## fit0 <- StartNewsamples(3, dmi, p.prior)
  fit0 <- run(StartNewsamples(5e2, dmi, p.prior, thin = 2))
  fit  <- fit0
  path <- "data/fit.RData"
  thin <- 4
  repeat {
    fit <- run(RestartSamples(5e2, fit, thin = thin))
    save(fit0, fit, file = path[1])
    rhat <- gelman(fit, verbose = TRUE)
    if (all(rhat$mpsrf < 1.1)) break
    thin <- thin * 2
  }
  cat("Done ", path[1], "\n")
  p0 <- plot(fit)
  p1 <- plot(fit, start = 101)
  p2 <- plot(fit0, start = 201)
  
  # gridExtra::grid.arrange(p0, p1, p2, ncol = 3)
  # d <- data.table(fit$data)
  # d[, .N, .(S)]
  
  ## Analysis -----------
  est <- summary(fit, start = 101, recovery = TRUE, ps = p.vector)
  ## str(est)
  round(est[,pnames], 2)
##                   A     B   t0 mean_v.s0.true mean_v.s1.true mean_v.s2.true
## True           0.75  1.25 0.15           1.50           1.80           2.10
## 2.5% Estimate  0.71  0.96 0.15           1.33           1.59           1.92
## 50% Estimate   0.86  1.09 0.18           1.47           1.72           2.05
## 97.5% Estimate 0.98  1.25 0.20           1.60           1.86           2.18
## Median-True    0.11 -0.16 0.03          -0.03          -0.08          -0.05
##                mean_v.s3.true mean_v.s4.true mean_v.s5.true mean_v.s6.true
## True                     2.40           2.70           3.00           3.30
## 2.5% Estimate            2.20           2.53           2.86           3.29
## 50% Estimate             2.33           2.66           2.99           3.42
## 97.5% Estimate           2.45           2.78           3.11           3.56
## Median-True             -0.07          -0.04          -0.01           0.12
##                mean_v.s7.true mean_v.s8.true mean_v.s9.true mean_v.st.true
## True                     3.60           3.90           4.20           4.50
## 2.5% Estimate            3.49           3.69           4.17           4.45
## 50% Estimate             3.63           3.83           4.33           4.62
## 97.5% Estimate           3.76           3.97           4.50           4.80
## Median-True              0.03          -0.07           0.13           0.12
##                mean_v.s0.false mean_v.s1.false mean_v.s2.false mean_v.s3.false
## True                      0.89            1.11            1.67            0.80
## 2.5% Estimate             0.74            0.85            1.46            0.33
## 50% Estimate              0.90            1.01            1.60            0.54
## 97.5% Estimate            1.06            1.17            1.75            0.75
## Median-True               0.01           -0.10           -0.07           -0.26
##                mean_v.s4.false mean_v.s5.false mean_v.s6.false mean_v.s7.false
## True                      1.10            0.26            1.15            0.99
## 2.5% Estimate             1.01           -0.48            0.68            0.68
## 50% Estimate              1.20           -0.12            0.92            0.95
## 97.5% Estimate            1.37            0.22            1.16            1.19
## Median-True               0.10           -0.38           -0.23           -0.04
##                mean_v.s8.false mean_v.s9.false mean_v.st.false
## True                      1.67            0.47            1.06
## 2.5% Estimate             1.36           -1.52            0.60
## 50% Estimate              1.57           -0.40            1.00
## 97.5% Estimate            1.78            0.39            1.36
## Median-True              -0.10           -0.87           -0.06

  
