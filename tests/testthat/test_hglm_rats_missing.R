require(ggdmc); require(data.table); require(ggplot2); require(gridExtra)

  rm(list = ls())
  setwd("/media/yslin/KIWI/Documents/Bayesian_rats_GA/")
  tmp <- dget("data/Rats_data.R")
  ## Control groups
  d <- data.frame(matrix(as.vector(tmp$Y), nrow = 30, byrow = TRUE))
  names(d) <- c(8, 15, 22, 29, 36)
  
  d[6:10,5] <- NA
  d[11:20,4:5] <- NA
  d[21:25,3:5] <- NA
  d[26:30,2:5] <- NA
  
  d$s <- factor(1:tmp$N)
  long <- melt(d, id.vars = c("s"), variable.name = "xfac",
               value.name = "RT")
  dplyr::tbl_df(long)
  long$X <- as.double(as.character(long$xfac)) - tmp$xbar
  long$S <- factor("x1")
  long$R <- factor("r1")
  d <- long[, c("s", "S", "R", "X", "RT")]
  
  p1 <- ggplot(d, aes(x = X, y = RT, group = s)) +
    geom_point() + geom_line() +
    theme_bw(base_size = 14) 
  plot(p1)
  dplyr::tbl_df(d)
  
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
  dmi <- BuildDMI(d[!is.na(d$RT),], model)

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
  loc <- c(a = 200, b = 6, tau = 3)
  sca <- c(a = 1e-4, b = 1e-3, tau = 1e-2)
  mu.prior  <- BuildPrior(
    dists = rep("tnorm2", npar),
    p1    = loc,
    p2    = sca,
    lower = c(NA, 0, 0),
    upper = rep(NA, npar))
  sigma.prior <- BuildPrior(
    dists = rep("gamma", npar),
    p1    = c(a = .01, b = .01, tau = .01),
    p2    = c(a = 1000,  b = 1000, tau = 1000),
    lower = c(0, 0, 0),
    upper = rep(NA, npar))
  prior <- list(p.prior, mu.prior, sigma.prior)
  head(StartNewhiersamples)
  fit <- run(StartNewhiersamples(500, dmi, start, prior, thin = 16))
  fit <- run(RestartHypersamples(500, fit, thin = 2))
  
  save(fit, file = "/media/yslin/KIWI/Documents/BUGS_Examples/vol1/Rats/data/ggdmc_missing.RData")
  plot(fit, hyper = TRUE, start = 1)  

  load("/media/yslin/KIWI/Documents/BUGS_Examples/vol1/Rats/data/ggdmc_missing.RData")
  est1 <- summary(fit, hyper = TRUE, type = 1, verbose = TRUE)
  round(est1$quantiles, 3)
  #              mean     sd    2.5%     25%     50%     75%   97.5%  Rhat n.eff
  # alpha0    101.177  3.779  93.780  98.670 101.200 103.700 108.600 1.001 15000
  # alpha.c   245.800  2.794 240.300 243.900 245.800 247.700 251.300 1.003  1100
  # beta.c      6.574  0.147   6.286   6.477   6.572   6.669   6.870 1.003   820
  # sigma       6.161  0.748   4.918   5.638   6.083   6.593   7.865 1.002  3700
  # alpha.tau   0.005  0.002   0.003   0.004   0.005   0.006   0.009 1.002  3200
  # beta.tau    9.540 74.115   1.505   2.676   3.620   5.044  13.601 1.092   260
  
  #            2.5%     25%    50%    75%     97.5%
  # a.h1     241.01  244.35  246.07  247.73  251.09
  # alpha.c  240.30  243.90  245.80  247.70  251.30
  # b.h1      6.362   6.526   6.605   6.688   6.844
  # beta.c    6.286   6.477   6.572   6.669   6.870
  
  # tau.h1    0.01   0.05   0.07   0.08   0.11
  
  # a.h2      0.003  0.004  0.005  0.006  0.009
  # alpha.tau 0.003  0.004  0.005  0.006  0.009
  
  # b.h2     1.640   2.757   3.595   4.479   7.241
  # beta.tau 1.505   2.676   3.620   5.044  13.601

  est3 <- summary(fit)
