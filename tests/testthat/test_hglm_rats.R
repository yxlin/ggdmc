require(ggdmc); require(data.table); require(ggplot2); require(gridExtra)

  rm(list = ls())
  setwd("/media/yslin/KIWI/Documents/Bayesian_rats_GA/")
  tmp <- dget("data/Rats_data.R")
  ## Control groups
  d <- data.frame(matrix(as.vector(tmp$Y), nrow = 30, byrow = TRUE))
  names(d) <- c(8, 15, 22, 29, 36)
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
  dmi <- BuildDMI(d, model)

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
  
  
  
  fit0 <- ggdmc:::init_newhier_start(500, dmi, start, prior, .001, 8, 9)
  fit <- run(fit0)
  thin <- 2
  fit <- run(RestartHypersamples(500, fit, thin = thin))
  # fit <- run(RestartHypersamples(500, fit, thin = 128))
  plot(fit, hyper = TRUE, start = 1)  
  plot(fit, hyper = TRUE, start = 101)  
  plot(fit, hyper = TRUE, start = 301)  
  

  pop.mdn <- c(a = 242.7, b = 6.189, tau = .02713186)
  pop.sd <- c(a = 14.37, b = 0.5077, tau = .04)
  pop.tau <- c(a = .005, b = 3.879, tau = .04)
  
  est1 <- summary(fit, hyper = TRUE, recover = TRUE, start = 101, 
                  ps = pop.mean,  type = 1, verbose = TRUE)
  est2 <- summary(fit, hyper = TRUE, recover = TRUE, start = 101, 
                  ps = pop.tau, type = 2, verbose = TRUE)
  # est3 <- summary(fit, verbose = TRUE)
  est3 <- summary(fit, verbose = FALSE)

  1/6.071^2
  #                 mean se_mean    sd   2.5%    25%    50%    75%  97.5% n_eff
  # alpha0        106.40    0.02  3.64  99.23 104.00 106.40 108.83 113.53 39013
  # mu_alpha      242.47    0.01  2.73 237.09 240.68 242.47 244.29 247.84 38781
  # mu_beta         6.18    0.00  0.11   5.97   6.11   6.18   6.26   6.40 40016
  # sigma_y         6.09    0.00  0.47   5.27   5.77   6.07   6.39   7.09 18710
  # sigma_alpha    14.59    0.01  2.05  11.20  13.15  14.37  15.79  19.14 37148
  # sigma_beta      0.52    0.00  0.09   0.36   0.45   0.51   0.57   0.72 28120
  # sigmasq_alpha 217.16    0.34 63.19 125.39 172.97 206.54 249.38 366.48 35391
  
  ##              mean     sd    2.5%     25%     50%     75%   97.5%  Rhat n.eff
  ## alpha0    106.527  3.644  99.390 104.100 106.500 108.900 113.702 1.001  4400
  ## alpha.c   242.681  2.766 237.100 240.900 242.700 244.500 248.100 1.001  5900
  ## beta.c      6.189  0.108   5.974   6.118   6.189   6.259   6.404 1.001  6700
  ## sigma       6.103  0.472   5.273   5.768   6.071   6.397   7.117 1.001 12000
  ## alpha.tau   0.005  0.001   0.003   0.004   0.005   0.006   0.008 1.001 12000
  ## beta.tau    4.153  1.565   1.952   3.078   3.879   4.922   8.026 1.001  4400
  ## deviance  967.271 14.901 940.800 956.800 966.300 976.600 999.700 1.001  4100
