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
  # sigma_y         6.09    0.00  0.47   5.27   5.77   6.07   6.39   7.09 18710
  
  pstart <- BuildPrior(
    dists = c("tnorm", "tnorm", "constant"),
    # p1    = c(a = 242,  b = 6.19, tau = 6.09),
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
  plot(lstart, ps = c(a = 242,  b = 6.19, tau = .027))
  plot(sstart, ps = c(a = 14, b = .49, tau = .01))
  start <- list(pstart, lstart, sstart)
  
  p.prior  <- BuildPrior(
    dists = c("tnorm2", "tnorm2", "tnorm"),
    p1    = c(a = NA,  b = NA, tau = NA),
    p2    = rep(NA, 3),
    lower = c(NA, 0, 0),
    upper = rep(NA, npar))
  mu.prior  <- BuildPrior(
    dists = c("tnorm2", "tnorm2", "gamma"),
    p1    = c(a = 0, b = 0, tau = .01),
    p2    = c(a = 1e-3, b = 1e-3, tau = 100),
    lower = c(0, 0, NA),
    upper = rep(NA, npar))
  tau.prior <- BuildPrior(
    dists = c("gamma", "gamma", "gamma"),
    p1    = c(a = .01, b = .01, tau = .01),
    p2    = c(a = 100,  b = 100, tau = 100),
    lower = c(0, 0, 0),
    upper = rep(NA, npar))
  print(p.prior)
  print(mu.prior)
  print(tau.prior)
  prior <- list(p.prior, mu.prior, tau.prior)
  
  fit0 <- ggdmc:::init_newhier_start(500, dmi, start, prior, .001, 8, 9)
  
  fit <- fit0
  hyper <- attr(fit, 'hyper')
  hyper$phi[[1]][,,1]
  hyper$phi[[2]][,,1]
  # hyper$h_summed_log_prior[1,]
  # hyper$h_log_likelihoods[1,]
  fit[[1]]$theta[,,1]
  fit[[2]]$theta[,,1]
  fit[[3]]$theta[,,1]
  
  fit[[1]]$theta[,,500]
  fit[[2]]$theta[,,500]
  fit[[3]]$theta[,,500]
  
  fit[[1]]$summed_log_prior[1,]
  fit[[1]]$log_likelihoods[1,]
  fit[[2]]$summed_log_prior[1,]
  fit[[2]]$log_likelihoods[1,]
  
  hyper$n.pars
  
  fit <- run(fit0)
  fit <- run(RestartHypersamples(500, fit, thin = 8))
  fit <- run(RestartHypersamples(500, fit, thin = 32))
  plot(fit, hyper = TRUE, start = 1)  
  plot(fit, hyper = TRUE, start = 101)  
  plot(fit, hyper = TRUE, start = 301)  
  
  1/sqrt(.056)
  1/sqrt(256.9)
  
  # plot(fit, hyper = FALSE, pll= FALSE, start = 1)  

  pop.mean <- c(a = 242.7, b = 6.185, tau = .02714)
  pop.scale <- c(a = 14.14, b = 0.5077, tau = .04)
  pop.tau <- c(a = .005, b = 3.879, tau = .04)
  est1 <- summary(fit, recovery = TRUE, ps = pop.mean, start = 301,
                  hyper = TRUE, type = 1, verbose = TRUE, digits = 4)
  truesd <- c(7.09, 6.39, 6.07, 5.77, 5.27)
  round(1/truesd^2, 3)
  
  est2 <- summary(fit, recovery = TRUE, ps = pop.tau, start = 101,
                  hyper = TRUE, type = 2, verbose = TRUE, digits = 4)
  
  est3 <- summary(fit)
  
  

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
  

  ## Sampling random-effect-----------
  fit0 <- ggdmc:::init_newhier_hstart(500, dmi, p.prior, pp.prior, .001, 2, 9,
                                      pop.mean, pop.scale)
  fit <- fit0
  class(fit0) <- "model"
  
  fit0[[1]]$theta[,,1]
  p0 <- plot(fit, hyper = TRUE)
  p0 <- plot(fit, hyper = TRUE, start = 101)
  
  
  # load(path)
  fit0 <- run(StartNewHypersamples(5e2, dmi, p.prior, pp.prior))
  fit  <- fit0
  thin <- 16
  repeat {
    fit <- run(RestartHypersamples(5e2, fit, thin = thin), hpm = .15, pm = .15)
    save(fit0, fit, file = path[1])
    rhat <- hgelman(fit, verbose = TRUE)
    if (all(rhat < 1.2)) break
    thin <- thin * 2
  }
  cat("Done ", path[1], "\n")
  setwd("/media/yslin/KIWI/Documents/ggdmc/")

  pop.mean <- c(1, 2, NA)
  pop.scale <- c(1:3)
  names(pop.mean) <- c('a', 'b', 'tau')
  # p0 <- plot(fit[[17]], den = TRUE, pll = FALSE)
  est1 <- summary(fit, hyper = TRUE, type = 1, verbose = F)
  round(est1$quantiles[,3], 2)
  # a.h1   b.h1   a.h2   b.h2 
  # 244.66   5.61   0.01   5.62 
  est2 <- summary(fit, verbose = FALSE)
  
  # a        b         tau
  # 1    243.1048 5.550598 0.120267304
  # 2    251.0300 6.240905 0.019611783
  # 3    254.5366 5.850082 0.011873713
  # 4    235.2547 4.796330 0.115399564
  # 5    235.1550 6.041915 0.071275270
  # 6    253.2140 5.663689 0.107865014
  # 7    231.7877 5.462674 0.084789366
  # 8    251.5453 5.828589 0.027363451
  # 9    257.0002 5.516868 0.002774478
  # 10   221.7492 5.292539 0.301262998
  # 11   261.0061 6.128857 0.026122457
  # 12   232.0981 5.593089 0.024510670
  # 13   245.7816 5.638280 0.069414263
  # 14   270.1866 6.021765 0.023402579
  # 15   245.6612 4.863838 0.142257835
  # 16   248.5643 5.386953 0.187815957
  # 17   235.6219 5.737261 0.049328918
  # 18   243.5725 5.302713 0.208180140
  # 19   257.1620 5.883405 0.049343504
  # 20   244.9125 5.539280 0.250891621
  # 21   252.1294 5.898792 0.084715352
  # 22   227.8890 5.309259 0.348106763
  # 23   231.4689 5.222650 0.096292489
  # 24   248.4044 5.348010 1.409227960
  # 25   238.3314 6.329836 0.050602463
  # 26   257.3844 5.986492 0.044588809
  # 27   257.8077 5.366971 0.147622576
  # 28   246.1707 5.335298 0.070301701
  # 29   220.4941 5.125831 0.128028656
  # 30   244.7251 5.611564 0.107879176
  # Mean 244.7916 5.595811 0.146037228
  # 
  244.66 - 22 * 5.61  
  1/sqrt(0.146037228)
  