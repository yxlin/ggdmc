require(ggdmc); require(data.table); require(ggplot2); require(gridExtra)
require(testthat)
context("hlogit")

test_that("hlogit", {
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
  pnames <- GetPNames(model)
  npar <- length(GetPNames(model))

  pop.location <- c(a0 = -.55, a1 = .08, a2 = -.81, a3 = 1.35, b = .27)
  pop.scale    <- c(a0 = .19, a1 = .30, a2 = .41, a3 = .26, b = .15)
  ntrial <- 30
  nsub <- 12
  pop.prior  <-BuildPrior(
    dists = rep("tnorm", npar),
    p1    = pop.location,
    p2    = pop.scale,
    lower = c(NA, NA, NA, NA, 0),
    upper = rep(NA, npar))
  plot(pop.prior, ps = pop.location)
  dat <- simulate(model, nsub = nsub, nsim = ntrial, prior = pop.prior)
  dmi <- BuildDMI(dat, model)
  head(dat)
  tail(dat)
  dplyr::tbl_df(dat)
  ps <- attr(dat, "parameters")
  round(ps, 2)

  pstart  <- BuildPrior(
    dists = c(rep("tnorm", npar)),
    p1    = pop.location,
    p2    = pop.scale * 10,
    lower = rep(NA, npar),
    upper = rep(NA, npar))
  lstart <- BuildPrior(
    dists = c(rep("tnorm", npar-1), NA),
    p1    = pop.location,
    p2    = pop.scale * 10,
    lower = rep(NA, npar),
    upper = rep(NA, npar))
  sstart <- BuildPrior(
    dists = c(rep(NA, npar-1), "gamma"),
    p1    = c(a0 = NA,  a1 = NA, a2 = NA, a3 = NA, b = .001),
    p2    = c(NA, NA, NA, NA, 1000),
    lower = rep(NA, npar),
    upper = rep(NA, npar))
  start <- list(pstart, lstart, sstart)
  plot(pstart)
  plot(lstart)
  plot(sstart)
  rprior(pstart)
  rprior(lstart)
  rprior(sstart)

  pprior  <- BuildPrior(
    dists = c(rep("constant", npar-1), "tnorm"),
    p1    = c(a0 = NA,  a1 = NA, a2 = NA, a3 = NA, b = 0),
    p2    = c(a0 = 0,  a1 = 0, a2 = 0, a3 = 0, b = NA),
    lower = c(rep(NA, npar-1), 0),
    upper = rep(NA, npar))
  lprior  <- BuildPrior(
    dists = c(rep("tnorm2", npar-1), NA),
    p1    = c(a0 = 0, a1 = 0, a2 = 0, a3 = 0, b = NA),
    p2    = c(a0 = 1e-6, a1 = 1e-6, a2 = 1e-6, a3 = 1e-6, b = NA),
    lower = rep(NA, npar),
    upper = rep(NA, npar))
  sprior <- BuildPrior(
    dists = c(rep(NA, npar-1), "gamma"),
    p1    = c(a0 = NA, a1 = NA, a2 = NA, a3 = NA, b = .001),
    p2    = c(a0 = NA, a1 = NA, a2 = NA, a3 = NA, b = 1000),
    lower = rep(0, npar),
    upper = rep(NA, npar))
  prior <- list(pprior, lprior, sprior)
  print(pstart)
  print(pprior)

  attr(dmi[[1]], "colnames")
  dimnames(dmi[[1]])
  names(attributes(dmi[[1]]))
  names(dmi[[1]])

  # DT <- data.table(d)
  # DT[, .N, .(s)]
  (nsub <- length(dmi))
  # tmp0 <- ggdmc:::init_logit_test(100, dmi[[1]], pstart, pprior, .001, 1, 15)
  args(Starthlogit)
  tmp0 <- Starthlogit(500, dmi, start, prior, 1, 15, .001)
  tmp1 <- tmp0[[1]]
  str(tmp1$theta)
  tmp1$theta[,,1]

  tmp1 <- print(lprior)
  tmp2 <- print(sprior)

  nlpar <- sum(!is.na(tmp1$dist))
  nspar <- sum(!is.na(tmp2$dist))

  nmc <- 500
  nchain <- 15
  location <- scale <-  array(dim = c(nchain, npar, nmc))
  hlp <- hll <- matrix(-Inf, nmc, nchain)
  theta <- GetTheta0(tmp0)  ## nsub x npar x nchain
  theta[,,1]
  str(theta)

  i <- 1
  print(lprior)
  print(sprior)
  for(i in 1:nchain) {
    (lvec <- rprior(lstart))
    (svec <- rprior(sstart))
    llik <- sumlogpriorNV(lvec, lprior[names(lvec)])
    slik <- sumlogpriorNV(svec, sprior[names(svec)])
    hlp[1,i] <- llik+slik

    location[i,,1] <- lvec
    scale[i,,1]    <- svec

    for(j in 1:nsub) {
      res <- res + ggdmc:::dprior_(thetak[j,], tmpp$dist, lvec, svec, tmpp$lower, tmpp$upper, tmpp$lg)
      res <- ggdmc:::sumlogprior(thetak[j,], tmpp$dist, lvec, svec, tmpp$lower,
                                 tmpp$upper, tmpp$lg)
    }

    hll.row(0).col(i) = sumloghlike(theta.slice(i).cols(pidx_fin),
            dist_pp(pidx_fin), lvec(pidx_fin), svec(pidx_fin),
            lower(pidx_fin), upper(pidx_fin), lg(pidx_fin));

  }


  fit0 <- run(tmp0)
  fit <- run(RestartSamples(5e2, fit0, thin = 8))


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


