require(ggdmc); require(testthat); require(ggplot2); require(data.table)
context("Location prior")

test_that("Location prio", {
  rm(list = ls())
  setwd("/media/yslin/KIWI/Documents/ggdmc")
  npar <- 3
  ## Test rprior_vec
  location.prior <- BuildPrior(
    p1    = c(a = 0, b = 0, s = 0),
    p2    = c(100, .001, .001),
    lower = c(NA, NA, NA),
    upper = rep(NA, npar),
    dists = c("tnorm", "constant", "constant"))
  
  scale.prior <- BuildPrior(
    p1    = c(a = 0, b = 100, s = 100),
    p2    = c(a = 100, b = .001, s = .001),
    lower = c(NA, NA, NA),
    upper = rep(NA, npar),
    dists = c("unif", "constant", "constant"))
  
  p.prior  <-BuildPrior(
    dists = c(a = "tnorm", b = "tnorm", s = "unif"),
    p1    = c(a = 0, b = 10.5, s = 0), ## unimportant in hierarchical
    p2    = c(a = 1, b = 5, s = 5),    ## unimportant in hierarchical
    lower = c(NA, NA, 0),
    upper = rep(NA, npar))
  str(p.prior)
  
  plot(p.prior)
  plot(location.prior)
  plot(scale.prior)
  print(scale.prior)
  print(location.prior)
  pnames <- names(location.prior)
  pp.prior <- list(location.prior, scale.prior)

  pdists <- character(npar)
  pp1 <- pp2 <- plower <- pupper <- numeric(npar)
  islogp <- logical(npar)
  ldists <- character(npar)
  lp1 <- lp2 <- llower <- lupper <- numeric(npar)
  islogl <- logical(npar)
  sdists <- character(npar)
  sp1 <- sp2 <- slower <- supper <- numeric(npar)
  islogs <- logical(npar)
  
  for(i in 1:npar) {
    pdists[i] <- attr(p.prior[[i]], "dist")
    pp1[i] <- p.prior[[i]][[1]]
    pp2[i] <- p.prior[[i]][[2]]
    plower[i] <- p.prior[[i]][[3]]
    pupper[i] <- p.prior[[i]][[4]]
    islogl[i] <- p.prior[[i]][[5]]
    
    ldists[i] <- attr(location.prior[[i]], "dist")
    lp1[i] <- location.prior[[i]][[1]]
    lp2[i] <- location.prior[[i]][[2]]
    llower[i] <- location.prior[[i]][[3]]
    lupper[i] <- location.prior[[i]][[4]]
    islogl[i] <- location.prior[[i]][[5]]
    
    sdists[i] <- attr(scale.prior[[i]], "dist")
    sp1[i] <- scale.prior[[i]][[1]]
    sp2[i] <- scale.prior[[i]][[2]]
    slower[i] <- scale.prior[[i]][[3]]
    supper[i] <- scale.prior[[i]][[4]]
    islogs[i] <- scale.prior[[i]][[5]]
  }

  ggdmc:::rprior_vec(ldists, lp1, lp2, llower, lupper)
  ggdmc:::rprior_vec(sdists, sp1, sp2, slower, supper)
  ggdmc:::rprior_vec(ldists, lp1, lp2, llower, lupper)
  
  tmp0 <- ggdmc:::RestorePrior(lp1, lp2, llower, lupper, islogl, ldists, pnames)
  str(tmp0)
  class(tmp0) <- "prior"
  plot(tmp0) ## does not work because name conflict mean sd vs p1 p2
  print(tmp0)

  # ggdmc:::dprior_(pvec, ldists, lp1, lp2, llower, lupper, islogl)
  # dnorm(pvec[1], lp1[1], lp2[1], log = TRUE)
  
})



