require(ggdmc); require(data.table); require(ggplot2); require(gridExtra)
require(testthat)
context("check_BuildDMI")

test_that("check_BuildDMI", {
  rm(list = ls())
  setwd("/media/yslin/KIWI/Documents/ggdmc/")
  model <- BuildModel(
    p.map     = list(alpha = "1", beta = "1", epsilon = "1"),
    match.map = NULL,
    factors   = list(X = c("x1")),
    responses = "r1",
    constants = NULL,
    type      = "glm")
  p.vector <- c(alpha = 242.7, beta = 6.185, epsilon = 6.086)
  ntrial <- 1000
  # dat <- simulate(model, nsim = ntrial, ps = p.vector)

  x <- rep(c(8, 15, 22, 29, 36), each = ntrial/5)
  mu <- p.vector[1] + p.vector[2] * x
  y  <- rnorm(ntrial, mu, p.vector[3])
  dat <- data.frame(X = x, R = rep("r1", ntrial), RT = y)
  dplyr::tbl_df(dat)
  dmi <- BuildDMI(dat, model)
  
  res <- ggdmc:::check_BuildDMI(dat, model)

  ## Case 2  
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
  res <- ggdmc:::check_BuildDMI(dat, model)
  res
  
  ## Case3, Set-up PLBA 1f Model -----------
  rm(list = ls())
  wkdir <- "/media/yslin/KIWI/Documents/ppda_paper/"
  setwd(wkdir);
  load("data/holmes2016.rda");
  source("R/utils.R")

  ## Extract switch time using switch trials.
  subjects <- levels(Sw$s)
  ns       <- length(levels(Sw$s))
  model           <- vector("list", ns)
  names(model)    <- subjects
  constants        <- numeric(length(4:18))
  names(constants) <- c(paste0("swt.B", 4:18))
  
  for (i in subjects) {
    Swi <- Sw[Sw$s==i, ]
    constants[1:length(constants)] <- tapply(Swi$ST, Swi$BL, unique)
    new.constants <- c(constants, sd_v = 1)
    model[[i]]   <- ggdmc::BuildModel(
      p.map     = list(A = "1", B = "1", mean_v = "M", mean_w = "M", sd_v = "1",
                       rD = "1", t0 = "1", swt = "BL"),
      match.map = list(M = list(RL = "right", LR = "left")),
      factors   = list(BL = paste0("B", 4:18), S = c("LR", "RL")),
      constants = new.constants,
      responses = c("right", "left"),
      type      = "plba1_gpu")   ## type selects LBA 1f model, using GPU
  }

  ## Mismatch data frame
  p.vector <- c(alpha = 242.7, beta = 6.185, epsilon = 6.086)
  ntrial <- 100
  x <- rep(c(8, 15, 22, 29, 36), each = ntrial/5)
  mu <- p.vector[1] + p.vector[2] * x
  y  <- rnorm(ntrial, mu, p.vector[3])
  dat <- data.frame(X = x, R = rep("r1", ntrial), RT = y)
  dplyr::tbl_df(dat)
  
  res <- ggdmc:::check_BuildDMI(dat, model)
  
  dat$s <- factor(1:nrow(dat))
  res <- ggdmc:::check_BuildDMI(dat, model)
  
  tx2$S <- factor(ifelse(tx2$S == "L", "LR",
                  ifelse(tx2$S == "R",  "RL",
                  ifelse(tx2$S == "RL", "RL",
                  ifelse(tx2$S == "LR", "LR", NA))))) ## NA is a check
  dat <- data.frame(tx2[, .(s, BL, S, R, RT)])
  res <- ggdmc:::check_BuildDMI(dat, model)
  res$issm

})


