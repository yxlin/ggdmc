require(ggdmc); require(testthat)
context("rlba_norm")

test_that("make_v", {
  rm(list = ls())
  mean_v <- matrix(c(2.4, 2.2)); mean_v
  drifts <- ggdmc:::make_v(10, mean_v, 1, TRUE); drifts
})

test_that("Make R", {
  mean_v <- matrix(c(2.4, 2.2), nrow = 2); mean_v
  drifts <- ggdmc:::make_v(10, mean_v, 1, TRUE)
  A  <- 1.2
  b  <- 2.7
  t0 <- .2
  st0 <- 0

  is.vector(mean_v)
  is.matrix(mean_v)

  dat0 <- ggdmc:::make_r(drifts, A, b, t0, st0, TRUE); dat0
  dat0 <- ggdmc:::make_r(drifts, A, b, t0, st0, FALSE); dat0
  dat1 <- data.frame(rca(10, A, b, t0, mean_v, 1, 0, 0.34)); dat1
  colnames(dat1) <- c("RT", "R")
  dat1$R <- factor(dat1$R)
  dplyr::tbl_df(dat1)

  x <- seq(0, 3, .001)
  dat2 <- n1PDF_cnorm(x, A, b, t0, mean_v, 1, 0, 0.34, 1048576, .01);
  plot(dat2[,1])
  ## ggdmc:::density_cnorm_pda
  ## ggdmc:::n1PDF_gpu

  # dat1 <- cam::MakeRs(drifts, A, b, t0, st0, TRUE); dat1
  # dat1 <- cam::MakeRs(drifts, A, b, t0, st0, FALSE); dat1
  # dat2 <- cam::makeR(drifts, A, b, t0, st0, TRUE); dat2
  # dat2 <- cam::makeR(drifts, A, b, t0, st0, FALSE); dat2

  ## rlba_norm ----
  n <- 10
  posdrift <- TRUE
  sd_v <- 1
  dat3 <- ggdmc:::rlba_norm(n, A, b, mean_v, sd_v, t0, st0, posdrift, TRUE); dat3
  dat3 <- ggdmc:::rlba_norm(n, A, b, mean_v, sd_v, t0, st0, posdrift, FALSE); dat3

})

