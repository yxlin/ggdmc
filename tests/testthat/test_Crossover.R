require(ggdmc); require(testthat); require(ggplot2); require(data.table)
context("Crossover")

test_that("Crossover", {
  rm(list = ls())
  setwd("/media/yslin/KIWI/Documents/ggdmc")
  model <- BuildModel(
    p.map     = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1",
                     st0 = "1"),
    match.map = list(M = list(s1 = 1, s2 = 2)),
    factors   = list(S = c("s1", "s2")),
    constants = c(st0 = 0, sd_v = 1),
    responses = c("r1", "r2"),
    type      = "norm"
    )
  
  ggdmc:::CrossoverDMCHyperchains_blocked
  
})



