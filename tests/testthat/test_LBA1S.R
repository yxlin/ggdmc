require(ggdmc); require(testthat); require(ggplot2); require(data.table)
context("LBA 1e4")

test_that("LBA 1e4", {
  rm(list = ls())
  setwd("/media/yslin/KIWI/Documents/ggdmc")
  model <- BuildModel(
    p.map     = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1",
                     st0 = "1"),
    match.map = list(M = list(s1 = 1, s2 = 2)),
    factors   = list(S = c("s1", "s2")),
    constants = c(st0 = 0, sd_v = 1),
    responses = c("r1", "r2"),
    type      = "norm")


  p.map     = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1",
                   st0 = "1")
  match.map = list(M = list(s1 = 1, s2 = 2))
  factors   = list(S = c("s1", "s2"))
  constants = c(st0 = 0, sd_v = 1)
  responses = c("r1", "r2")
  type      = "norm"
  posdrift   = TRUE
  verbose    = TRUE
  
  p.vector <- c(A = .75, B = 1.25, t0 = .15, mean_v.true = 2.5, mean_v.false = 1.5)
  ntrial <- 1e4
  dat <- simulate(model, nsim = ntrial, ps = p.vector)
  dmi <- BuildDMI(dat, model)
  # n1.order   <- attr(model, "n1.order")
  # row.names(n1.order)
  # dimnames(model)[[1]]
  
  
  p.prior <- BuildPrior(
    dists = c("tnorm", "tnorm", "beta", "tnorm", "tnorm"),
    p1    = c(A = 1, B = 1, t0 = 1, mean_v.true = 1, mean_v.false = 1),
    p2    = c(1,  1,  1, 1, 1),
    lower = c(rep(0, 3),  rep(NA, 2)),
    upper = c(rep(NA, 2), 1, rep(NA, 2)))
  plot(p.prior, ps = p.vector)

  # res <- ggdmc:::check_BuildDMI(dat, model)
  # subject_models <- res$issm
  # modeli <- res$model
  # fnams <- names(attr(modeli, "factors"))
  # res <- ggdmc:::check_BuildDMI(dat, model)
  # subject_models <- res$issm
  # modeli <- res$model
  # fnams <- names(attr(modeli, "factors"))
  # cells <- apply(dat[, c(fnams, "R")], 1, paste, collapse = ".")
  # str(cells)
  
  
  # cell.index <- vector("list", dim(model)[1])
  # dim1 <- dimnames(model)[[1]]
  # names(cell.index) <- dim1
  # str(cell.index)
  
  # for ( j in names(cell.index) ) cell.index[[j]] <- cells %in% j
  # str(cell.index)

  ## Sampling ---------
  setwd("/media/yslin/KIWI/Documents/ggdmc_lesson")
  path <- c("data/Lesson3/ggdmc_3_2_LBA1S_1e4_tmp.rda")
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
  p2 <- plot(fit0, start = 201)
  
  # gridExtra::grid.arrange(p0, p1, p2, ncol = 3)
  # d <- data.table(fit$data)
  # d[, .N, .(S)]
  
  ## Analysis -----------
  est <- summary(fit, start = 201, recovery = TRUE, ps = p.vector, verbose = TRUE)
  #                    A    B mean_v.false mean_v.true   t0
  # True            0.75 1.25         1.50        2.50 0.15
  # 2.5% Estimate   0.62 1.14         1.45        2.46 0.13
  # 50% Estimate    0.74 1.25         1.51        2.50 0.15
  # 97.5% Estimate  0.83 1.38         1.57        2.55 0.17
  # Median-True    -0.01 0.00         0.01        0.00 0.00
  ##                   A    B   t0 mean_v.true mean_v.false
  ## True           0.75 0.25 0.20        2.50         1.50
  ## 2.5% Estimate  0.60 0.17 0.17        2.24         1.12
  ## 50% Estimate   0.75 0.25 0.20        2.59         1.48
  ## 97.5% Estimate 0.90 0.35 0.22        2.95         1.84
  ## Median-True    0.00 0.00 0.00        0.09        -0.02
  
  ##                    A     B mean_v.false mean_v.true   t0
  ## True            0.75  1.25         1.50        2.50 0.15
  ## 2.5% Estimate   0.56  1.08         1.37        2.40 0.12
  ## 50% Estimate    0.72  1.23         1.46        2.47 0.15
  ## 97.5% Estimate  0.84  1.41         1.55        2.53 0.18
  ## Median-True    -0.03 -0.02        -0.04       -0.03 0.00
  
  #                    A    B mean_v.false mean_v.true    t0
  # True            0.75 1.25         1.50        2.50  0.15
  # 2.5% Estimate   0.51 1.21         1.50        2.50  0.08
  # 50% Estimate    0.71 1.39         1.58        2.57  0.12
  # 97.5% Estimate  0.86 1.63         1.68        2.64  0.15
  # Median-True    -0.04 0.14         0.08        0.07 -0.03
  #                    A    B mean_v.false mean_v.true    t0
  # True            0.75 1.25         1.50        2.50  0.15
  # 2.5% Estimate   0.37 1.21         1.44        2.44  0.09
  # 50% Estimate    0.61 1.40         1.54        2.51  0.12
  # 97.5% Estimate  0.77 1.65         1.63        2.58  0.15
  # Median-True    -0.14 0.15         0.04        0.01 -0.03
  
  # A    B mean_v.false mean_v.true    t0
  # True            0.75 1.25         1.50        2.50  0.15
  # 2.5% Estimate   0.41 1.23         1.51        2.50  0.09
  # 50% Estimate    0.62 1.43         1.60        2.57  0.12
  # 97.5% Estimate  0.79 1.65         1.69        2.64  0.15
  # Median-True    -0.13 0.18         0.10        0.07 -0.03
  
  #                   A     B mean_v.false mean_v.true   t0
  # True           0.75  1.25         1.50        2.50 0.15
  # 2.5% Estimate  0.70  1.01         1.32        2.38 0.15
  # 50% Estimate   0.82  1.11         1.40        2.44 0.17
  # 97.5% Estimate 0.92  1.27         1.48        2.50 0.19
  # Median-True    0.07 -0.14        -0.10       -0.06 0.02
  #                   A     B mean_v.false mean_v.true   t0
  # True           0.75  1.25         1.50        2.50 0.15
  # 2.5% Estimate  0.67  1.04         1.40        2.42 0.14
  # 50% Estimate   0.80  1.17         1.47        2.48 0.17
  # 97.5% Estimate 0.92  1.33         1.54        2.54 0.19
  # Median-True    0.05 -0.08        -0.03       -0.02 0.02
  
  # A    B mean_v.false mean_v.true    t0
  # True            0.75 1.25         1.50        2.50  0.15
  # 2.5% Estimate   0.40 1.19         1.42        2.43  0.10
  # 50% Estimate    0.61 1.37         1.51        2.50  0.13
  # 97.5% Estimate  0.76 1.57         1.59        2.56  0.16
  # Median-True    -0.14 0.12         0.01        0.00 -0.02
  
})



