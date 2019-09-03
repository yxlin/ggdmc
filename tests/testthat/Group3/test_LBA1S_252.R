cat("\n-------------------- Testing LBA 1 Subject --------------------")

rm(list = ls())
model0 <- ggdmc:::BuildModel(
  p.map     = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1",
                   st0 = "1"),
  match.map = list(M = list(s1 = 1, s2 = 2)),
  factors   = list(S = c("s1", "s2")),
  constants = c(st0 = 0, sd_v = 1),
  responses = c("r1", "r2"),
  type      = "norm")
p.vector <- c(A = .75, B = 1.25, t0 = .15, mean_v.true = 2.5, mean_v.false = 1.5)
ntrial <- 100
dat0 <- ggdmc:::simulate.model(model0, nsim = ntrial, ps = p.vector)
dmi0 <- ggdmc:::BuildDMI(dat0, model0)

p.prior0 <- ggdmc:::BuildPrior(
  dists = c("tnorm", "tnorm", "beta", "tnorm", "tnorm"),
  p1    = c(A = 1, B = 1, t0 = 1, mean_v.true = 1, mean_v.false = 1),
  p2    = c(1,  1,  1, 1, 1),
  lower = c(rep(0, 3),  rep(NA, 2)),
  upper = c(rep(NA, 2), 1, rep(NA, 2)))


fit0 <- ggdmc:::StartNewsamples(dmi0, p.prior0, block=FALSE, pm0=0, pm1=0, nmc=2e2)
fit <- ggdmc:::run(fit0, block=FALSE, thin = 1)

fit <- ggdmc:::run(fit, block=FALSE, thin = 4)
fit <- ggdmc:::run(fit, block=FALSE, thin = 16)

fit1 <- ggdmc:::run(fit, block=FALSE, thin = 1, pm1=.05)
fit2 <- ggdmc:::run(fit1, block=FALSE, thin = 1)
hat <- ggdmc:::gelman(fit); hat
hat <- ggdmc:::gelman(fit1); hat
hat <- ggdmc:::gelman(fit2); hat

est <- ggdmc:::summary.model(fit, recovery = TRUE, ps = p.vector, verbose = TRUE)
est <- ggdmc:::summary.model(fit1, recovery = TRUE, ps = p.vector, verbose = TRUE)

p0 <- ggdmc:::plot.model(fit0, start = 21)
p1 <- ggdmc:::plot.model(fit)
p2 <- ggdmc:::plot.model(fit, pll = F, den = T)




model <- ggdmc252:::BuildModel(
    p.map     = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1",
                     st0 = "1"),
    match.map = list(M = list(s1 = 1, s2 = 2)),
    factors   = list(S = c("s1", "s2")),
    constants = c(st0 = 0, sd_v = 1),
    responses = c("r1", "r2"),
    type      = "norm")
dat <- ggdmc252:::simulate.model(model, nsim = ntrial, ps = p.vector)




dmi <- ggdmc252:::BuildDMI(dat, model)
class(dmi) <- c("dmi", "model", "data.frame")
dmi1 <- ggdmc:::BuildDMI(dat, model1)

attr(model, "is.r1")
attr(model1, "is.r1")

names(attributes(dmi))
attr(dmi, "cell.empty")
attr(dmi1, "cell.empty")
dim(model)

p.prior <- ggdmc252::BuildPrior(
    dists = c("tnorm", "tnorm", "beta", "tnorm", "tnorm"),
    p1    = c(A = 1, B = 1, t0 = 1, mean_v.true = 1, mean_v.false = 1),
    p2    = c(1,  1,  1, 1, 1),
    lower = c(rep(0, 3),  rep(NA, 2)),
    upper = c(rep(NA, 2), 1, rep(NA, 2)))


## Sampling ---------
fit0 <- ggdmc252:::StartNewsamples(200, dmi, p.prior)
fit  <- ggdmc252::run(fit0)
fit  <- ggdmc252::run( ggdmc252::RestartSamples(500, fit, thin = 1) )
hat  <- ggdmc252::gelman(fit); hat

est  <- ggdmc252:::summary.model(fit, recovery = TRUE, ps = p.vector, verbose = TRUE)

# p0 <- ggdmc252:::plot.model(fit)
# p1 <- ggdmc252:::plot.model(fit, pll = F, den = T)
npar <- length(ggdmc::GetPNames(model))



pdf(file = "LBA1S.pdf")

p2 <- ggdmc:::plot.model(fit0, start = 51)
p2 <- ggdmc:::plot.model(fit)
p2 <- ggdmc:::plot.model(fit, pll = F, den = T)
dev.off()

## Analysis -----------
est <- ggdmc:::summary.model(fit, recovery = TRUE, ps = p.vector, verbose = TRUE)

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

  #                    A    B mean_v.false mean_v.true    t0
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

  #                    A    B mean_v.false mean_v.true    t0
  # True            0.75 1.25         1.50        2.50  0.15
  # 2.5% Estimate   0.40 1.19         1.42        2.43  0.10
  # 50% Estimate    0.61 1.37         1.51        2.50  0.13
  # 97.5% Estimate  0.76 1.57         1.59        2.56  0.16
  # Median-True    -0.14 0.12         0.01        0.00 -0.02




