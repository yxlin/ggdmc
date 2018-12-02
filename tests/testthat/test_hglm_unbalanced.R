rm(list = ls())
setwd("/media/yslin/KIWI/Documents/ggdmc/")
loadPackages <-c("ggdmc", "ggplot2", "gridExtra",
                 "data.table")
sapply(loadPackages, require, character.only=TRUE)
##              mean     sd    2.5%     25%     50%     75%   97.5%  Rhat n.eff
## alpha0    106.527  3.644  99.390 104.100 106.500 108.900 113.702 1.001  4400
## alpha.c   242.681  2.766 237.100 240.900 242.700 244.500 248.100 1.001  5900
## beta.c      6.189  0.108   5.974   6.118   6.189   6.259   6.404 1.001  6700
## sigma       6.103  0.472   5.273   5.768   6.071   6.397   7.117 1.001 12000

## alpha.tau   0.005  0.001   0.003   0.004   0.005   0.006   0.008 1.001 12000
## beta.tau    4.153  1.565   1.952   3.078   3.879   4.922   8.026 1.001  4400
## deviance  967.271 14.901 940.800 956.800 966.300 976.600 999.700 1.001  4100

# 1/sqrt(.005)
# 1/sqrt(4.153)
## 14.14, .491
##### Simulation study
## Level 3 ---------
p1.abs <- BuildPrior(
  dists = c(a = "tnorm2", b = "tnorm2", tau = NA),
  p1    = c(a =  50, b = 0, sd = NA),
  p2    = c(a = 1e-4, b = 1e-2, sd = NA),
  lower = c(0, 0, NA),
  upper = c(NA, NA, NA))
set.seed(235)
(res1 <- rprior_vec(p1.abs))

p2.abs <- BuildPrior(
  dists = c(a = "tnorm", b = "tnorm", tau = NA),
  p1    = c(a = 0, b = 0, sd = NA),
  p2    = c(a = .01, b = 5, sd = NA),
  lower = c(0, 0, NA),
  upper = c(NA, NA, NA))
set.seed(237)
(res2 <- rprior_vec(p2.abs))

level3 <- c(unclass(p1.abs), unclass(p2.abs))
names(level3) <- c("p1.a", "p1.b", "p1.tau", "p2.a", "p2.b", "p2.tau")
class(level3) <- "prior"
pp.vector <- c(res1, res2)
round(pp.vector, 2)
plot(level3, ps = pp.vector)

## level 2 (p.prior)
(pop.location <- c(pp.vector[1:2], 0))
(pop.scale <- c(pp.vector[4:5], 1))
names(pop.location) <- c("a", "b", "tau")
names(pop.scale) <- c("a", "b", "tau")

pop.prior <- BuildPrior(
  dists = c(a = "tnorm2", b = "tnorm2", tau = "unif"),
  p1    = pop.location,    ## true location
  p2    = pop.scale,       ## true tau / scale
  lower = c(NA, 0, NA),
  upper = c(NA, NA, NA))
print(pop.prior)
plot(pop.prior)
## 1/sqrt(pop.scale)

## Level 1 ----------
## Randomly draw 5 a's (from 1 to ns)
## Use Gaussian GLM to draw trial-level data 
model <- BuildModel(
  p.map      = list(a = "1", b = "1", tau = "1"),
  match.map  = NULL,
  regressors = c(8, 15, 22, 29, 36),
  factors    = list(S = c("x1")),
  responses  = "r1",
  constants  = NULL,
  type       = "glm")
ns <- 30
ntrial <- 100
npar <- length(GetPNames(model))

ps <- GetParameterMatrix(model, ns = ns, prior = pop.prior, seed = 123)
dat <- simulate(model, nsim = ntrial, nsub = ns, ps = ps)
ps <- attr(dat, "parameters")
round(colMeans(ps), 2)
round(matrixStats::colSds(ps), 2)
range(ps[,3])
plot(ps[,3])


p0 <- ggplot(dat, aes( x = X, y = RT)) +
  geom_point() +
  facet_wrap(~s) +
  theme_bw(base_size = 14)
print(p0)

dmi <- BuildDMI(dat, model)

## 1/ c(a = 100, b = 20, tau = 5)^2
## set tnorm estimates sd???
p.prior  <-BuildPrior(
  dists = c("tnorm2", "tnorm2", "tnorm2"),
  p1    = c(a = 200, b = 5, tau = 0),
  p2    = 1 / c(a = 100, b = 20, tau = 5)^2,
  lower = c(NA, NA, 0),
  upper = rep(NA, npar))
plot(p.prior, ps = ps)
plot(p.prior, ps = pop.location)
print(p.prior)

mu.prior <- BuildPrior(
  dists = c(rep("tnorm2", 2), NA),
  p1    = c(pop.location[1:2], tau = NA),
  p2    = c(pop.scale[1:2], tau = NA),
  lower = c(NA, 0, NA),
  upper = rep(NA, npar))
plot(mu.prior, ps = pop.location)

tau.prior <- BuildPrior(
  dists = c("tnorm2", "tnorm2", NA),
  p1    = c(a = 0, b = 0, tau = NA),
  p2    = c(a = .01, b = .01, tau = NA),
  lower = c(0, 0, NA),
  upper = rep(NA, npar))
plot(tau.prior, ps = pop.scale)
pp.prior <- list(mu.prior, tau.prior)


## Sampling -----------
(lstart <- rprior(mu.prior))
(sstart <- rprior(tau.prior))
fit0 <- run(StartNewHypersamples(5e2, dmi, p.prior, pp.prior, 1, 9, .001))
fit <- fit0
thin <- 16
fit <- run(RestartHypersamples(5e2, fit, thin = thin), hpm0 = .05)
rhat <- hgelman(fit, verbose = TRUE)

## Analysis -----------
p0 <- plot(fit, hyper = TRUE)
p0 <- plot(fit, hyper = TRUE, start = 101)
p0 <- plot(fit, hyper = TRUE, pll = FALSE, den = TRUE, start = 101)

p0 <- plot(fit, pll = FALSE, den = TRUE, start = 1)
p1 <- plot(fit, start = 101)
p1 <- plot(fit)

1/sqrt(pop.scale)

est1 <- summary(fit, hyper = TRUE, recover = TRUE, ps = pop.location[1:2], type = 1, verbose = TRUE)
est2 <- summary(fit, hyper = TRUE, recover = TRUE, ps = pop.scale[1:2], type = 2, verbose = TRUE)
est3 <- summary(fit, recover = TRUE, ps = ps, verbose = FALSE)
est4 <- summary(fit, verbose = FALSE)

range(est4[1:ns,3])
plot(est4[1:ns,3])
hist(est4[1:ns,1], breaks = "fd")
hist(est4[1:ns,2], breaks = "fd")
hist(est4[1:ns,3], breaks = "fd")
