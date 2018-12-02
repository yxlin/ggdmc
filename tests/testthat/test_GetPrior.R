rm(list = ls())
setwd("/media/yslin/KIWI/Documents/ggdmc/")
model <- BuildModel(
  p.map      = list(a = "1", b = "1", sd = "1"),
  match.map  = NULL,
  regressors = c(8, 15, 22, 29, 36),
  factors    = list(S = c("x1")),
  responses  = "r1",
  constants  = NULL,
  type       = "glm")
npar <- length(GetPNames(model))
p.prior  <-BuildPrior(
  dists = c("tnorm", "tnorm", "constant"),
  p1    = c(a = 200, b = 10, sd = 0),
  p2    = c(a = 100, b = 5, sd = 100),
  lower = c(NA, NA, 0),
  upper = rep(NA, npar))
plot(p.prior)

## Level 3 ---------
p1.abs <- BuildPrior(
  dists = c(a = "tnorm", b = "tnorm", sd = NA),
  p1    = c(a =  0, b = 0, sd = 0),
  p2    = c(a = 300, b = 10, sd = .01),
  lower = c(0, 0, NA),
  upper = c(NA, NA, NA))
(res1 <- rprior_vec(p1.abs))

p2.abs <- BuildPrior(
  dists = c(a = "tnorm", b = "tnorm", sd = NA),
  p1    = c(a = 0, b = 0, sd = 100),
  p2    = c(a = 100, b = 10, sd = .01),
  lower = c(0, 0, NA),
  upper = c(NA, NA, NA))
(res2 <- rprior_vec(p2.abs))

level3 <- c(unclass(p1.abs), unclass(p2.abs))
names(level3) <- c("p1.a", "p1.b", "p1.sd", "p2.a", "p2.b", "p2.sd")
class(level3) <- "prior"
pp.vector <- c(res1, res2)
round(pp.vector, 2)
plot(level3, ps = pp.vector)

## level 2 (p.prior)
(pop.location <- c(pp.vector[1:2], 0))
(pop.scale <- c(pp.vector[4:5], 100))
names(pop.location) <- c("a", "b", "sd")
names(pop.scale) <- c("a", "b", "sd")

p1 <- c(pop.location[1:2], sd = NA)
p2 <-  c(pop.scale[1:2], sd = NA)
lower <- c(-Inf, 0, NA)
upper <- c(Inf, Inf, NA)
islog <- rep(1, 3)
dists0 <- c(1, 1, NA)
pnames <- GetPNames(model)

# tmp0 <- ggdmc:::RestorePrior(p.prior, p1, p2, lower, upper ,islog, dists0, pnames)
# print(tmp0)
# unclass(tmp0)
# plot(tmp0)

tmp1 <- ggdmc:::GetPrior_test(p.prior, dists0, p1, p2, lower, upper, islog)
print(p.prior)
