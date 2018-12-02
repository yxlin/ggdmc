rm(list = ls())
setwd("/media/yslin/KIWI/Documents/BUGS_Examples/vol1/Rats/")
loadPackages <-c("R2OpenBUGS", "coda", "ggdmc", "ggplot2", "gridExtra",
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
  p2    = pop.scale,       ## true scale
  lower = c(NA, 0, NA),
  upper = c(NA, NA, NA))
print(pop.prior)
plot(pop.prior)
1/sqrt(pop.scale)

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

1/sqrt(ps[,3])
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

1/ c(a = 100, b = 20, tau = 5)^2
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

# p2.abs <- BuildPrior(
#   dists = c(a = "tnorm", b = "tnorm", tau = NA),
#   p1    = c(a = 0, b = 0, sd = NA),
#   p2    = c(a = .01, b = 5, sd = NA),
#   lower = c(0, 0, NA),
#   upper = c(NA, NA, NA))
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

tmp0 <- ggdmc:::init_newhier_hstart(100, dmi, p.prior, pp.prior, .001, 1, 9,
                                    lstart, sstart)
tmp0 <- ggdmc:::init_newhier_glm(100, dmi, p.prior, pp.prior, .001, 1, 9)

fit0 <- run(StartNewHypersamples(5e2, dmi, p.prior, pp.prior, 1, 9, .001))
fit <- fit0
thin <- 1
fit <- run(RestartHypersamples(5e2, fit, thin = thin))
rhat <- hgelman(fit, verbose = TRUE)

## Analysis -----------
p0 <- plot(fit, hyper = TRUE)
p0 <- plot(fit, hyper = TRUE, start = 101)
p0 <- plot(fit, hyper = TRUE, pll = FALSE, den = TRUE, start = 1)

p0 <- plot(fit, pll = FALSE, den = TRUE, start = 1)
p1 <- plot(fit, start = 101)
p1 <- plot(fit)

1/sqrt(pop.scale)

est1 <- summary(fit, hyper = TRUE, recover = TRUE, ps = pop.location[1:2], type = 1, verbose = TRUE)
est2 <- summary(fit, hyper = TRUE, recover = TRUE, ps = pop.scale[1:2], type = 2, verbose = TRUE)

est3 <- summary(fit, recover = TRUE, ps = ps, verbose = FALSE)
est4 <- summary(fit, verbose = FALSE)

range(est4[1:50,3])
plot(est4[1:50,3])
hist(est4[1:50,1], breaks = "fd")
hist(est4[1:50,2], breaks = "fd")
hist(est4[1:50,3], breaks = "fd")
# hyper <- attr(fit, "hyper")
# hyper$phi[[1]][,,1]
# hyper$phi[[2]][,,1]
# fit3[[1]]$theta[,,1]
# str(hyper)
# 
# i <- 1
# hyper$h_log_likelihoods[i,]   ## nmc x nchain
# hyper$h_summed_log_prior[i,]  ## nmc x nchain
# 
# fit[[1]]$log_likelihoods[i,]  ## nmc x nchain
# fit[[2]]$log_likelihoods[i,]  ## nmc x nchain
# 
# fit[[1]]$summed_log_prior[i,] ## nmc x nchain







# fit0 <- StartNewsamples(5e2, dmi[[1]], p.prior)
fit0 <- run(StartNewsamples(5e2, dmi[[1]], p.prior))
fit0 <- run(StartManynewsamples(5e2, dmi, p.prior, 1, 9, .001))


tmp1 <- ggdmc:::rprior_vec_(ldists, lp1, lp2, llower, lupper);
tmp2 <- ggdmc:::rprior_vec_(sdists, sp1, sp2, slower, supper);




tmp1 <- ggdmc:::init_newnonhier_glm(10, dmi, tmp0, .001, 1, 9)
tmp2 <- ggdmc:::init_new_glm(100, tmp0, dmi[[1]], .001, 1, 9)


print(tmp0)
plot(tmp0)
rprior(tmp0, 1)


plot_prior(1, tmp0)
plot_prior(2, tmp0)
plot_prior(3, tmp0)
dprior(pop.location, p.prior)
ggdmc:::rprior_vec_(c(1, 1, NA), c(pop.location[1:2], sd = NA),
                    c(pop.scale[1:2], sd = NA), c(NA, 0, NA), c(NA, NA, NA))
print(mu.prior)




lstart <- rprior_vec(mu.prior)
sstart <- rprior_vec(sigma.prior)

ldists <- c(1, 1, NA)
lp1 <- c(pop.location[1:2], sd = 0)
lp2 <- c(pop.scale[1:2], sd = 100)
llower <- c(-Inf, 0, NA)
lupper <- c(Inf, Inf, NA)
llog <- c(1, 1, 1)
print(mu.prior)

sdists = c(1, 1, NA)
sp1    = c(a = 100, b = 5, sd = NA)
sp2    = c(a = 50,  b = 5, sd = NA)
slower = c(0, 0, NA)
supper = c(Inf, Inf, NA)
slog <- c(1, 1, 1)
print(sigma.prior)


tmp0 <- ggdmc:::sumloghprior(lstart, sstart, ldists, sdists, lp1, sp1,
                             lp2, sp2, llower, slower, lupper, supper, llog, slog);
tmp1 <- ggdmc:::sumlogprior(lstart, ldists, lp1, lp2, llower, lupper, llog) 
tmp2 <- ggdmc:::sumlogprior(sstart, sdists, sp1, sp2, slower, supper, slog)

tmp3 <- ggdmc:::dprior_(lstart, ldists, lp1, lp2, llower, lupper, llog)
sum(tmp3)

theta0 <- ggdmc:::GetTheta0(fit) ## // thetas: nsub x npar x nchain

hll.row(0).col(i) = sumloghlike(thetas.slice(i), pdists, lstart, sstart,
                                plower, pupper, plog);

pdists <- c(1, 1, 5)
llower <- c(NA, NA, NA)
lupper <- rep(NA, 3)


tmp4 <- ggdmc:::sumloghlike(theta0[,,1], pdists, lstart, sstart, llower, lupper, slog)



setwd("/media/yslin/KIWI/Documents/ggdmc_lesson/")
path <- c("data/Lesson4/ggdmc_4.0_hglm_tmp.rda")


fit0 <- StartNewHypersamples(5e2, dmi, p.prior, pp.prior))

fit0 <- run(StartNewHypersamples(5e2, dmi, p.prior, pp.prior))
fit  <- fit0
thin <- 4
repeat {
  fit <- run(RestartHypersamples(5e2, fit, thin = thin))
  save(fit0, fit, file = path[1])
  rhat <- hgelman(fit, verbose = TRUE)
  if (all(rhat$mpsrf < 1.1)) break
  thin <- thin * 2
}
cat("Done ", path[1], "\n")
setwd("/media/yslin/KIWI/Documents/ggdmc/")



setwd("/media/yslin/KIWI/Documents/BUGS_Examples/vol1/Rats/")
tmp <- dget("data/dataBUGS.R")
d <- data.frame(matrix(as.vector(tmp$Y), nrow = 30, byrow = TRUE))
ns <- nrow(d)
na <- ncol(d)
names(d) <- as.character(c(5, 15, 22, 29, 36))
d$s <- factor(1:ns)
long <- melt(d, id.vars = c("s"), variable.name = "x",
             value.name = "RT")
dplyr::tbl_df(long)
dat <- long
dat$X <- as.numeric(as.character(long$x)) - tmp$xbar

dat$S <- "x1"
dat$R <- "r1"
head(dat)
dmi <- BuildDMI(dat, model)

sd(dat$RT)
npar <- length(GetPNames(model))

# pop.mean <- c(a = 242.7, b = 6.185, sd = 10)
# pop.scale <- c(a = 100, b = 5, sd = 50)

print(p.prior)
print(pp.prior[[1]])



res <- 0
for(i in 1:ns) {
  res <- res + DIC(fit[[1]])  
}
res





## Load data ---------------------------------------------------------
## system("cat blockerDataBugs.R")

fit0 <- run(StartNewHypersamples(5e2, dmi, p.prior, pp.prior))

fit  <- fit0
thin <- 2
p0 <- plot(fit, hyper = TRUE)
p0 <- plot(fit, hyper = TRUE, start = 101)


## Set-up BUGS model -------------------------------------------------
model <- function() {
  ## Likelihood
  for( i in 1 : N ) {
    for( j in 1 : T ) {
      Y[i, j] ~ dnorm(mu[i, j], tau.c)
      mu[i, j] <- alpha[i] + beta[i] * (x[j] - xbar)
      
      ## Other monitoring items
      ##             culmative.Y[i, j] <- cumulative(Y[i, j], Y[i, j])
      ##             post.pv.Y[i, j] <- post.p.value(Y[i, j])
      ##             prior.pv.Y[i, j] <- prior.p.value(Y[i, j])
      ##             replicate.post.Y[i, j] <- replicate.post(Y[i, j])
      ##             pv.post.Y[i, j] <- step(Y[i, j] - replicate.post.Y[i, j])
      ##             replicate.prior.Y[i, j] <- replicate.prior(Y[i, j])
      ##             pv.prior.Y[i, j] <- step(Y[i, j] - replicate.prior.Y[i, j])
      
    }
    ## Prior
    alpha[i] ~ dnorm(alpha.c, alpha.tau)
    beta[i] ~ dnorm(beta.c, beta.tau)
  }
  
  ## Hyper-prior
  tau.c ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau.c)
  alpha.c ~ dnorm(0.0, 1.0E-6)
  alpha.tau ~ dgamma(0.001, 0.001)
  beta.c ~ dnorm(0.0, 1.0E-6)       ## growth rates
  beta.tau ~ dgamma(0.001, 0.001)
  alpha0 <- alpha.c - xbar * beta.c ## initial weight
}

## tempdir()
wkdir <- "/media/yslin/KIWI/Documents/BUGS_Examples/vol1/Rats/models"
model.file <- file.path(wkdir, "model.txt")
write.model(model, model.file)

## Load data ---------------------------------------------------------
## system("cat blockerDataBugs.R")
dat <- dget("data/dataBUGS.R")
ini <- dget("data/inits.R")

N <- dat$N
T <- dat$T
Y <- matrix(as.vector(dat$Y), nrow = 30, byrow = TRUE)
xbar <- dat$xbar
x <- dat$x
data <- list("N", "T", "Y", "xbar", "x")
p.vector <- c("alpha0",  "alpha.c", "beta.c", "sigma", "alpha.tau", "beta.tau")

inits <- function() {
  list(alpha.c = ini$alpha.c,
       beta.c = ini$beta.c,
       alpha.tau = ini$alpha.tau,
       beta.tau = ini$beta.tau,
       tau.c = ini$tau.c)
}

fit <- R2OpenBUGS::bugs(data, inits, p.vector, n.iter = 10000,
                        working.directory = wkdir,
                        model.file, n.chains = 3, n.burnin = 5000)

print(fit, digits = 3)
## Inference for Bugs model at "/media/yslin/KIWI/Documents/BUGS_Examples/vol1/Rats/models/model.txt", 
## Current: 3 chains, each with 10000 iterations (first 5000 discarded)
## Cumulative: n.sims = 15000 iterations saved
##              mean     sd    2.5%     25%     50%     75%   97.5%  Rhat n.eff
## alpha0    106.527  3.644  99.390 104.100 106.500 108.900 113.702 1.001  4400
## alpha.c   242.681  2.766 237.100 240.900 242.700 244.500 248.100 1.001  5900
## beta.c      6.189  0.108   5.974   6.118   6.189   6.259   6.404 1.001  6700
## sigma       6.103  0.472   5.273   5.768   6.071   6.397   7.117 1.001 12000
## alpha.tau   0.005  0.001   0.003   0.004   0.005   0.006   0.008 1.001 12000
## beta.tau    4.153  1.565   1.952   3.078   3.879   4.922   8.026 1.001  4400
## deviance  967.271 14.901 940.800 956.800 966.300 976.600 999.700 1.001  4100
## 
## For each parameter, n.eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
## 
## DIC info (using the rule, pD = Dbar-Dhat)
## pD = 54.300 and DIC = 1022.000
## DIC is an estimate of expected predictive error (lower deviance is better).
