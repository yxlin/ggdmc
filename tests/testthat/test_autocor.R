require(ggdmc); require(ggplot2); library(gridExtra)
setwd("/media/yslin/KIWI/Documents/ggdmc")
rm(list = ls())

ac2 <- function (x, nLags) {
  X <- matrix(NA, ncol = nLags, nrow = length(x))
  X[, 1] <- x ## 1st column is the posterior log-likelihood
  for (i in 2:nLags) { ## 2nd to the nLags-th column push the value of posterior log-likelihood
    ## one 1, 2, etc step back and let the tail gone (i.e., not used)
    X[, i] <- c(rep(NA, i - 1), x[1:(length(x) - i + 1)])
  }

  # X <- list(Lag = 1:nLags, Autocorrelation = cor(X, use = "pairwise.complete.obs")[, 1])
  X <- data.frame(Lag = 1:nLags, Autocorrelation = cor(X, use = "pairwise.complete.obs")[, 1])
  return(X)
}

model <- BuildModel(
  p.map     = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1",
                   st0 = "1"),
  match.map = list(M = list(s1 = 1, s2 = 2)),
  factors   = list(S = c("s1", "s2")),
  constants = c(st0 = 0, sd_v = 1),
  responses = c("r1", "r2"),
  type      = "norm")

p.vector <- c(A = .75, B = 1.25, t0 = .15, mean_v.true = 2.5, mean_v.false = 1.5)
print(model, p.vector)

ntrial <- 1e2
dat <- simulate(model, nsim = ntrial, ps = p.vector)
dmi <- BuildDMI(dat, model)

p.prior <- BuildPrior(
  dists = c("tnorm", "tnorm", "beta", "tnorm", "tnorm"),
  p1    = c(A = 1, B = 1, t0 = 1, mean_v.true = 1, mean_v.false = 1),
  p2    = c(1,  1,  1, 1, 1),
  lower = c(rep(0, 3),  rep(NA, 2)),
  upper = c(rep(NA, 2), 1, rep(NA, 2)))

sam0 <- run(StartNewsamples(5e2, dmi, p.prior))
sam  <- sam0
sam <- run(RestartSamples(5e2, sam, thin = 1))
rhats <- gelman(sam, verbose = TRUE)
p0 <- plot(sam)
p1 <- autocor(sam, nsubchain = 3)
p1 <- autocor(sam, nsubchain = 4)
p1 <- autocor(sam, nsubchain = 5)
d <- ConvertChains(sam, 1, 5e2, FALSE)


tmp0 <- ggmcmc::ac(d[Chain == 1 & Parameter == "A"]$value, 50)
tmp1 <- ac2(d[Chain == 1 & Parameter == "A"]$value, 50)
tmp2 <- ggdmc:::ac_(d[Chain == 1 & Parameter == "A"]$value, 50)
all(tmp0 == tmp1)
all(tmp0$Autocorrelation == tmp2[,1])
cbind(tmp0$Autocorrelation, tmp2[,1])
tmp0$Autocorrelation - tmp2[,1]
all.equal(tmp0$Autocorrelation, tmp2[,1])

library(rbenchmark)
res <- benchmark(r1 = ggmcmc::ac(d[Chain == 1 & Parameter == "A"]$value, 50),
                 r2 = ac2(d[Chain == 1 & Parameter == "A"]$value, 50),
                 r3 = ggdmc:::ac_(d[Chain == 1 & Parameter == "A"]$value, 50),
                 replications = 10)
print(res[,1:4])

nLags <- 50
wc.ac0 <- d[, .SD[, .(Lag = 1:nLags,
                      Autocorrelation = ggdmc:::ac_(value, nLags))],
            .(Parameter, Chain)]

library(dplyr)
wc.ac1 <- d %>% dplyr::group_by(Parameter, Chain) %>% dplyr::do(ggmcmc::ac(.$value,
                                                                  nLags))


tmp0 <- wc.ac0[wc.ac$Chain == 1 & wc.ac$Parameter == "A", ]$Autocorrelation
tmp1 <- wc.ac1[wc.ac$Chain == 1 & wc.ac$Parameter == "A", ]$Autocorrelation
all.equal(tmp0, tmp1)
dplyr::tbl_df(wc.ac0)
wc.ac1

