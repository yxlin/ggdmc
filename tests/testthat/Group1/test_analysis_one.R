cat("\n-------------------- Testing analysis_one (Wiener)--------------------")

rm(list = ls())

## Set up model ----
cat("\nBuildModel report: \n")
model <- BuildModel(
  p.map     = list(a = "1", v="1", z="1", d="1", sz="1", sv="1", t0="1", st0="1"),
  match.map = list(M = list(s1 = "r1", s2 = "r2")),
  factors   = list(S = c("s1", "s2")),
  responses = c("r1","r2"),
  constants = c(st0 = 0, d = 0, sv = 0, sz = 0),
  type      = "rd")

npar <- length(GetPNames(model))
p.vector <- c(a=1, v=1.5, z=0.5, t0=.15)
dat <- simulate(model, nsim = 100, ps = p.vector)
dmi <- BuildDMI(dat, model)

p.prior <- BuildPrior(
  dists = rep("tnorm", npar),
  p1=c(a=1, v=0, z=1, t0=1),
  p2=c(a=1, v=2, z=1, t0=1),
  lower = c(0, -5, rep(0, 2)),
  upper = rep(NA, npar))

## Sampling and check model ----
cat("Starting a new model fit: \n")
fit <- run(StartNewsamples(dmi, p.prior, block = FALSE), block = FALSE)
res <- gelman(fit, verbose=TRUE)

pdf(file = "analysis_one.pdf")
p0 <- plot(fit)
p1 <- plot(fit, den = TRUE)
p2 <- plot(fit, pll=FALSE)
p3 <- plot(fit, pll=FALSE, den = TRUE)
dev.off()

## Analysis -----------
cat("Reporting the result of a model fit: \n")
est <- summary(fit, recovery = TRUE, ps = p.vector, verbose = TRUE)
tmp1 <- theta2mcmclist(fit)
tmp2 <- theta2mcmclist(fit, start = 10, end = 90)
tmp3 <- theta2mcmclist(fit, split = TRUE)

cat("\nWhich chains were selected: \n")
tmp4 <- theta2mcmclist(fit, subchain = TRUE)
tmp5 <- theta2mcmclist(fit, subchain = TRUE, nsubchain = 4)
tmp6 <- theta2mcmclist(fit, thin = 2)
