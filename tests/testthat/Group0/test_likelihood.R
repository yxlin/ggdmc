cat("\n-------------------- Testing likelihood --------------------")

## Wiener  ----------
cat("\nWiener diffusion model: \n")

rm(list = ls())
model <- BuildModel(
  p.map     = list(a = "1", v="1", z="1", d="1", sz="1", sv="1", t0="1", st0="1"),
  match.map = list(M = list(s1 = "r1", s2 = "r2")),
  factors   = list(S = c("s1", "s2")),
  responses = c("r1","r2"),
  constants = c(st0 = 0, d = 0, sv = 0, sz = 0),
  type      = "rd")

p.vector <- c(a=1, v=1.5, z=0.5, t0=.15)
dat <- simulate(model, nsim = 10, ps = p.vector)
dmi <- BuildDMI(dat, model)
res0 <- likelihood(p.vector, dmi)
cat( round(res0, 2), "\n")

## DDM  ----------
cat("\nDecision diffusion model: \n")

model <- BuildModel(
  p.map     = list(a = "1", v = "1", z = "1", d = "1", sz = "1", sv = "1",
                   t0 = "1", st0 = "1"),
  match.map = list(M = list(s1 = "r1", s2 = "r2")),
  factors   = list(S = c("s1", "s2")),
  responses = c("r1", "r2"),
  constants = c(st0 = 0, d = 0),
  type      = "rd")
p.vector <- c(a = 1, v = 1.2, z = .38, sz = .25, sv = .2, t0 = .15)

dat <- simulate(model, nsim = 10, ps = p.vector)
dmi <- BuildDMI(dat, model)
res0 <- likelihood(p.vector, dmi)
cat( round(res0, 2), "\n")

## LBA --------
cat("\nLBA model: \n")
model <- BuildModel(
  p.map     = list(A = "1", B = "R", t0 = "1", mean_v = c("D", "M"),
                   sd_v = "M", st0 = "1"),
  match.map = list(M = list(s1 = 1, s2 = 2)),
  factors   = list(S = c("s1", "s2"), D = c("d1", "d2")),
  constants = c(sd_v.false = 1, st0 = 0),
  responses = c("r1", "r2"),
  type      = "norm")

p.vector <- c(A=.51, B.r1=.69, B.r2=.88, t0=.24, mean_v.d1.true=1.1,
              mean_v.d2.true=1.0, mean_v.d1.false=.34, mean_v.d2.false=.02,
              sd_v.true=.11)

dat <- simulate(model, nsim = 10, ps = p.vector)
dmi <- BuildDMI(dat, model)
res0 <- likelihood(p.vector, dmi)
cat( round(res0, 2) , "\n")


## LBA --------
cat("\nLikelihood example 1: \n")

 model <- BuildModel(
 p.map     = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1",
             st0 = "1"),
 match.map = list(M = list(s1 = 1, s2 = 2)),
 factors   = list(S = c("s1", "s2")),
 constants = c(st0 = 0, sd_v = 1),
 responses = c("r1", "r2"),
 type      = "norm")
 p.vector <- c(A = .25, B = .35,  t0 = .2, mean_v.true = 1, mean_v.false = .25)
 dat <- simulate(model, 1e3,  ps = p.vector)
 dmi <- BuildDMI(dat, model)
 den <- likelihood(p.vector, dmi)


cat("\nLikelihood example 2: \n")
 model <- BuildModel(
 p.map     = list(a = "1", v = "1", z = "1", d = "1", t0 = "1", sv = "1",
             sz = "1", st0 = "1"),
 constants = c(st0 = 0, d = 0),
 match.map = list(M = list(s1 = "r1", s2 = "r2")),
 factors   = list(S = c("s1", "s2")),
 responses = c("r1", "r2"),
 type      = "rd")

 p.vector <- c(a = 1, v = 1, z = 0.5, sz = 0.25, sv = 0.2, t0 = .15)
 dat <- simulate(model, 1e2, ps = p.vector)
 dmi <- BuildDMI(dat, model)
 den <- likelihood (p.vector, dmi)

