cat("\n-------------------- Testing rprior --------------------")

## rprior  ----------
rm(list = ls())

p.prior <- BuildPrior(
 dists = c("tnorm", "tnorm", "beta", "tnorm", "beta", "beta"),
 p1    = c(a = 1, v = 0, z = 1, sz = 1, sv = 1, t0 = 1),
 p2    = c(a = 1, v = 2, z = 1, sz = 1, sv = 1, t0 = 1),
 lower = c(0,-5, NA, NA, 0, NA),
 upper = c(2, 5, NA, NA, 2, NA))

ggdmc:::rprior_mat(p.prior, 1)
ggdmc:::rprior_mat(p.prior, 2)

rprior(p.prior, 9)




