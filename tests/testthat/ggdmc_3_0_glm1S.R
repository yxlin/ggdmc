## Sampling ------------
require(ggdmc)
setwd("/media/yslin/KIWI/Documents/ggdmc_lesson")
rm(list = ls())
ntrial <- 1000
p.vector <- c(alpha = 242.7, beta = 6.185, epsilon = 6.086)
x <- rep(c(8, 15, 22, 29, 36), each = ntrial/5)
mu <- p.vector[1] + p.vector[2] * x
y  <- rnorm(ntrial, mu, p.vector[3])

npar <- 3
p.prior  <-BuildPrior(
  dists = rep("tnorm", 3),
  p1    = c(alpha = 200, beta = 0, epsilon = 1),
  p2    = c(alpha = 100, beta = 5, epsilon = 5),
  lower = rep(NA, 3),
  upper = rep(NA, 3))
plot(p.prior, ps = p.vector)

A  <- rep(c("a1", "a2"), each = ntrial/2)
dat <- data.frame(A = A, R = rep("r1", ntrial), X = x, RT = y)
dplyr::tbl_df(dat)
fit.lm <- lm(RT ~ X, data = dat)
summary(fit.lm)
coef(fit.lm)
model <- BuildModel(
  p.map     = list(alpha = "1", beta = "1", epsilon = "1"),
  match.map = NULL,
  factors   = list(X = c("x1")),
  responses = "r1",
  constants = NULL,
  type      = "glm")


pnames   <- names(attr(model, "p.vector"))
allpar   <- attr(model, "all.par")
parnames <- attr(model, "par.names")
type     <- attr(model, "type")
n1       <- attr(model, "n1.order")
resp     <- attr(model, "responses")
cell     <- ggdmc:::check_cell(1, model)
isr1     <- ggdmc:::check_rd(type, model)

out <- ggdmc:::p_df(p.vector, cell, pnames, allpar, parnames,
                          model, type, dimnames(model)[[1]], dimnames(model)[[2]],
                          dimnames(model)[[3]], isr1, n1, TRUE)


TableParameters(p.vector, 1, model, FALSE)
TableParameters(p.vector, 2, model, FALSE)

dmi <- BuildDMI(dat, model)

setwd("/media/yslin/KIWI/Documents/ggdmc_lesson/")
path <- c("data/ggdmc_glm.rda")
fit0 <- run(StartNewsamples(5e2, dmi, p.prior))
fit  <- fit0
thin <- 4
repeat {
  fit <- run(RestartSamples(5e2, fit, thin = thin))
  save(fit0, fit, file = path[1])
  rhat <- gelman(fit, verbose = TRUE)
  if (all(rhat$mpsrf < 1.1)) break
  thin <- thin * 2
}
cat("Done ", path[1], "\n")
setwd("/media/yslin/KIWI/Documents/ggdmc/")


p0 <- plot(fit0)
p1 <- plot(fit0, start = 201)

library(gridExtra)
grid.arrange(p0, p1, ncol = 2)

p1 <- plot(fit)
p2 <- plot(fit, pll = FALSE, den = TRUE)
est <- summary(fit, recover = TRUE, ps = p.vector, verbose = TRUE)
est
# p3 <- ggdmc::autocor(sam)
# p4 <- pairs(sam)

coef(fit.lm)
