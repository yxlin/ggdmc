cat("\n-------------------- Testing tnorm --------------------")

## rtn example
dat1 <- rtnorm(1e5, 0, 1, 0, Inf)

## dtn example
x <- seq(-5, 5, length.out = 1e3)
dat1 <- dtnorm(x, 0, 1, -2, 2, 0)

## ptn example
x <- seq(-10, 10, length.out = 1e2)
mean <- 0
sd <- 1
lower <- 0
upper <- 5
dat1 <- ptnorm(x, 0, 1, 0, 5, lg = TRUE)
