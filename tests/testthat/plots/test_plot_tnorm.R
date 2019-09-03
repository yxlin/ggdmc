cat("\n-------------------- Testing plot tnorm --------------------")

## rtn example
pdf(file = "tnorm.pdf")

dat1 <- rtnorm(1e5, 0, 1, 0, Inf)
hist(dat1, breaks = "fd", freq = FALSE, xlab = "",
     main = "Truncated normal distributions")

## dtn example
x <- seq(-5, 5, length.out = 1e3)
dat1 <- dtnorm(x, 0, 1, -2, 2, 0)
plot(x, dat1, type = "l", lwd = 2, xlab = "", ylab= "Density",
     main = "Truncated normal distributions")

## ptn example
x <- seq(-10, 10, length.out = 1e2)
mean <- 0
sd <- 1
lower <- 0
upper <- 5
dat1 <- ptnorm(x, 0, 1, 0, 5, lg = TRUE)
