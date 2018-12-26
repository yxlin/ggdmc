invlogit <- function(x) { 1 / (1+exp(-x)) }

dat <- data.frame(X = c(1.6907, 1.7242, 1.7552, 1.7842, 1.8113, 1.8369, 1.8610, 1.8839),
                  N = c(59, 60, 62, 56, 63, 59, 62, 60),
                  S = c(6, 13, 18, 28, 52, 53, 61, 60))
dat$F <- dat$N - dat$S
dat$logit <- log(dat$F/ dat$S)

m1 <- glm(cbind(S, F) ~ X, family = binomial, data = dat)
arm::display(m1)
## glm(formula = cbind(S, F) ~ X, family = binomial, data = dat)
##             coef.est coef.se
## (Intercept) -60.72     5.18 
## X            34.27     2.91 
## ---
##   n = 8, k = 2
##   residual deviance = 11.2, null deviance = 284.2 (difference = 273.0)
qchisq(.95, 6)

## Fitted values for y
tmp1 <- boot::inv.logit(predict(m1)) * dat$N
tmp2 <- stats::plogis(predict(m1)) * dat$N
tmp3 <- ggdmc:::invlogit(predict(m1)) * dat$N
tmp4 <- ggdmc:::invlogit2(predict(m1)) * dat$N
all.equal(tmp1, tmp3[,1])
cbind(tmp1, tmp2, tmp3[,1], tmp4[,1])
all(tmp1 == tmp3[,1])

require(microbenchmark)
res <- microbenchmark(
  boot::inv.logit(predict(m1)),
  stats::plogis(predict(m1)),
                      ggdmc:::invlogit(predict(m1)),
  ggdmc:::invlogit2(predict(m1))
)
res
mat1 <- as.matrix(dat[,c("S", "N")])
ggdmc:::test_mat(mat1)
