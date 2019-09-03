## Two accumulators -------------

mean_v <- matrix(c(2.4, 2.2)); mean_v
A <- 1.2
b <- 2.7
t0 <- .2
sd_v <- c(1, 1)
x <- seq(0, 3, .4);
posdrift <- TRUE
nv <- length(mean_v)

res0 <- ggdmc252:::n1PDFfixedt0(x, A, b, mean_v, sd_v, t0, posdrift)
res1 <- ggdmc:::n1PDFfixedt0(x, A = rep(A, nv), b = rep(b, nv),
                             mean_v=mean_v[,1], sd_v=sd_v,
                             t0=rep(t0, nv), posdrift)
all(res0 == res1)

cbind(res0[,1], res1[,1])
all.equal(res0[,1], res1[,1])

rtdists:::n1PDFfixedt0

 require(microbenchmark)
 res <- microbenchmark(
   res0 <- ggdmc252:::n1PDFfixedt0(x, A, b, mean_v, sd_v, t0, posdrift),
   res1 <- ggdmc:::n1PDFfixedt0(x, mean_v[,1], sd_v, A, b, t0, posdrift),
   times = 1000L)

res

## Three accumulators -------------
mean_v <- matrix(c(2.4, 2.2, 1.5)); mean_v
A <- 1.2
b <- 2.7
t0 <- .2
sd_v <- c(1, 1, 1.5)
x <- seq(0, 3, .4);
posdrift <- TRUE
nv <- 3

res0 <- ggdmc252:::n1PDFfixedt0(x, A, b, mean_v, sd_v, t0, posdrift)
res1 <- ggdmc:::n1PDFfixedt0(x, A, b, mean_v=mean_v[,1], sd_v=sd_v,
                             t0=t0, posdrift)
all(res0 == res1)



## Bad parameters

## pmat: A        b       t0   mean_v sd_v st0
##  -0.9817   2.2475  -0.1214   2.5698   1.0000        0
##  -0.9817   2.2475  -0.1214  -0.9132   1.0000        0
A <- -0.9817
b <- 2.2475
t0 <- -0.1214
mean_v <- c(2.5698, -0.9132)
sd_v <- c(1,1)
x <- seq(0, 3, .4);
posdrift <- TRUE


res0 <- ggdmc252:::n1PDFfixedt0(x, A, b, mean_v, sd_v, t0, posdrift)
res1 <- ggdmc:::n1PDFfixedt0(x, A, b, mean_v=mean_v, sd_v=sd_v,
                             t0=t0, posdrift)

all(res0==res1)


mean_v <- matrix(c(2.4, 2.2, 1.5)); mean_v
nv <- length(mean_v)
b <- 2.7
sd_v <- c(1, 1, 1.5)
posdrift <- TRUE

A <- seq(0.01, 5, .01)
b <- seq(0.01, 5, .01)
t0 <- seq(.01, .2, .01)
dt <- seq(0, 3, .4);


for(i in 1:length(A))
{
  for(j in 1:length(b))
  {
    for(k in 1:length(t0))
    {
      dt <- dt+t0[k]
      res0 <- ggdmc252:::n1PDFfixedt0(dt, A[i], b[j], mean_v, sd_v, t0[k], posdrift)
      res1 <- ggdmc:::n1PDFfixedt0(dt, A[i], b[j], mean_v[,1], sd_v, t0[k], posdrift)

      if( !all(res0[,1]==res1[,1]) )
      {
        cat("not equal, i: " , i, "\n")
        stop("")
      }
    }
  }
}


round(res0, 2)



# b = A + B




# den3 <- ggdmc::n1PDFfixedt0_pda(x, A, b, mean_v, sd_v, t0, n, h, debug)
# den4 <- gpda::n1PDF(x, A, b, mean_v, sd_v, t0, n, h, debug = debug)
# plot(x, res0[,1], type="l")
# lines(x, res1, lwd = 2)
# lines(x, res2[,1], lwd = 2)
#
# lines(x, den3, lwd= 3, col = wesanderson::wes_palettes$BottleRocket[1])
# lines(x, den4, lwd= 3, col = wesanderson::wes_palettes$BottleRocket[2])
#
# library(rbenchmark)
# res <- benchmark(r1 = ggdmc::n1PDFfixedt0(x, A, b, mean_v, sd_v, t0),
#                  r2 = rtdists::n1PDF(x, A, b, mean_v = as.vector(mean_v), sd_v = sd_v, t0, silent = T),
#                  r3 = ggdmc::n1PDFfixedt0_pda(x, A, b, mean_v, sd_v, t0, n, h, debug),
#                  r4 = gpda::n1PDF(x, A, b, mean_v, sd_v, t0, n, h, debug = debug),
#                  replications = 10)
# print(res[,1:4])
#
#
# x1 <- ggdmc::rlba_norm(n, A, b, t0, mean_v, sd_v)
# x2 <- ggdmc::rlba(n, A, b, t0, mean_v, sd_v)
#
# x1RT1 <- x1[x1[,2] == 1, 1]
# x1RT2 <- x1[x1[,2] == 2, 1]
# x2RT1 <- x2[x2[,2] == 1, 1]
# x2RT2 <- x2[x2[,2] == 2, 1]
#
# trim <- 5
# x1RT1T <- x1RT1[x1RT1 < trim]
# x1RT2T <- x1RT2[x1RT2 < trim]
# x2RT1T <- x2RT1[x2RT1 < trim]
# x2RT2T <- x2RT2[x2RT2 < trim]
#
# par(mfrow=c(2,2))
# hist(x1RT1T, breaks = "fd", freq = FALSE)
# hist(x1RT2T, breaks = "fd", freq = FALSE)
# hist(x2RT1T, breaks = "fd", freq = FALSE)
# hist(x2RT2T, breaks = "fd", freq = FALSE)
#
# par(mfrow=c(1, 1))
# hist(x1RT1T, breaks = "fd", freq = FALSE, col = "lightgrey")
# hist(x2RT1T, breaks = "fd", freq = FALSE, add = TRUE, col = "blue")
# hist(x1RT2T, breaks = "fd", freq = FALSE)
# hist(x2RT2T, breaks = "fd", freq = FALSE, add = TRUE, col = "blue")
#
#
#
# drifts <- clba::make_v(n, A, b, t0, mean_v, sd_v, st0, posdrift, return.ttf)
# dat0 <- clba::make_r(drifts, b, A, n_v, t0, st0, n, return.ttf)
# dat1 <- clba:::make.r(drifts, b, A, n_v, t0, st0, n, return.ttf)
# dat3 <- clba:::rlba.norm(n, A, b, t0, mean_v, sd_v, st0, posdrift, return.ttf)
