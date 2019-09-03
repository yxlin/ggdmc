## Basic test ----
fptpdf <- function(z,x0max,chi,v,sdv) {
  if (x0max==0) return( (chi/z^2)*dnorm(chi/z,mean=v,sd=sdv))
  zs=z*sdv ; zu=z*v ; chiminuszu=chi-zu
  chizu=chiminuszu/zs ; chizumax=(chiminuszu-x0max)/zs
  (v*(pnorm(chizu)-pnorm(chizumax)) +
      sdv*(dnorm(chizumax)-dnorm(chizu)))/x0max
}

fptcdf <- function(z,x0max,chi,v,sdv) {
  if (x0max==0) return(pnorm(chi/z,mean=v,sd=sdv,lower.tail=F))
  zs=z*sdv ; zu=z*v ; chiminuszu=chi-zu ; xx=chiminuszu-x0max
  chizu=chiminuszu/zs ; chizumax=xx/zs
  tmp1=zs*(dnorm(chizumax)-dnorm(chizu))
  tmp2=xx*pnorm(chizumax)-chiminuszu*pnorm(chizu)
  1+(tmp1+tmp2)/x0max
}

## PDF --------
mean_v <- 2.4
A <- 1.2
b <- 2.7
t0 <- .2
sd_v <- 1
st0 <- 0
posdrift <- TRUE


pnorm(mean_v/sd_v)

res0 <- ggdmc252:::fptpdf(.3, A, b, mean_v, sd_v, t0, posdrift)
res1 <- ggdmc:::fptpdf(.3, A, b, mean_v, sd_v, t0, posdrift)
res2 <- rtdists:::dlba_norm_core(.3, A, b, t0, mean_v, sd_v)
res3 <- fun(.3, A, b, t0, mean_v, sd_v)
# rtdists:::pnormP
# rtdists:::dnormP

res0
res1
res2
res3


res0[,1]==res1[,1]

A <- .02
b <- .01
mean_v <- c(2.4, 2.2, 1.5)
sd_v   <- c(1,     1, 1.5)
t0 <- .01

res0 <- ggdmc252:::fptpdf(0, A, b, mean_v[1], sd_v[1], t0, posdrift)
res1 <- ggdmc:::fptpdf(0, A, b, mean_v[1], sd_v[1], t0, posdrift)
res2 <- fptpdf(0-t0, x0max=A, chi=b, v=mean_v[1], sdv=sd_v[1])
res3 <- rtdists:::dlba_norm_core(0, A, b, t0, mean_v[1], sd_v[1])
# rtdists:::dlba_norm
res0[,1]
res1[,1]
res2
res3


res0 <- ggdmc252:::fptcdf(0, A, b, mean_v[1], sd_v[1], t0, posdrift)
res1 <- ggdmc:::fptcdf(0, A, b, mean_v[1], sd_v[1], t0, posdrift)
res2 <- fptcdf(0-t0, x0max=A, chi=b, v=mean_v[1], sdv=sd_v[1])
res3 <- rtdists:::plba_norm_core(0, A, b, t0, mean_v[1], sd_v[1])

res0[,1]
res1[,1]
res2
res3



rt <- seq(0, 10, .01)
res0 <- ggdmc252:::fptpdf(rt, A, b, mean_v[1], sd_v[1], t0, posdrift)
res1 <- ggdmc:::fptpdf(rt, A, b, mean_v[1], sd_v[1], t0, posdrift)
res2 <- fptpdf(rt-t0, x0max=A, chi=b, v=mean_v[1], sdv=sd_v[1])
all(res0[,1]==res1[,1])
all(res1[,1]==res2)

all.equal(res0[,1], res1[,1])
all.equal(res1[,1], res2)

head(cbind(res0, res1, res2, rt))
tail(cbind(res0, res1, res2, rt))

rt <- seq(.1, 5, .01)
mean_v <- seq(0, 5, .1)
A <- seq(0, 5, .1)
b <- 2.7
t0 <- .1
sd_v <- 1


for(i in 1:length(mean_v))
{
  for(j in 1:length(A))
  {
    res0 <- ggdmc252:::fptpdf(rt, A[j], b, mean_v[i], sd_v, t0, posdrift)
    res1 <- ggdmc:::fptpdf(rt, A[j], b, mean_v[i], sd_v, t0, posdrift)
    test <- all(res0[,1]==res1[,1])

    if(!test) cat("[", mean_v[i], " ", A[j], " ", b, " ", t0, " ", sd_v,
                  "]", " results in a different likelihood\n")
  }
}




# rt <- seq(0, 5, .01)
# mean_v <- seq(0, 5, .1)
# A <- seq(0, 5, .1)
# b <- seq(0, 5, .1)
# t0 <- seq(0, 2, .1)
# sd_v <- 1
#
# for(i in 1:length(mean_v))
# {
#   for(j in 1:length(A))
#   {
#     for(k in 1:length(b))
#     {
#       for(l in 1:length(t0))
#       {
#         res0 <- ggdmc252:::fptpdf(rt, A[j], b[k], mean_v[i], sd_v, t0[l],
#                                   posdrift)
#         res1 <- ggdmc:::fptpdf(rt, A[j], b[k], mean_v[i], sd_v, t0[l],
#                                posdrift)
#         test <- all(res0[,1]==res1[,1])
#         if(!test) cat("[", mean_v[i], " ", A[j], " ", b[k], " ", t0[l], " ", sd_v[m],
#                       "]", " results in a different likelihood\n")
#       }
#     }
#   }
# }


## CDF --------

mean_v <- 2.4
A <- 1.2
b <- 2.7
t0 <- .2
sd_v <- 1
st0 <- 0
posdrift <- FALSE


res0 <- ggdmc252:::fptcdf(.6, A, b, mean_v, sd_v, t0, posdrift)
res1 <- ggdmc:::fptcdf(.6, A, b, mean_v, sd_v, t0, posdrift)
res2 <- fptcdf(z=.6-t0,  x0max=A, chi=b, v=mean_v, sdv=sd_v)
res0[,1]==res1[,1]
all.equal(res1[,1], res2)
cbind(res1[,1], res2)


rt <- seq(0, 10, .01)
res0 <- ggdmc252:::fptcdf(rt, A, b, mean_v, sd_v, t0, posdrift)
res1 <- ggdmc:::fptcdf(rt, A, b, mean_v, sd_v, t0, posdrift)
res2 <- fptcdf(rt-t0, x0max=A, chi=b, v=mean_v, sdv=sd_v)

all(res0[,1]==res1[,1])
all(res1[,1]==res2)
all.equal(res1[,1], res2)

head(cbind(res0, res1, res2, rt), 22)
tail(cbind(res0, res1, res2, rt))

head(cbind(res0, res1, res2, rt))
tail(cbind(res0, res1, res2, rt))

all.equal(res0[,1], res1[,1])
all.equal(res0[,1], res2)


rt <- seq(.2, .6, .01)
mean_v <- seq(1.2, 5, .01)
for(j in 1:length(mean_v))
{
  res0 <- ggdmc252:::fptcdf(rt, A, b, mean_v[j], sd_v, t0, posdrift)
  res1 <- ggdmc:::fptcdf(rt, A, b, mean_v[j], sd_v, t0, posdrift)

  if( !all.equal(res0[,1], res1[,1]) )
  {
    cat("Not the same \n")
  }
}






