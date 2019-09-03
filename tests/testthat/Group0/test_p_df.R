cat("\n-------------------- Testing p_df --------------------")

## Wiener  ----------
rm(list = ls())
model <- BuildModel(
  p.map     = list(a = "1", v="1", z="1", d="1", sz="1", sv="1", t0="1", st0="1"),
  match.map = list(M = list(s1 = "r1", s2 = "r2")),
  factors   = list(S = c("s1", "s2")),
  responses = c("r1","r2"),
  constants = c(st0 = 0, d = 0, sv = 0, sz = 0),
  type      = "rd")

p.vector <- c(a=1, v=1.5, z=0.5, t0=.15)

type     <- attr(model, "type")
pnames   <- names(attr(model, "p.vector"))
(parnames <- attr(model, "par.names"))
dim0 <- dimnames(model)[[1]]
dim1 <- dimnames(model)[[2]]
dim2 <- dimnames(model)[[3]]
(allpar   <- attr(model, "all.par"));
(isr1     <- ggdmc:::check_rd(type, model))
n1idx    <- attr(model, "n1.order"); n1idx
n1order <- TRUE

res0 <- TableParameters(p.vector, 1, model, FALSE)
res0

for(i in 1:length(dim0))
{
  cell <- dim0[i]
  res0 <- ggdmc:::p_df(p.vector, cell, type, pnames, parnames, dim0, dim1, dim2,
                       allpar, model, isr1, n1idx, n1order)
  print(res0)
}

## DDM  ----------
model <- BuildModel(
  p.map     = list(a = "1", v = "1", z = "1", d = "1", sz = "1", sv = "1",
                   t0 = "1", st0 = "1"),
  match.map = list(M = list(s1 = "r1", s2 = "r2")),
  factors   = list(S = c("s1", "s2")),
  responses = c("r1", "r2"),
  constants = c(st0 = 0, d = 0),
  type      = "rd")
p.vector <- c(a = 1, v = 1.2, z = .38, sz = .25, sv = .2, t0 = .15)

type     <- attr(model, "type")
pnames   <- names(attr(model, "p.vector"))
(parnames <- attr(model, "par.names"))
dim0 <- dimnames(model)[[1]]
dim1 <- dimnames(model)[[2]]
dim2 <- dimnames(model)[[3]]
(allpar   <- attr(model, "all.par"));
(isr1     <- ggdmc:::check_rd(type, model))
n1idx    <- attr(model, "n1.order")
n1order <- TRUE

res0 <- TableParameters(p.vector, 1, model, FALSE)
res0

for(i in 1:length(dim0))
{
  cell <- dim0[i]
  res0 <- ggdmc:::p_df(p.vector, cell, type, pnames, parnames, dim0, dim1, dim2,
                       allpar, model, isr1, n1idx, n1order)
  print(res0)
}


## LBA --------
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

type     <- attr(model, "type")
pnames   <- names(attr(model, "p.vector"))
(parnames <- attr(model, "par.names"))
dim0 <- dimnames(model)[[1]]
dim1 <- dimnames(model)[[2]]
dim2 <- dimnames(model)[[3]]
(allpar   <- attr(model, "all.par"));
(isr1     <- ggdmc:::check_rd(type, model))
n1idx    <- attr(model, "n1.order")
n1order <- TRUE

res0 <- TableParameters(p.vector, 1, model, FALSE)
res0

for(i in 1:length(dim0))
{
  cell <- dim0[i]
  res0 <- ggdmc:::p_df(p.vector, cell, type, pnames, parnames, dim0, dim1, dim2,
                       allpar, model, isr1, n1idx, n1order)
  print(res0)
}



