cat("\n-------------------- Testing GetParameterMatrix --------------------")

rm(list = ls())
## Used in simulate.model and simulate_many
## Does not use p_df
model <- BuildModel(
  p.map     = list(a ="1", v = "1",z = "1", d = "1", sz = "1", sv = "1",
                   t0 = "1", st0 = "1"),
  match.map = list(M = list(s1 = "r1", s2 ="r2")),
  factors   = list(S = c("s1", "s2")),
  constants = c(st0 = 0, d = 0),
  responses = c("r1", "r2"),
  type      = "rd")

p.prior <- BuildPrior(
  dists = c("tnorm", "tnorm", "beta", "beta", "tnorm", "beta"),
  p1    = c(a = 1, v = 0, z = 1, sz = 1, sv = 1, t0 = 1),
  p2    = c(a = 1, v = 2, z = 1, sz = 1, sv = 1, t0 = 1),
  lower = c(0, -5, NA, NA, 0, NA),
  upper = c(2,  5, NA, NA, 2, NA))

## Example 1: Randomly generate 2 sets of true parameters from
## p.prior
GetParameterMatrix(model, 2, p.prior)
##            a         v         z        sz       sv        t0
## [1,] 1.963067  1.472940 0.9509158 0.5145047 1.344705 0.0850591
## [2,] 1.512276 -1.995631 0.6981290 0.2626882 1.867853 0.1552828

## Example 2: Use a user-selected true parameters
true.vector  <- c(a=1, v=1, z=0.5, sz=0.2, sv=1, t0=.15)
GetParameterMatrix(model, 2, NA, true.vector)
##      a v   z  sz sv   t0
## [1,] 1 1 0.5 0.2  1 0.15
## [2,] 1 1 0.5 0.2  1 0.15
GetParameterMatrix(model, 2, ps = true.vector)

## Example 3: When a user enters arbritary sequence of parameters.
## The user enters sv before sz.
## GetParameterMatrix will rearrange the sequence.
true.vector  <- c(a=1, v=1, z=0.5, sv=1, sz = .2, t0=.15)
GetParameterMatrix(model, 2, NA, true.vector)
##      a v   z  sz sv   t0
## [1,] 1 1 0.5 0.2  1 0.15
## [2,] 1 1 0.5 0.2  1 0.15



