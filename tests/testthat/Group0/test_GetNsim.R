cat("-------------------- Testing GetNsim --------------------")

  model <- BuildModel(
    p.map     = list(A = "1", B = "R", t0 = "1", mean_v = "M", sd_v = "M",
                     st0 = "1"),
    match.map = list(M = list(s1 = 1, s2 = 2)),
    constants = c(sd_v.false = 1, st0 = 0),
    factors   = list(S = c("s1","s2")),
    responses = c("r1", "r2"),
    type      = "norm")


  #######################30
  ## Example 1
  #######################30
  GetNsim(model, ns = 2, n = 1)
  #      [,1] [,2]
  # [1,]    1    1
  # [2,]    1    1

  #######################30
  ## Example 2
  #######################30
  n <- matrix(c(1:2), ncol = 1)
  #      [,1]
  # [1,]    1  ## subject 1 has 1 response for each cell
  # [2,]    2  ## subject 2 has 2 responses for each cell

  GetNsim(model, ns = 2, n = n)
  #      [,1] [,2]
  # [1,]    1    1
  # [2,]    2    2

  #######################30
  ## Example 3
  #######################30
  n <- matrix(c(1:2), nrow = 1)
  #      [,1] [,2]
  # [1,]    1    2
  GetNsim(model, ns = 2, n = n)
  #     [,1] [,2]
  # [1,]   1    2 ## subject 1 has 1 response for cell 1 and 2 responses for cell 2
  # [2,]   1    2 ## subject 2 has 1 response for cell 1 and 2 responses for cell 2

  #######################30
  ## Example 4
  #######################30
  n <- matrix(c(1:4), nrow=2)
  #      [,1] [,2]
  # [1,]    1    3
  # [2,]    2    4
  ggdmc::GetNsim(model, ns = 2, n = n)
  #      [,1] [,2]
  # [1,]    1    3 ## subject 1 has 1 response for cell 1 and 3 responses for cell 2
  # [2,]    2    4 ## subject 2 has 2 responses for cell 1 and 4 responses for cell 2




