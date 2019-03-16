cat("-------------------- Testing BuildModel --------------------")

  rm(list = ls())

  pmap1 <- list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1",
                st0 = "1")
  pmap2 <- list(A = "S", B = "1", t0 = "1", mean_v = "M", sd_v = "1",
                st0 = "1")
  pmap3 <- list(A = "S", B = "D", t0 = "1", mean_v = "M", sd_v = "1",
                st0 = "1")
  pmap4 <- list(A.S = "S", B = "D", t0 = "1", mean_v = "M", sd_v = "1",
                st0 = "1")
  pmap5 <- list(A.S = "S", B = "D", t0 = "1", mean_v = "M", sd_v = "1",
                st0.S = "1")
  pmap6 <- list(A = "S", B = "D", t0 = "1", mean_v = "M", sd_v = "1",
                st0.S = "1")
  factors1 <- list(A = "1")
  factors2 <- list(S = c("s1", "s2"), F = c("f1", "f2"))
  factors3 <- list(S = c("s1", "s2"), F = c("f1", "f2", "f2"))
  factors4 <- list(S = c("s1", "s2"), F = c("f1", "f2", "f3"), R = c("r1", "r2"))
  factors5 <- list(S = c("s1", "s2"), F = c("f1", "f2", "f3"), D = c("d1", "d2"),
                   R = c("r1", "r2"))
  match.map1 <- list(M = list(s1 = 1, s2 = 2))
  constants <- NULL
  type <- "norm"
  responses <- c("r1", "r2")

  ## Normal functioning
  mapinfo1 <- ggdmc:::check_BuildModel(pmap1, responses, factors1, match.map1,
                                       constants, type)
  ## No need to supply factors
  mapinfo1 <- ggdmc:::check_BuildModel(pmap1, responses, NULL, match.map1,
                                       constants, type)
  ## Error: Must supply p.map
  expect_error(

  mapinfo1 <- ggdmc:::check_BuildModel(NULL, responses, factors1, match.map1,
                                       constants, type)
  )

  ## Error: Must supply responses
  expect_error(

  mapinfo1 <- ggdmc:::check_BuildModel(pmap1, NULL, factors1, match.map1,
                                       constants, type)
  )
  ## Error: All factors levels must be unqiue
  expect_error(

  mapinfo1 <- ggdmc:::check_BuildModel(pmap1, responses, factors3, match.map1,
                                       constants, type)
  )
  ## Error: Do not use 1, s, or R as a factor name
  expect_error(
  mapinfo1 <- ggdmc:::check_BuildModel(pmap1, responses, factors4, match.map1,
                                       constants, type)
  )
  expect_error(
  mapinfo1 <- ggdmc:::check_BuildModel(pmap1, responses, factors5, match.map1,
                                       constants, type)
  )
  ## Error: Dots not allowed in p.map names, fix: A.S
  expect_error(
  mapinfo1 <- ggdmc:::check_BuildModel(pmap4, responses, factors1, match.map1,
                                       constants, type)
  )
  ## Error: Dots not allowed in p.map names, fix: st0.S
  expect_error(
  mapinfo1 <- ggdmc:::check_BuildModel(pmap5, responses, factors1, match.map1,
                                       constants, type)
  )
  expect_error(
    mapinfo1 <- ggdmc:::check_BuildModel(pmap6, responses, factors1, match.map1,
                                       constants, type)
  )

  hasM <- sapply(pmap1, function(x){ any(x == "M")})

  ## Error: match.map is NULL, but found M factor
  expect_error(
  mapinfo <- ggdmc:::check_BuildModel(
    p_map     = pmap1,
    match_map = NULL,
    factors   = list(S = c("s1", "s2")),
    constants = c(st0 = 0, sd_v = 1),
    responses = c("r1", "r2"),
    type      = "norm")
  )

  pmap7 <- list(A = "1", B = "1", t0 = "1", mean_v = "1", sd_v = "1",
                st0 = "1")
  expect_error(
  model <- BuildModel(
    p.map     = pmap7,
    match_map = NULL,
    factors   = list(S = c("s1", "s2")),
    constants = c(st0 = 0, sd_v = 1),
    responses = c("r1", "r2"),
    type      = "norm")
  )



