require(ggdmc); require(testthat); require(ggplot2); require(data.table)
context("Post-predictive")

test_that("Post-predictive", {
  rm(list = ls())
  setwd("/media/yslin/KIWI/Documents/ggdmc")
  model <- BuildModel(
    p.map     = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1",
                     st0 = "1"),
    match.map = list(M = list(s1 = 1, s2 = 2)),
    factors   = list(S = c("s1", "s2")),
    constants = c(st0 = 0, sd_v = 1),
    responses = c("r1", "r2"),
    type      = "norm")


  # p.map     = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1",
  #                  st0 = "1")
  # match.map = list(M = list(s1 = 1, s2 = 2))
  # factors   = list(S = c("s1", "s2"))
  # constants = c(st0 = 0, sd_v = 1)
  # responses = c("r1", "r2")
  # type      = "norm"
  # posdrift   = TRUE
  # verbose    = TRUE
  
  p.vector <- c(A = .75, B = 1.25, t0 = .15, mean_v.true = 2.5, mean_v.false = 1.5)
  ntrial <- 5e2
  dat <- simulate(model, nsim = ntrial, ps = p.vector)
  dmi <- BuildDMI(dat, model)
  # n1.order   <- attr(model, "n1.order")
  # row.names(n1.order)
  # dimnames(model)[[1]]

  p.prior <- BuildPrior(
    dists = c("tnorm", "tnorm", "beta", "tnorm", "tnorm"),
    p1    = c(A = 1, B = 1, t0 = 1, mean_v.true = 1, mean_v.false = 1),
    p2    = c(1,  1,  1, 1, 1),
    lower = c(rep(0, 3),  rep(NA, 2)),
    upper = c(rep(NA, 2), 1, rep(NA, 2)))
  plot(p.prior, ps = p.vector)

  # res <- ggdmc:::check_BuildDMI(dat, model)
  # subject_models <- res$issm
  # modeli <- res$model
  # fnams <- names(attr(modeli, "factors"))
  # res <- ggdmc:::check_BuildDMI(dat, model)
  # subject_models <- res$issm
  # modeli <- res$model
  # fnams <- names(attr(modeli, "factors"))
  # cells <- apply(dat[, c(fnams, "R")], 1, paste, collapse = ".")
  # str(cells)
  
  
  # cell.index <- vector("list", dim(model)[1])
  # dim1 <- dimnames(model)[[1]]
  # names(cell.index) <- dim1
  # str(cell.index)
  
  # for ( j in names(cell.index) ) cell.index[[j]] <- cells %in% j
  # str(cell.index)

  ## Sampling ---------
  # setwd("/media/yslin/KIWI/Documents/ggdmc_lesson")
  path <- c("data/Lesson3/ggdmc_3_2_LBA1S_1e4_tmp.rda")
  fit0 <- run(StartNewsamples(5e2, dmi, p.prior))
  fit  <- fit0
  thin <- 1
  repeat {
    fit <- run(RestartSamples(5e2, fit, thin = thin))
    save(fit0, fit, file = path[1])
    rhat <- gelman(fit, verbose = TRUE)
    if (all(rhat$mpsrf < 1.1)) break
    thin <- thin * 2
  }
  cat("Done ", path[1], "\n")
  setwd("/media/yslin/KIWI/Documents/ggdmc")
  p0 <- plot(fit)
  p1 <- plot(fit, start = 101)
  p2 <- plot(fit0, start = 201)
  
  # gridExtra::grid.arrange(p0, p1, p2, ncol = 3)
  # d <- data.table(fit$data)
  # d[, .N, .(S)]
  
  ## Analysis -----------
  est <- summary(fit, start = 1, recovery = TRUE, ps = p.vector, verbose = TRUE)
  pp <- ggdmc:::predict_one(fit, xlim = c(0, 5))
  range(pp$RT)

  dat$C <- ifelse(dat$S == "s1"  & dat$R == "r1",  TRUE,
           ifelse(dat$S == "s2" & dat$R == "r2", TRUE,
           ifelse(dat$S == "s1"  & dat$R == "r2", FALSE,
           ifelse(dat$S == "s2" & dat$R == "r1",  FALSE, NA))))

  pp$C <- ifelse(pp$S == "s1"  & pp$R == "r1",  TRUE,
          ifelse(pp$S == "s2" & pp$R == "r2", TRUE,
          ifelse(pp$S == "s1"  & pp$R == "r2", FALSE,
          ifelse(pp$S == "s2" & pp$R == "r1",  FALSE, NA))))
  
  dat0 <- dat
  dat0$reps <- NA
  dat0$type <- "Data"
  pp$reps <- factor(pp$reps)
  pp$type <- "Simulation"
  head(dat0)
  head(pp)
  tmp0 <- rbind(dat0, pp)
  
  dplyr::tbl_df(tmp0)
  p1 <- ggplot(tmp0, aes(RT, color = reps, size = type)) +
    geom_freqpoly(binwidth = .05) +
    scale_size_manual(values = c(1, .3)) +
    scale_color_grey(na.value = "black") +
    theme(legend.position = "none") +
    facet_grid(S ~ C)
  print(p1)
  

  p1 <- ggplot(tmp0, aes(RT, color = reps, size = type, shape = C)) +
    geom_freqpoly(binwidth = .05) +
    scale_size_manual(values = c(1, .3)) +
    scale_color_grey(na.value = "black") +
    coord_cartesian(xlim = c(0, 2)) +
    theme(legend.position = "none") +
    facet_grid(S ~ .)
  print(p1)
  
  
  setwd("/media/yslin/KIWI/Documents/DMC-180405/")
  source ("dmc/dmc.R")
  load_model ("LBA","lba_B.R")
  setwd("/media/yslin/KIWI/Documents/ggdmc")
  pp <- post.predict.dmc(fit)
  
  plot.pp.dmc(pp)
  
  pp$RE <- rep(1:100, each = 200)

  bin <- seq(min(dat$RT) - .1, max(dat$RT) + .1, .1)
  
  
  head(pp)
  str(pp)
  
  setwd("/media/yslin/KIWI/Documents/Fixed-Cue-vs-Varied-Cue/")
  source("R/functions/summarise.R")
  cvg <- summarySEwithin(data.frame(pp), 
                          measurevar = 'RT', 
                          withinvar   =c('S','C'), 
                          na.rm = FALSE, conf.interval = .95, .drop=TRUE)
  dvg <- summarySEwithin(data.frame(dat), 
                         measurevar = 'RT', 
                         withinvar   =c('S','C'), 
                         na.rm = FALSE, conf.interval = .95, .drop=TRUE)
  
cvg
dvg
head(dat)
head(pp)

  p <- ggplot(cvg, aes(x = S, y = RT, group = C)) +
    geom_ribbon(aes(ymin = RT-se, ymax = RT+se, 
                    fill = C, colour = C), alpha=.2 )  +
    geom_line(data = dvg, aes(colour = C, linetype = C), size=1.5)  
    # geom_point(aes(shape = C, colour = C), size=4) 
  p
  
  

    geom_text(aes(x, y, label = lab), size=5.5, data = textdf) +
    theme_bw()
  
  p1 <- ggplot(pp, aes(RT, color = C)) +
    geom_freqpoly(binwidth = .05) +
    coord_cartesian(xlim = c(0, 5)) +
    facet_grid(S ~ .)
  print(p1)
  
  pp <- post.predict.dmc(fit)
  str(pp)
  plot.pp.dmc(pp)
  
  
  # And cdf
  plot.pp.dmc(pp,"cdf")

})



