context("Post-predictive")

test_that("Post-predictive", {
  rm(list = ls())
  predict_one <- function(object, npost = 100, rand = TRUE, factors = NA,
                          xlim = NA, seed = NULL)
  {
    model <- attributes(object$data)$model
    facs <- names(attr(model, "factors"))
    class(object$data) <- c("data.frame", "list")

    if (!is.null(factors))
    {
      if (any(is.na(factors))) factors <- facs
      if (!all(factors %in% facs))
        stop(paste("Factors argument must contain one or more of:",
                   paste(facs, collapse=",")))
    }

    resp <- names(attr(model, "responses"))
    ns   <- table(object$data[,facs], dnn = facs)
    npar   <- object$n.pars
    nchain <- object$n.chains
    nmc    <- object$nmc
    ntsample <- nchain * nmc
    pnames   <- object$p.names

    thetas <- matrix(aperm(object$theta, c(3,2,1)), ncol = npar)

    colnames(thetas) <- pnames

    if (is.na(npost)) {
      use <- 1:ntsample
    } else {
      if (rand) {
        use <- sample(1:ntsample, npost, replace = F)
      } else {
        use <- round(seq(1, ntsample, length.out = npost))
      }
    }

    npost  <- length(use)
    posts   <- thetas[use, ]
    nttrial <- sum(ns) ## number of total trials

    v <- lapply(1:npost, function(i) {
      ggdmc:::simulate_one(model, n = ns, ps = posts[i,], seed = seed)
    })
    out <- data.table::rbindlist(v)

    reps <- rep(1:npost, each = nttrial)
    out <- cbind(reps, out)

    if (!any(is.na(xlim)))
    {
      out <- out[RT > xlim[1] & RT < xlim[2]]
    }

    attr(out, "data") <- object$data
    return(out)
  }

  model <- BuildModel(
    p.map     = list(a = "1", v="1", z="1", d="1", sz="1", sv="1", t0="1", st0="1"),
    match.map = list(M = list(s1 = "r1", s2 = "r2")),
    factors   = list(S = c("s1", "s2")),
    responses = c("r1","r2"),
    constants = c(st0 = 0, d = 0, sv = 0, sz = 0),
    type      = "rd")

  npar <- length(GetPNames(model))
  p.vector <- c(a=1, v=1.5, z=0.5, t0=.15)
  dat <- simulate(model, nsim = 1e3, ps = p.vector)
  dmi <- BuildDMI(dat, model)

  p.prior <- BuildPrior(
    dists = rep("tnorm", npar),
    p1=c(a=1, v=0, z=1, t0=1),
    p2=c(a=1, v=2, z=1, t0=1),
    lower = c(0, -5, rep(0, 2)),
    upper = rep(NA, npar))

  ## Fit model ----
  fit <- run(StartNewsamples(dmi, p.prior))
  rhat <- gelman(fit, verbose = TRUE)

  pdf("Postpredictive.pdf")
  p0 <- plot(fit)
  p0 <- plot(fit, pll=F, den=T)

  ## Analysis -----------
  est <- summary(fit, recovery = TRUE, ps = p.vector, verbose = TRUE)
  pp <- predict_one(fit, xlim = c(0, 5))

  dat$C <- ifelse(dat$S == "s1" & dat$R == "r1", TRUE,
           ifelse(dat$S == "s2" & dat$R == "r2", TRUE,
           ifelse(dat$S == "s1" & dat$R == "r2", FALSE,
           ifelse(dat$S == "s2" & dat$R == "r1", FALSE, NA))))

  pp$C <- ifelse(pp$S == "s1" & pp$R == "r1", TRUE,
          ifelse(pp$S == "s2" & pp$R == "r2", TRUE,
          ifelse(pp$S == "s1" & pp$R == "r2", FALSE,
          ifelse(pp$S == "s2" & pp$R == "r1", FALSE, NA))))

  dat0 <- dat
  dat0$reps <- NA
  dat0$type <- "Data"
  pp$reps <- factor(pp$reps)
  pp$type <- "Simulation"
  tmp0 <- rbind(dat0, pp)


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
})



