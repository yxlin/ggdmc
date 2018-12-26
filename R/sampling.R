### Sampling  ----------------------------------------------------------------

##' @export
Startlogit <- function(nmc, data = NULL, start = NULL, prior = NULL,
                       thin = 1, nchain = NULL, rp = .001) {
  model  <- CheckDMI(data, prior, NULL, nchain)
  pnames <- GetPNames(model)
  npar   <- length(pnames)
  if (is.null(nchain)) nchain <- 3*npar

  datnames <- names(data)
  Xnames   <- datnames[!datnames %in% c("R", "N", "Y")]
  Ynames   <- c("Y", "N")

  X_ <- data.matrix(data[, Xnames])

  attr(data, "X") <- cbind(X_, X_[,1]*X_[,2])
  attr(data, "Y") <- data.matrix(data[, Ynames])

  out <- init_logit(nmc, data, start, prior, rp, thin, nchain)
  return(out)
}

##' @export
Startlogits <- function(nmc, data = NULL, start = NULL, prior = NULL,
                        thin = 1, nchain = NULL, rp = .001) {
  if (is.null(prior) || is.null(start)) stop("Neither p nor start priors supplied")
  npar <- length(prior)
  if (is.null(nchain)) nchain <- 3*npar

  model1 <- attr(data[[1]], "model")
  pnames <- GetPNames(model1)
  type <- attr(model1, "type")
  out <- vector("list", length = length(data))

  for(i in 1:length(data)) {
    out[[i]] <- Startlogit(nmc, data[[i]], start, prior, thin, nchain, rp)
    dimnames(out[[i]]$theta) <- list(NULL, pnames, NULL)
    class(out[[i]]) <- c("model", "list")
  }
  cat("\n")
  return(out)


}

##' @export
Starthlogit <- function(nmc, data = NULL, start = NULL, prior = NULL,
                         thin = 1, nchain = NULL, rp = .001) {

  nchain <- CheckHyperDMI(data, nchain)  ## If nchain=NULL, change it to default
  if (is.null(prior)) stop("No prior")   ## Check priors
  if (is.null(start)) stop("No start prior")
  model1 <- attr(data[[1]], "model")
  pnames <- GetPNames(model1)
  ispok <- pnames %in% names(prior[[1]])
  islok <- pnames %in% names(prior[[2]])
  issok <- pnames %in% names(prior[[3]])
  if (!all(ispok)) stop("p.prior incompatible with model")
  if (!all(islok)) stop("location prior incompatible with model")
  if (!all(issok)) stop("scale prior incompatible with model")
  type <- attr(model1, "type")

  out <- Startlogits(nmc, data, start[[1]], prior[[1]], thin, nchain, rp)

  names(out) <- names(data)
  class(out) <- c("model", "list")
  return(out)
}


# nmc <- 5e2
# data <- dmi
# start <- start
# prior <- p.prior
# thin <- 1
# nchain <- NULL
# rp <- .001

##' @export
Start_glm <- function(nmc, data = NULL, start = NULL, prior = NULL,
                     thin = 1, nchain = NULL, rp = .001) {

  model  <- CheckDMI(data, prior, NULL, nchain)
  pnames <- GetPNames(model)
  npar   <- length(pnames)
  if (is.null(nchain)) nchain <- 3*npar

  datnames <- colnames(data)
  Xnames   <- datnames[!datnames %in% c("R", "S", "Y")]
  Ynames   <- c("Y")
  attr(data, "X") <- data.matrix(data[, Xnames])
  attr(data, "Y") <- data.matrix(data[, Ynames])
  out <- ggdmc:::init_new_glm(nmc, data, start, prior, rp, thin, nchain)
  return(out)

}

# nmc <- 100
# data <- dmi
# prior <- p.prior
##' @export
Start_glms <- function(nmc, data = NULL, start = NULL, prior = NULL,
                      thin = 1, nchain = NULL, rp = .001) {
  if (is.null(prior) || is.null(start)) stop("Neither p nor start priors supplied")
  npar <- length(prior)
  if (is.null(nchain)) nchain <- 3*npar
  if (is.null(data[[1]]$X)) { stop("No regressors"); }

  model1 <- attr(data[[1]], "model")
  pnames <- GetPNames(model1)
  type <- attr(model1, "type")
  out <- vector("list", length = length(data))

  for(i in 1:length(data)) {
    out[[i]] <- Start_glm(nmc, data[[i]], start, prior, thin, nchain, rp)
    dimnames(out[[i]]$theta) <- list(NULL, pnames, NULL)
    class(out[[i]]) <- c("model", "list")
  }
  cat("\n")
  return(out)
}

# nmc <- 10
# data <- dmi
# start <- start
# prior <- prior
# thin <- 1
# nchain <- 9
# rp <- .001
##' @export
Start_hglm <- function(nmc, data = NULL, start = NULL, prior = NULL,
                      thin = 1, nchain = NULL, rp = .001) {

    nchain <- CheckHyperDMI(data, nchain)  ## If nchain=NULL, change it to default
    if (is.null(prior)) stop("No prior")   ## Check priors
    if (is.null(start)) stop("No start prior")

    model1 <- attr(data[[1]], "model")
    pnames <- GetPNames(model1)

    pprior <- prior[[1]]
    lprior <- prior[[2]]
    sprior <- prior[[3]]
    pstart <- start[[1]]
    lstart <- start[[2]]
    sstart <- start[[3]]
    ispok <- pnames %in% names(pprior)
    islok <- pnames %in% names(lprior)
    issok <- pnames %in% names(sprior)
    if (!all(ispok)) stop("p.prior incompatible with model")
    if (!all(islok)) stop("location prior incompatible with model")
    if (!all(issok)) stop("scale prior incompatible with model")

    npar <- length(pnames)
    nsub <- length(data)

    out <- Start_glms(nmc, data, pstart, pprior, thin, nchain, rp)
    theta <- GetTheta0(out)

    lower <- upper <- dist <- lg <- numeric(npar)
      for(i in 1:npar) {
        v <- pprior[[i]]
        dist[i] <- attr(v, "dist")
        lower[i] <- v$lower
        upper[i] <- v$upper
        lg[i] <- v$lg
      }

      ## containers
      array_size <- c(nchain, npar, nmc)
      location <- array(NA, dim = array_size)
      scale <- array(NA, dim = array_size)
      hlp <- matrix(-Inf, nrow = nmc, ncol = nchain)
      hll <- matrix(-Inf, nrow = nmc, ncol = nchain)

      # i <- 1
    ## fill in containers
    for(i in 1:nchain) {
        lvec <- rprior(lstart)
        svec <- rprior(sstart)
        llik <- sumlogpriorNV(lvec, lprior)
        slik <- sumlogpriorNV(svec, sprior)
        hlp[1, i] <- llik + slik
        location[i,,1] <- lvec
        scale[i,,1] <- svec
        tmp <- 0
        for(j in 1:nsub) {
            tmp <- tmp + ggdmc:::sumlogprior(theta[j,,i], dist, lvec, svec, lower, upper, lg)
        }
        hll[1, i] <- tmp
    }

    found_inf <- any(is.infinite(hll[1,]))
    if (found_inf) { cat("\nFound infinite in hyper log likelihoods") }

    phi <- list(location, scale)
    ppprior <- list(lprior, sprior)
    names(ppprior) <- c("location", "scale")
    hyper <- list(phi, hlp, hll, ppprior, 1, npar, pnames, rp, nmc, thin, nchain)
    names(hyper) <- c("phi", "h_summed_log_prior", "h_log_likelihoods", "pp.prior",
                      "start", "n.pars", "p.names", "rp", "nmc", "thin", "n.chains")
    attr(out, "hyper") <- hyper
    class(out) <- c("model", "list")
    return(out)
}

##' Initialise New Samples
##'
##' Set points radnomly based on supplied / default prior distributions,
##' either from \code{start}, or \code{prior} to generate over-dispersed
##' initial parameter values.
##'
##' @param nmc numbers of Monte Carlo samples.
##' @param data a data model instance
##' @param prior parameter prior distributions from \code{BuildPrior}.
##' @param ppprior hyper parameter prior distributions from \code{BuildPrior}.
##' This must be a set of location and scale hyper prior distributions.
##' @param thin thinning length.
##' @param nchain numbers of Markov chains. Default is 3 times the numbers of
##' model parameters.
##' @param rp DE-MCMC tuning parameter to generate random noise either from
##' uniform or Gaussian distribution.
##' @param samples a collection fo posterior samples.
##' @param add add more samples onto an existing samples
##' @examples
##' m1 <- BuildModel(
##'     p.map     = list(a = "1", v = "F", z = "1", d = "1", sz = "1", sv = "1",
##'                 t0 = "1", st0 = "1"),
##'     constants = c(st0 = 0, d = 0),
##'     match.map = list(M = list(s1 = "r1", s2 = "r2")),
##'     factors   = list(S = c("s1","s2"), F = c("f1", "f2")),
##'     responses = c("r1","r2"),
##'     type      = "rd")
##'
##' ## m1 is "model" class
##' class(m1)
##' ## [1] "model"
##'
##' pVec <- c(a = 1, v.f1 = 1, v.f2 = 1.5, z = .5, sz = .25, sv = .2, t0 = .15)
##' dat  <- simulate(m1, nsim = 1e2, ps = pVec)
##' str(dat)
##' ## 'data.frame':	400 obs. of  4 variables:
##' ## $ S : Factor w/ 2 levels "s1","s2": 1 1 1 1 1 1 1 1 1 1 ...
##' ## $ F : Factor w/ 2 levels "f1","f2": 1 1 1 1 1 1 1 1 1 1 ...
##' ## $ R : Factor w/ 2 levels "r1","r2": 1 1 1 2 1 1 1 1 2 1 ...
##' ## $ RT: num  0.26 0.255 0.572 0.25 0.518 ...
##'
##' dmi1 <- BuildDMI(dat, m1)
##' npar <- length(GetPNames(m1))
##'
##' p.prior <- BuildPrior(
##'    dists = rep("tnorm", npar),
##'    p1    = c(a=2,  v.f1=2.5, v.f2=1.25, z=.5, sz=.3, sv=1,  t0=.3),
##'    p2    = c(a=.5, v.f1=.5,  v.f2=.35,  z=.1, sz=.1, sv=.3, t0=.05),
##'    lower = c(0,-5, -5, 0, 0, 0, 0),
##'    upper = c(5, 7,  7, 2, 2, 2, 2))
##'
##' ## Set up a new DMC sample with 16 iteration. The default thin is 1
##' sam0 <- StartNewsamples(nmc = 16, data = dmi1, prior = p.prior)
##' sam0$nmc
##' ## [1] 16
##'
##' sam1 <- RestartSamples(16, sam0)
##' sam1$nmc
##' ## [1] 16
##' sam1 <- RestartSamples(16, sam1, add = TRUE)
##' sam1$nmc
##' ## [1] 32
##'
##' #####################28
##' ## Hierarchical      ##
##' #####################28
##' model <- BuildModel(
##'         p.map     = list(A = "1", B = "R", t0 = "1", mean_v = c("D", "M"),
##'           sd_v = "M", st0 = "1"),
##'         match.map = list(M = list(s1 = 1, s2 = 2)),
##'         factors   = list(S = c("s1", "s2"), D = c("d1", "d2")),
##'         constants = c(sd_v.false = 1, st0 = 0),
##'         responses = c("r1", "r2"),
##'         type      = "norm")
##'
##' ## Population distribution, rate effect on F
##' pop.mean <- c(A=.4, B.r1=.6, B.r2=.8, t0=.3,
##'   mean_v.d1.true  = 1.5,
##'   mean_v.d2.true  = 1.0,
##'   mean_v.d1.false = 0,
##'   mean_v.d2.false = 0,  sd_v.true = .25)
##' pop.scale <-c(A=.1, B.r1=.1, B.r2=.1, t0=.05,
##'   mean_v.d1.true  =.2,
##'   mean_v.d2.true  =.2,
##'   mean_v.d1.false =.2,
##'   mean_v.d2.false =.2,  sd_v.true = .1)
##'
##' pop.prior <- BuildPrior(
##'   dists = rep("tnorm", 9),
##'   p1 = pop.mean,
##'   p2 = pop.scale,
##'   lower = c(0,0,0,   .1, NA,NA,NA,NA, 0),
##'   upper = c(NA,NA,NA, 1, NA,NA,NA,NA, NA))
##'
##' dat <- simulate(model, nsub = 6, nsim = 30, prior = pop.prior)
##' dmi <- BuildDMI(dat, model)
##' p.prior <- BuildPrior(
##'   dists = rep("tnorm", 9),
##'   p1   = pop.mean,
##'   p2   = pop.scale*5,
##'   lower=c(0,0,0,   .1, NA,NA,NA,NA, 0),
##'   upper=c(NA,NA,NA,NA, NA,NA,NA,NA, NA)
##' )
##'
##' mu.prior <- BuildPrior(
##'   dists = rep("tnorm",  9),
##'   p1    = pop.mean,
##'   p2    = c(1,   1,  1,  1,   2,  2,  2, 2,  1),
##'   lower = c(0,   0,  0, .01, NA, NA, NA, NA, 0),
##'   upper = c(NA, NA, NA,  NA, NA, NA, NA, NA, NA)
##' )
##'
##' sigma.prior <- BuildPrior(
##'   dists = rep("beta", length(p.prior)),
##'   p1    = c(A = 1, B.r1 = 1, B.r2 = 1, t0 = 1, mean_v.d1.true = 1,
##'     mean_v.d2.true = 1, mean_v.d1.false = 1, mean_v.d2.false = 1,
##'     sd_v.true = 1),
##'   p2    = rep(1, 9))
##'
##' pp.prior <- list(mu.prior, sigma.prior)
##' hsam <- StartNewHypersamples(32, dmi, pop.prior, pp.prior)
##'
##' @export
StartNewsamples <- function(nmc, data = NULL, prior = NULL, thin = 1,
  nchain = NULL, rp = .001) {

  ## theta1 == NULL
  model <- CheckDMI(data, prior, NULL, nchain)
  npar  <- length(GetPNames(model))
  if (is.null(nchain)) nchain <- 3*npar

  if (is.null(data)) stop("Please supply data-model instance.")
  if (is.null(prior)) stop("Please supply prior distributions")
  if (length(prior) != npar) stop("Different numbers of parameters in p.vector and p.prior ")

  ## Temporary measures
  if (is.null(attr(data, "n.pda"))) attr(data, "n.pda") <- 2^14
  if (is.null(attr(data, "bw")))    attr(data, "bw")    <- .01
  if (is.null(attr(data, "debug"))) attr(data, "debug") <- 0
  if (is.null(attr(data, "gpuid"))) attr(data, "gpuid") <- 0
  ncore <- 1
  debug <- FALSE

  type <- attr(model, "type")

  if (type == "glm") {
    if (is.null(data$X)) {
      message("Add dummy X column");
      data$X <- rep(1, nrow(data))
    }
    out <- init_new_glm(nmc, prior, data, rp, thin, nchain)

  } else {
    out <- init_new(nmc, prior, data, rp, thin, nchain, ncore, debug)
  }

  class(out) <- c("model", "list")
  cat("\n")
  return(out)
}


##' @rdname StartNewsamples
##' @export
RestartSamples <- function(nmc, samples = NULL, thin = NULL, rp = .001,
  add = FALSE) {

  hyper <- attr(samples, "hyper")
  # model <- CheckSamples(samples, p.prior, NULL)
  # npar  <- length(GetPNames(model))
  if (is.null(samples)) stop("Use StartNewsamples")
  if (!is.null(hyper)) stop("Use RestartHypersamples")
  if (is.null(thin)) {
    thin <- samples$thin
  } else if (thin < 0 ) {
    stop("thin must be a positive integer")
  } else {
    cat("Reset thin to ", thin, "\n")
  }

  if (is.null(attr(samples$data, "n.pda"))) attr(samples$data, "n.pda") <- 2^14
  if (is.null(attr(samples$data, "bw")))    attr(samples$data, "bw")    <- .01
  if (is.null(attr(samples$data, "debug"))) attr(samples$data, "debug") <- 0
  if (is.null(attr(samples$data, "gpuid"))) attr(samples$data, "gpuid") <- 0

  if (add) {
    out <- init_add(nmc, samples, rp, thin) # add onto existing one
  } else {
    out <- init_old(nmc, samples, rp, thin) # start afresh
  }

  class(out) <- c("model", "list")
  cat("\n")
  return(out)
}

##' @rdname StartNewsamples
##' @export
StartManynewsamples <- function(nmc, data = NULL, prior = NULL, thin = 1,
                                nchain = NULL, rp = .001) {
  npar <- length(prior)
  if (is.null(nchain)) nchain <- 3*npar
  if (is.null(prior)) stop("No p.prior no new samples")

  if (is.null(attr(data[[1]], "bw"))) {
    message("No GPU attributes. Default bw = .001, using GPU 0.")
    for(i in 1:length(data)) {
      attr(data[[i]], "n.pda") <- 1e4
      attr(data[[i]], "bw") <- .001
      attr(data[[i]], "gpuid") <- 0
    }
  }

  model1 <- attr(data[[1]], "model")
  type <- attr(model1, "type")
  if (type == "glm") {
    if (is.null(data[[1]]$X)) { stop("No regressors"); }
    out <- init_newnonhier_glm(nmc, data, prior, rp, thin, nchain)
  } else {
    out <- init_newnonhier(nmc, data, prior, rp, thin, nchain)
  }

  for(i in 1:length(out)) {
    pnames <- GetPNames(attr(out[[i]]$data, "model"))
    dimnames(out[[i]]$theta) <- list(NULL, pnames, NULL)
    class(out[[i]]) <- c("model", "list")
  }
  cat("\n")
  return(out)
}

##' @rdname StartNewsamples
##' @export
RestartManysamples <- function(nmc, samples = NULL, thin = NULL,
  rp = .001, add = FALSE) {

  if (is.null(samples)) stop("Use StartManynewsamples")
  if (is.null(thin)) thin <- samples[[1]]$thin

  if (add) {
    out <- init_addnonhier(nmc, samples, rp, thin)
  } else {
    out <- init_oldnonhier(nmc, samples, rp, thin)
  }

  for(i in 1:length(out)) {
    pnames <- GetPNames(attr(out[[i]]$data, "model"))
    dimnames(out[[i]]$theta) <- list(NULL, pnames, NULL)
    class(out[[i]]) <- c("model", "list")
  }
  cat("\n")
  return(out)
}

##' @rdname StartNewsamples
##' @export
StartNewHypersamples <- function(nmc, data = NULL, prior = NULL,
  ppprior = NULL, thin = 1, nchain = NULL, rp = .001) {

  nchain <- CheckHyperDMI(data, nchain)    ## If nchain=NULL, change it to default
  if (is.null(prior)) stop("No p.prior")   ## Check priors
  if (is.null(ppprior)) stop("No pp.prior")
  if (!is.list(ppprior)) stop("pp.prior must be a list")
  if (length(ppprior[[1]]) < length(ppprior[[2]]))
    stop("Location priors must have as many or more elements than scale priors")

  if (is.null(attr(data[[1]], "bw"))) {
    message("No GPU attributes. Default bw = .01, using GPU 0.")
    for(i in 1:length(data)) {
      attr(data[[i]], "n.pda") <- 16384
      attr(data[[i]], "bw") <- .01
      attr(data[[i]], "gpuid") <- 0
    }
  }

  model1 <- attr(data[[1]], "model")
  pnames <- GetPNames(model1)
  isppriorok <- pnames %in% names(prior)
  islpriorok <- pnames %in% names(ppprior[[1]])
  isspriorok <- pnames %in% names(ppprior[[2]])
  if (!all(isppriorok)) stop("p.prior incompatible with model")
  if (!all(islpriorok)) stop("location prior incompatible with model")
  if (!all(isppriorok)) stop("scale prior incompatible with model")

  type <- attr(model1, "type")
  if (type == "glm") {
    if (is.null(data[[1]]$X)) {
      stop("No regressor");
    }
    out <- init_newhier_glm(nmc, data, prior, ppprior, rp, thin, nchain)

  } else {
    out <- init_newhier(nmc, data, prior, ppprior, rp, thin, nchain)
  }

  names(out) <- names(data)
  class(out) <- c("model", "list")
  return(out)
}

##' @rdname StartNewsamples
##' @export
RestartHypersamples <- function(nmc, samples = NULL, thin = NULL, rp = .001,
                                add = FALSE) {

  if (is.null(samples)) stop("Use StartNewHypersamples")
  if (is.null(thin)) thin <- samples[[1]]$thin
  model1 <- attr(samples[[1]]$data, "model")
  pnames <- GetPNames(model1)
  type <- attr(model1, "type")

  if (add) {
    if (type == "glm") {
      out <- init_addhier_glm(nmc, samples, rp, thin) # add onto existing one
    } else {
      out <- init_addhier(nmc, samples, rp, thin) # add onto existing one
    }
  } else {
    if (type == "glm") {
      out <- init_oldhier_glm(nmc, samples, rp, thin)
    } else {
      out <- init_oldhier(nmc, samples, rp, thin)
    }
  }
  names(out) <- names(samples)
  class(out) <- c("model", "list")
  return(out)
}

##' @rdname StartNewsamples
##' @export
StartNewhiersamples <- function(nmc, data = NULL, start = NULL,
                                prior = NULL, thin = 1, nchain = NULL,
                                rp = .001)
{
  ppprior <- prior[2:3]

  nchain <- CheckHyperDMI(data, nchain)  ## If nchain=NULL, change it to default
  if (is.null(prior)) stop("No prior")   ## Check priors
  if (is.null(start)) stop("No start prior")
  if (length(ppprior[[1]]) < length(ppprior[[2]])) {
    stop("Location priors must have as many or more elements than scale priors")
  }

  if (is.null(attr(data[[1]], "bw"))) {
    message("No GPU attributes. Default bw = .01, using GPU 0.")
    for(i in 1:length(data)) {
      attr(data[[i]], "n.pda") <- 16384
      attr(data[[i]], "bw") <- .01
      attr(data[[i]], "gpuid") <- 0
    }
  }

  model1 <- attr(data[[1]], "model")
  pnames <- GetPNames(model1)
  isppriorok <- pnames %in% names(prior[[1]])
  islpriorok <- pnames %in% names(ppprior[[1]])
  isspriorok <- pnames %in% names(ppprior[[2]])
  if (!all(isppriorok)) stop("p.prior incompatible with model")
  if (!all(islpriorok)) stop("location prior incompatible with model")
  if (!all(isppriorok)) stop("scale prior incompatible with model")

  type <- attr(model1, "type")
  if (type == "glm") {
    if (is.null(data[[1]]$X)) {
      stop("No regressor");
    }
    out <- init_newhier_start(nmc, data, start, prior, rp, thin, nchain)

  } else {
    stop("Unknonw error")
  }

  names(out) <- names(data)
  class(out) <- c("model", "list")
  return(out)
}


##' Fit a Bayesian Model to a Single Participant
##'
##' Use either DE-MCMC or DGMC sampler to fit Bayesian model to a participant
##'
##' @param samples a initialized sample
##' @param report progress report intervel
##' @param pm probability of migration type 1
##' @param pm0 probability of migration type 2 (Turner et al., 2013)
##' @param qm probability of mutation
##' @param gammamult tuning parameter of the crossover sampler
##' @param ngroup number of independent groups
##' @param force PDA re-calculate interval
##' @param sampler a string indicating to use which sampler
##' @param slice a Boolean switch to debug blocked sampling
##' @return posterior samples
##' @export
run_one <- function(samples, report, pm, pm0, qm, gammamult, ngroup, force,
  sampler, slice) {

  ## message("Run one subject. Enforce ncore = 1.")
  ncore  <- 1
  force  <- MakeForce(samples, force)
  pnames <- GetPNames(attr(samples$data, "model"))
  nchain <- samples$n.chains
  type <- attr(attr(samples$data, "model"), "type")


  if (is.null(attr(samples$data, "n.pda"))) attr(samples$data, "n.pda") <- 2^14
  if (is.null(attr(samples$data, "bw")))    attr(samples$data, "bw")    <- .01
  if (is.null(attr(samples$data, "debug"))) attr(samples$data, "debug") <- 0
  if (is.null(attr(samples$data, "gpuid"))) attr(samples$data, "gpuid") <- 0

  if (sampler == "DGMC") {
    message("DGMC")

    out <- run_dgmc(samples, force, report, pm, qm, gammamult, ncore, ngroup)
  } else if (sampler == "DE-MCMC") {
    message("DE-MCMC")

    if (type == "glm") {
      out <- run_glm(samples, force, report, pm, pm0, gammamult, ncore)
    } else if (type == "logit") {
      out <- run_glm(samples, force, report, pm, pm0, gammamult, ncore)
    } else {
      out <- run_dmc(samples, force, report, pm, pm0, gammamult, ncore, slice)
    }

  } else {
    stop ("Sampler yet implemented")
  }

  dimnames(out$theta) <- list(NULL, pnames, NULL)
  return(out)
}

##' Fit a Bayesian Model to multiple Participants
##'
##' Use either DE-MCMC or DGMC sampler to fit independent Bayesian model to
##' many participants.
##'
##' @param samples a initialized samples list. Each element should contain
##' samples for a participant.
##' @param report progress report intervel
##' @param ncore number of CPU cores
##' @param pm probability of migration
##' @param pm0 probability of migration
##' @param qm probability of mutation
##' @param gammamult turning parameter for crossover sampler
##' @param ngroup number of independent groups
##' @param force PDA re-calculate interval
##' @param sampler a string indicating to use which sampler
##' @param slice a Boolean switch to do slice sampling
##' @importFrom parallel parLapply
##' @importFrom parallel makeCluster
##' @importFrom parallel mclapply
##' @importFrom parallel stopCluster
##' @export
run_many <- function(samples, report, ncore, pm, pm0, qm, gammamult, ngroup,
  force, sampler, slice) {

  force <- MakeForce(samples[[1]], force)
  for(i in 1:length(samples)) {
    if (is.null(attr(samples[[i]]$data, "n.pda"))) {
      attr(samples[[i]]$data, "n.pda") <- 2^14
    }
    if (is.null(attr(samples[[i]]$data, "bw"))) {
      attr(samples[[i]]$data, "bw")    <- .01
    }
    if (is.null(attr(samples[[i]]$data, "debug"))) {
      attr(samples[[i]]$data, "debug") <- 0
    }
    if (is.null(attr(samples[[i]]$data, "gpuid"))) {
      attr(samples[[i]]$data, "gpuid") <- 0
    }
  }

  if (get_os() == "windows" & ncore > 1) {
    cl  <- parallel::makeCluster(ncore)
    if (sampler == "DGMC") {
      message("Run many subjects using DGMC in parallel")

      out <- parallel::parLapply(cl, samples, run_dgmc, force, report, pm, qm,
        gammamult, ncore, ngroup)
    } else if (sampler == "DE-MCMC") {
      message("Run many subjects using DE-MCMC in parallel")

      out <- parallel::parLapply(cl, samples, run_dmc, force, report, pm, pm0,
        gammamult, ncore, slice)
    } else {
      stop("Sampler unknown")
    }
    stopCluster(cl)

  } else if (ncore > 1) {
    if (sampler == "DGMC") {
      message("Run many subjects using DGMC in parallel")

      out <- parallel::mclapply(samples, run_dgmc, force, report, pm, qm,
        gammamult, ncore, ngroup, mc.cores=ncore)
    } else if (sampler == "DE-MCMC") {
      message("Run many subjects using DE-MCMC in parallel")

      out <- parallel::mclapply(samples, run_dmc, force, report, pm, pm0,
        gammamult, ncore, slice, mc.cores = ncore)
    } else {
      stop("Sampler unknown")
    }

  } else {
    if (sampler == "DGMC") {
      message("Run many subject using DGMC with lapply")

      out <- lapply(samples, run_dgmc, force, report, pm, qm, gammamult, ncore,
        ngroup)
    } else if (sampler == "DE-MCMC") {
      message("Run many subject using DE-MCMC lapply")

      out <- lapply(samples, run_dmc, force, report, pm, pm0, gammamult, ncore, slice)
    } else {
      stop ("Sampler unknown")
    }
  }

  for(i in 1:length(out)) {
    dimnames(out[[i]]$theta) <- list(NULL, out[[i]]$p.names, NULL)
  }
  return(out)
}

run_many_glm <- function(samples, report, ncore, pm, pm0, qm, gammamult, ngroup,
                     force, sampler) {

  force <- MakeForce(samples[[1]], force)
  for(i in 1:length(samples)) {
    if (is.null(attr(samples[[i]]$data, "n.pda"))) {
      attr(samples[[i]]$data, "n.pda") <- 2^14
    }
    if (is.null(attr(samples[[i]]$data, "bw"))) {
      attr(samples[[i]]$data, "bw")    <- .01
    }
    if (is.null(attr(samples[[i]]$data, "debug"))) {
      attr(samples[[i]]$data, "debug") <- 0
    }
    if (is.null(attr(samples[[i]]$data, "gpuid"))) {
      attr(samples[[i]]$data, "gpuid") <- 0
    }
  }

  if (get_os() == "windows" & ncore > 1) {
    cl  <- parallel::makeCluster(ncore)
    if (sampler == "DGMC") {
      message("Run many subjects using DGMC in parallel")

      out <- parallel::parLapply(cl, samples, run_dgmc, force, report, pm, qm,
                                 gammamult, ncore, ngroup)
    } else if (sampler == "DE-MCMC") {
      message("Run many subjects using DE-MCMC in parallel")

      out <- parallel::parLapply(cl, samples, run_glm, force, report, pm, pm0,
                                 gammamult, ncore)
    } else {
      stop("Sampler unknown")
    }
    stopCluster(cl)

  } else if (ncore > 1) {
    if (sampler == "DGMC") {
      message("Run many subjects using DGMC in parallel")

      out <- parallel::mclapply(samples, run_dgmc, force, report, pm, qm,
                                gammamult, ncore, ngroup, mc.cores=ncore)
    } else if (sampler == "DE-MCMC") {
      message("Run many subjects using DE-MCMC in parallel")

      out <- parallel::mclapply(samples, run_glm, force, report, pm, pm0,
                                gammamult, ncore, mc.cores = ncore)
    } else {
      stop("Sampler unknown")
    }

  } else {
    if (sampler == "DGMC") {
      message("Run many subject using DGMC with lapply")

      out <- lapply(samples, run_dgmc, force, report, pm, qm, gammamult, ncore,
                    ngroup)
    } else if (sampler == "DE-MCMC") {
      message("Run many subject using DE-MCMC lapply")

      out <- lapply(samples, run_glm, force, report, pm, pm0, gammamult, ncore)
    } else {
      stop ("Sampler unknown")
    }
  }

  for(i in 1:length(out)) {
    dimnames(out[[i]]$theta) <- list(NULL, out[[i]]$p.names, NULL)
  }
  return(out)
}

run_many_logit <- function(samples, report, ncore, pm, pm0, gammamult, force) {
  force <- MakeForce(samples[[1]], force)
  if (ncore > 1) {
    message("Run many subjects using DE-MCMC in parallel (logit)")
    out <- parallel::mclapply(samples, run_glm, force, report, pm, pm0,
                              gammamult, ncore, mc.cores = ncore)
  } else {
    message("Run many subjects using DE-MCMC with lapply (logit)")
    out <- lapply(samples, run_glm, force, report, pm, pm0,
                  gammamult, ncore)
  }
  return(out)
}
##' Run model fits
##'
##' This function fit a hierarchical or a fixed-effect model, using Bayeisan
##' sampling.  We use pMCMC, with a suite of DE-MCMC, DGMC, and simply,
##' crossover (i.e., DE-MC), mutation, or migration operators. Note that
##' the latter two operators essentially are random-walk Metroplolis, so they
##' will be very inefficient, if been applied alone, even with our fast C++
##' implementation.
##'
##' @param samples a sample list generated by calling DMC's samples.dmc.
##' @param report how many iterations to return a report
##' @param ncore parallel core for run_many
##' @param pm probability of migration
##' @param pm0 probability of old migration
##' @param qm probability of mutation
##' @param hpm probability of migration at the hyper level
##' @param hpm0 probability of old migration at the hyper level
##' @param hqm probability of mutation at the hyper level
##' @param gammamult a tuning parameter, affecting the size of jump
##' @param ngroup number of distributed groups
##' @param force set force to FALSE for turning off recalculation of PDA.
##' Set it  as an integer between 1 and 10, forcing to re-calculate new
##' likelihood, every e.g., 1, 2, 3 step.
##' @param sampler which sampler to run MCMC, "DE-MCMC" or "DGMC"
##' @param slice use for debugging blocked sampling
##'
##' @export
run <- function(samples, report = 1e2, ncore = 1, pm = 0, pm0 = 0, qm = 0,
                hpm = 0, hpm0 = 0, hqm = 0, gammamult = 2.38, ngroup = 5,
                force = FALSE,
  sampler = "DE-MCMC", slice = FALSE) {

  hyper <- attr(samples, "hyper")

  if (!is.null(hyper)) {   ## hierarchical model
    if (is.null(attr(samples[[1]]$data, "bw"))) {
      message("No GPU attributes. Default bandwidth = .01, using GPU 0.")
      for(i in 1:length(samples)) {
        attr(samples[[i]]$data, "bw") <- .01
        attr(samples[[i]]$data, "gpuid") <- 0
      }
    }

    if (sampler == "DE-MCMC") {
      type <- attr(attr(samples[[1]]$data, "model"), "type")

      if (type == "norm" ) checklba(samples[[1]])
      message("DE-MCMC; hierarchical modeling")

      if (type == "glm") {
        message("Run hglm")
        out <- run_hglm(samples, report, pm, pm0, hpm, hpm0, gammamult)
      } else {
        out <- run_hyper_dmc(samples, report, pm, pm0, hpm, hpm0, gammamult,
                             ncore, FALSE)
      }


    } else if (sampler == "DGMC") {
      if (attr(attr(samples[[1]]$data, "model"), "type") == "norm" ) checklba(samples[[1]])

      message("DGMC; hierarchical modeling")
      out <- run_hyper_dgmc(samples, report, pm, hpm, qm, hqm, gammamult,
        ngroup, ncore)
    } else {
      out <- NULL
      message("Unknown sampler?")
    }

    # hyper <- attr(out, "hyper")
    # notna <- lapply(hyper$pp.prior, function(x) {
    #          sapply(x, function(y) { !is.na(attr(y, "dist")) } )
    #       }
    # )

    pnames <- hyper$p.names
    # location_names <- pnames[notna$location]
    # scale_names <- pnames[notna$scale]

    phi1_tmp <- attr(out, "hyper")$phi[[1]]
    phi2_tmp <- attr(out, "hyper")$phi[[2]]
    pnames <- GetPNames(attr(samples[[1]]$data, "model"))
    dimnames(phi1_tmp) <- list(NULL, pnames, NULL)
    dimnames(phi2_tmp) <- list(NULL, pnames, NULL)
    attr(out, "hyper")$phi[[1]] <- phi1_tmp
    attr(out, "hyper")$phi[[2]] <- phi2_tmp

    for(i in 1:length(out)) {
      pnames <- GetPNames(attr(out[[i]]$data, "model"))
      dimnames(out[[i]]$theta) <- list(NULL, pnames, NULL)
      class(out[[i]]) <- c("model", "list")
    }

  } else if (any(names(samples) == "theta")) { ## One subject

    if (attr(attr(samples$data, "model"), "type") == "norm" ) checklba(samples)

    out <- run_one(samples, report, pm, pm0, qm, gammamult, ngroup, force,
      sampler, slice)
    pnames <- GetPNames(attr(out$data, "model"))
    dimnames(out$theta) <- list(NULL, pnames, NULL)

  } else {  ## Multiple independent subjects

    type <- attr(attr(samples[[1]]$data, "model"), "type")
    if (type == "norm" ) checklba(samples[[1]])
    if (type == "glm" ) {
      out <- run_many_glm(samples, report, ncore, pm, pm0, qm, gammamult, ngroup,
                      force, sampler)
    } else if (type == "logit") {
      out <- run_many_logit(samples, report, ncore, pm, pm0, gammamult, force)

    } else {
      out <- run_many(samples, report, ncore, pm, pm0, qm, gammamult, ngroup,
                      force, sampler, slice)
    }

    for(i in 1:length(out)) {
      pnames <- GetPNames(attr(out[[i]]$data, "model"))
      dimnames(out[[i]]$theta) <- list(NULL, pnames, NULL)
      class(out[[i]]) <- c("model", "list")
    }
  }

  class(out) <- c("model", "list")
  cat("\n")
  return(out)
}


##' @rdname run
##' @export
CheckConverged <- function(samples) {
  stuck <- isstuck(samples, verbose = FALSE, cut = 10)
  flat  <- isflat(samples, p1 = 1/3, p2 = 1/3,
    cut_location = 0.25, cut_scale = Inf, verbose = FALSE)
  mix  <- ismixed(samples, cut = 1.05, verbose = FALSE)
  size <- iseffective(samples, minN = 512, nfun = "mean", FALSE)
  isstuck <- TRUE
  if (stuck == 0) isstuck <- FALSE

  out <- c(isstuck, flat, mix, size)
  names(out) <- c("Stuck", "Flat", "Mix", "ES")
  return(out)
}

