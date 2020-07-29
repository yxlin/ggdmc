### Model ------------------------------------------------------------
##' An S4 class of the process model.
##'
##' The class is to represent a process model, e.g., a DDM, a LBA model, a PM
##' model, or a CDDM.
##'
##' @slot model A 3-D model array. Dimension one stores the combinations of the
##' factor levels and response types (when discrete), dimension two stores
##' parameters, and dimension three stores response types.
##' @slot all.par all parameters
##' @slot p.vector parameter vector, excluding constant parameters
##' @slot par.names parameter names / labels
##' @slot type model type
##' @slot factors a list of factors and their levels
##' @slot responses response types
##' @slot constants constant parameters
##' @slot posdrift a Boolean switch indicating whether drift rates must be
##' positive
##' @slot n1.order node 1 ordering. This is only for the LBA model
##' @slot match.cell an indicator matrix storing whether a particular trial
##' matches a cell
##' @slot match.map a mapping mechanism for calculating whether a trial matches
##' a positive boundary / accumulator or a negative boundary / accumulator.
##' @slot dimnames dimension names of the model array
##' @slot pnames parameter names
##' @slot npar number of parameters
##' @export
setClass("model", slots = c(
             model = "array",
             all.par    = "numeric",
             p.vector   = "numeric",
             par.names  = "character",
             type       = "character",
             factors    = "list",
             responses  = "character",
             constants  = "numeric",
             posdrift   = "logical",
             n1.order   = "matrix",
             match.cell = "logical",
             match.map  = "list",
             dimnames   = "list",
             pnames     = "character",
             npar       = "numeric"))

### Data-model Instance  ------------------------------------------------------
##' An S4 class of the Data-model Instance
##'
##' The class is to represent a data-model instance, which joins a model object
##' with a data frame. The process of BuildDMI also generates cell.index and
##' cell.empty.
##'
##' @slot data A data frame storing the would-be fit data set
##' @slot model A 3-D model array.  Dimension one stores the combinations
##' of the factor levels and response types, dimension two stores parameters,
##' and dimension three stores response types.
##' @slot cell.index A ncell-element list. Each element represents one cell.
##' Each element stores \code{nobs} Boolean indicators, showing whether a
##' particular observation belongs to this cell.
##' @slot cell.empty A ncell-element logical vector, indicating whether this
##' cell has no observation.
##' @export
setClass("dmi", slots = c(
           data       = "data.frame",
           model      = "model",
           cell.index = "list",
           cell.empty = "logical"))

### Prior ------------------------------------------------------------
##' An S4 class to represent an object storing prior distributions
##'
##' @slot npar the number of parameters
##' @slot pnames the names of parameters
##' @slot priors a list storing the location parameter, scale parameter, upper
##' bound, lower bound, log indicator (0=FALSE, 1=TRUE), distribution type and
##' transform information.
##' @export
setClass("prior", slots =  c(
  npar   = "numeric",
  pnames = "character",
  priors = "list"))

##' The Parameter Names in a Prior Object
##'
##' Extract parameter names from a prior object. This function extends the
##' \code{names} funciton in the \code{base} package.
##'
##' @param x a prior object.
##'
##' @return a string vector
##' @export
##' @docType methods
##' @rdname names-methods
setMethod("names", "prior", function (x) { slot(x, "pnames") })

##' Generate Random Numbers
##'
##' Random number generation based on a prior object
##'
##' @param x a prior object.
##' @param n number of observations
##' @param ... Additional argument passing via dot dot dot.
##' @examples
##' p.prior <- BuildPrior(
##'  dists = c("tnorm", "tnorm", "beta", "tnorm", "beta", "beta"),
##'  p1    = c(a = 1, v = 0, z = 1, sz = 1, sv = 1, t0 = 1),
##'  p2    = c(a = 1, v = 2, z = 1, sz = 1, sv = 1, t0 = 1),
##'  lower = c(0,-5, NA, NA, 0, NA),
##'  upper = c(2, 5, NA, NA, 2, NA))
##'
##' rprior(p.prior, 9)
##' ##               a           v         z         sz        sv         t0
##' ## [1,] 0.97413686  0.78446178 0.9975199 -0.5264946 0.5364492 0.55415052
##' ## [2,] 0.72870190  0.97151662 0.8516604  1.6008591 0.3399731 0.96520848
##' ## [3,] 1.63153685  1.96586939 0.9260939  0.7041254 0.4138329 0.78367440
##' ## [4,] 1.55866180  1.43657110 0.6152371  0.1290078 0.2957604 0.23027759
##' ## [5,] 1.32520281 -0.07328408 0.2051155  2.4040387 0.9663111 0.06127237
##' ## [6,] 0.49628528 -0.19374770 0.5142829  2.1452972 0.4335482 0.38410626
##' ## [7,] 0.03655549  0.77223432 0.1739831  1.4431507 0.6257398 0.63228368
##' ## [8,] 0.71197612 -1.15798082 0.8265523  0.3813370 0.4465184 0.23955415
##' ## [9,] 0.38049166  3.32132034 0.9888108  0.9684292 0.8437480 0.13502154
##'
##' @docType methods
##' @rdname rprior-methods
##' @export
setGeneric("rprior", function(x, ... ) {
  warning("Class ", class(x), " not defined for rprior")
  return(NULL)
})

##' @rdname rprior-methods
setMethod("rprior", "prior", function (x, n = 1) { rprior_mat(x, n) })

### Posterior ---------------------------------------------------------
##' An S4 class to represent an object storing posterior samples at the
##' participant level. Posterior samples storing both the participant and the
##' hyper lever are represented by an S4 class hyper
##'
##' @slot theta posterior samples for one-participant fit.
##' @slot summed_log_prior summed log prior likelihoods.
##' @slot log_likelihoods log likelihoods
##' @slot dmi a S4 object of data model instance
##' @slot prior a S4 prior object
##' @slot start the index of starting sample
##' @slot npar number of parameters
##' @slot pnames parameter names
##' @slot nmc number of Monte Carlo samples
##' @slot thin thinning length
##' @slot nchain number of Markov chains
##' @seealso \code{\link{hyper-class}}
##' @export
setClass("posterior",
         slot= c(
           theta = "array",
           summed_log_prior = "matrix",
           log_likelihoods  = "matrix",
           dmi    = "dmi",
           prior  = "prior",
           start  = "numeric",
           npar   = "numeric",
           pnames = "character",
           nmc    = "numeric",
           thin   = "numeric",
           nchain = "numeric"))

### Hyper ----------------------------------------------------
##' An S4 class to represent an object storing posterior samples at the
##' participant and hyper level
##'
##' @slot phi_loc posterior samples for the location parameters
##' @slot phi_sca posterior samples for the scale parameters
##' @slot summed_log_prior summed log prior likelihoods for phi.
##' @slot log_likelihoods log likelihoods for phi
##' @slot prior_loc a S4 prior object for the location parameters
##' @slot prior_sca a S4 prior object for the scale parameters
##' @slot start the index of starting sample
##' @slot npar number of parameters
##' @slot pnames parameter names
##' @slot nmc number of Monte Carlo samples
##' @slot thin thinning length
##' @slot nchain number of Markov chains
##' @slot individuals a list storing posterior samples for each individual participant
##' @slot snames names of individual participants
##' @seealso \code{\link{posterior-class}}
##' @export
setClass("hyper",
         slot = c(phi_loc = "array",
                  phi_sca = "array",
                  summed_log_prior = "matrix",
                  log_likelihoods  = "matrix",
                  prior_loc = "prior",
                  prior_sca = "prior",
                  start = "numeric",
                  npar = "numeric",
                  pnames = "character",
                  nmc = "numeric",
                  thin = "numeric",
                  nchain = "numeric",
                  individuals = "list",
                  snames = "character")
)


#### Automatic Sampling  ---------------------------------------
##' Convergence Diagnosis
##'
##' These functions test whether Markov chains are converged .
##'
##' \code{isstuck} tests whether a chain hovers around a region significantly
##' deviates from other its peers.
##'
##' \code{PickStuck} calculate each chain separately for the mean (across
##' MC samples) of posterior log likelihood. If the difference of the means and
##' the median (across chains) of the mean of posterior log likelihood
##' is greater than the value set in \code{cut}, chains are considered stuck.
##' The default value for \code{cut} is 10. The user should consider their
##' situatin to set the cut value.
##'
##' \code{unstick} removes stuck chains from posterior samples (not well tested).
##'
##' \code{ismixed} tests whether the potential scale reduction factor for a
##' model fit is lower than a criterion, defined by \code{cut}.
##'
##' \code{iseffective} testes whether posterior samples are enough adjusted
##' autocorrelation.
##'
##' \code{CheckConverged} is a wrapper function running the four checking
##' functions, \code{isstuck}, \code{isflat}, \code{ismixed} and \code{iseffective}.
##'
##' @param x posterior samples
##' @param hyper whether x are hierarhcial samples
##' @param cut a criterion for deciding whether chains get stuck
##' (\code{isstuck}); whether chains are not flat
##' (using median or IQR \code{isflat}); whether chains are well mixed
##' \code{ismixed}.
##' @param start start to evaluate from which iteration.
##' @param end end at which iteration for evaeuation.
##' @param verbose a boolean switch to print more information
##' @param digits print how many digits. Default is 2
##'
##' @param p1 the range of the head of MCMC chains
##' @param p2 the range of the tail of the MCMC chains
##' @param cut_scale Use IQR to decide whether chains are not flat
##' @param ... other additional arguments
##'
##' @return \code{PickStuck} gives an index vector; \code{unstick} gives a
##' posterios samples.
##'
##' @examples
##' model <- BuildModel(
##' p.map     = list(a = "1", v="1", z="1", d="1", sz="1", sv="1", t0="1",
##'                  st0="1"),
##' match.map = list(M = list(s1 = "r1", s2 = "r2")),
##' factors   = list(S = c("s1", "s2")),
##' responses = c("r1","r2"),
##' constants = c(st0 = 0, d = 0, sv = 0, sz = 0),
##' type      = "rd")
##'
##' npar <- model@npar
##' pop.mean  <- c(a=2,   v=4, z=0.5, t0=0.3)
##' pop.scale <- c(a=0.5, v=.5, z=0.1, t0=0.05)
##' pop.prior <- BuildPrior(
##'   dists = rep("tnorm", npar),
##'   p1    = pop.mean,
##'   p2    = pop.scale,
##'   lower = c(0,-5,  0, 0),
##'   upper = c(5, 7,  1, 1))
##'
##' dat <- simulate(model, nsub = 8, nsim = 30, prior = pop.prior)
##' dmi <- BuildDMI(dat, model)
##' ps <- attr(dat, "parameters")
##'
##' p.prior <- BuildPrior(
##'   dists = rep("tnorm", npar),
##'   p1    = pop.mean,
##'   p2    = pop.scale*5,
##'   lower = c(0,-5, 0, 0),
##'   upper = c(5, 7, 1, 1))
##'
##' mu.prior <- BuildPrior(
##'   dists = rep("tnorm", npar),
##'   p1    = pop.mean,
##'   p2    = pop.scale*5,
##'   lower = c(0,-5,  0, 0),
##'   upper = c(5, 7,  1, 1))
##'
##' sigma.prior <- BuildPrior(
##'   dists = rep("beta", npar),
##'   p1    = c(a=1, v=1, z=1, t0=1),
##'   p2    = rep(1, npar),
##'   upper = rep(1, npar))
##'
##' ## Note the names are important
##' priors <- list(pprior=p.prior, location=mu.prior, scale=sigma.prior)
##'
##' \dontrun{
##' Fit hierarchical model ----##
##' fit0 <- StartNewsamples(dmi, priors)
##' fit  <- run(fit0)
##'
##' PickStuck(fit, hyper=TRUE)
##' PickStuck(fit@individuals[[1]])
##' PickStuck(fit)
##'
##' tmp <- PickStuck(fit, hyper=TRUE, verbose=T)
##' tmp <- PickStuck(fit@individuals[[1]], verbose=T)
##' tmp <- PickStuck(fit, verbose=T)
##' isstuck(fit0@individuals[[1]])
##' isstuck(fit@individuals[[1]])
##' isstuck(fit, hyper = TRUE)
##'
##' tmp <- isflat(fit@individuals[[1]])
##' tmp <- isflat(fit@individuals[[1]], verbose = TRUE)
##'
##' tmp <- isflat(fit@individuals[[1]], cut_scale = .25)
##' tmp <- isflat(fit@individuals[[1]], cut_scale = .25, verbose = TRUE)
##'
##' ## Test unstick
##' fit0 <- StartNewsamples(dmi, priors, nmc=50)
##' fit  <- run(fit0, nmc=200)
##' bad <- PickStuck(fit@individuals[[1]], verbose=T)
##' chain_removed <- unstick_one(fit@individuals[[1]], bad)
## 'plot(fit@individuals[[1]])
##' plot(tmp)
##' }
##' @docType methods
##' @rdname PickStuck-methods
##' @export
setGeneric("PickStuck", function(x, ... ) {
  warning("Class ", class(x), " is not defined for PickStuck")
  return(NULL)
})

### PickStuck  ----------------
##' @importFrom matrixStats rowMeans2
##' @importFrom stats median
##' @rdname PickStuck-methods
setMethod("PickStuck", "posterior", function (x, cut = 10, start = 1, end = NA,
                                              verbose = FALSE, digits = 2)
{
  if (is.na(end)) end <- slot(x, "nmc")
  if (end <= start) stop("End must be greater than start")
  iter   <- start:end
  nchain <- slot(x, "nchain")

  ## nchain x nmc
  lp  <- x@summed_log_prior[,iter] +  x@log_likelihoods[,iter]
  alp <- matrixStats::rowMeans2(lp) ## average across samples
  names(alp) <- 1:nchain
  dev <- -( sort(alp) - median(alp) )
  bad <- as.numeric(names(dev)[dev > cut])

  if (verbose)
  {
    cat("Deviation of mean chain log-likelihood from median of means\n")
    print(round(dev, digits))
    cat("Bad chains (deviance >", cut, "): ")
    if (length(bad) == 0) { cat("None\n") } else { cat(bad, "\n") }
  }
  return(bad)

})

##' @rdname PickStuck-methods
setMethod("PickStuck", "list", function (x, cut = 10, start = 1,
                                         end = NA, verbose = FALSE, digits = 2)
{
  sam1 <- x[[1]]
  if (is.na(end)) end <- slot(sam1, "nmc")
  if (end <= start) stop("end must be greater than start")
  bad <- sapply(x, PickStuck, cut, start, end, verbose, digits)
  return(bad)
})

##' @importFrom matrixStats rowMeans2
##' @importFrom stats median
##' @rdname PickStuck-methods
setMethod("PickStuck", "hyper", function(x, hyper = TRUE, cut = 10, start = 1,
                                         end = NA, verbose = FALSE, digits = 2)
{
  # stucks <- PickStuck(x, cut, start, end, FALSE, digits)

  if (is.na(end)) end <- x@nmc

  if (hyper) {
    if (end <= start) stop("End must be greater than start")
    iter <- start:end
    nchain <- x@nchain

    pll <- x@log_likelihoods[, iter] + x@summed_log_prior[, iter]
    mean.ll <- matrixStats::rowMeans2(pll)
    names(mean.ll) <- 1:nchain
    dev <- -(sort(mean.ll) - median(mean.ll))
    bad <- as.numeric(names(dev)[dev > cut])

    if (verbose) {
      message("Deviation of mean chain log-likelihood of phi from median of means\n")
      print(round(dev,digits))
      cat("Bad chains (deviance >", cut, "): ")
      if (length(bad) == 0) { cat("None\n") } else { cat(bad, "\n") }
    }
    return(bad)


  } else {
    names(x@individuals) <- x@snames
    ## x@individuals is a class list, so the following will call PickStuck.list
    out <- PickStuck(x@individuals, cut, start, end, digits)
  }
  return(out)
})


### isstuck ----------------
##' @docType methods
##' @rdname PickStuck-methods
##' @export
setGeneric("isstuck", function(x, ... ) {
  warning("Class ", class(x), " not defined for isstuck")
  return(NULL)
})

##' @rdname PickStuck-methods
setMethod("isstuck", "posterior", function(x, cut = 10,
                                           start = 1, end = NA, verbose = FALSE)
{
  stucks <- PickStuck(x, cut, start, end, FALSE)
  out <- length(stucks) != 0
  if (verbose)
  {
    message("Found stuck chains: ", appendLF = FALSE)
    if (out) cat (stucks,"\n") else message("none\n")
  }
  return(out)
})

##' @rdname PickStuck-methods
setMethod("isstuck", "list", function(x, cut=10, start=1, end=NA,
                                      verbose=FALSE, digits=2)
{
  if(is.na(end)) end <- x[[1]]@nmc
  snames <- names(x)
  stucks <- PickStuck(x, cut, start, end, digits)
  v <- sapply(stucks, function(xx){ length(xx) == 0 })

  if (verbose) {
    if (all(v)) {
      message("passed")
    } else {
      message("List of participants showing stuck chains: ", appendLF = FALSE)
      cat(snames[!v])
    }
  }
  return(all(v))
})

##' @rdname PickStuck-methods
setMethod("isstuck", "hyper", function(x, hyper=TRUE, cut=10, start=1, end=NA,
                                       verbose=FALSE, digits=2)
{
  if(is.na(end)) end <- x@nmc
  if (hyper) {

    stucks <- PickStuck(x, hyper, cut, start, end, FALSE, digits)
    out <- length(stucks) != 0
    if (verbose)
    {
      message("Found stuck chains: ", appendLF = FALSE)
      if (out) cat (stucks,"\n") else message("none")
    }
  } else {
    snames <- x@snames

    stucks <- PickStuck(x@individuals, cut, start, end, digits)
    v <- sapply(stucks, function(xx){ length(xx) == 0 })
    if (verbose) {
      if (all(v)) {
        message("passed")
      } else {
        message("List of participants showing stuck chains: ", appendLF = FALSE)
        cat(snames[!v])
      }
    }
    out <- all(v)
  }

  return(out)
})


##' @docType methods
##' @rdname PickStuck-methods
##' @export
setGeneric("isflat", function(x, ... ) {
  warning("Class ", class(x), " not defined for isflat")
  return(NULL)
}
)

### isflat  ------------------------
##' @rdname PickStuck-methods
setMethod("isflat", "posterior", function(x, p1 = 1/3, p2 = 1/3, cut = 0.25,
                                          cut_scale = Inf, verbose = FALSE,
                                          digits = 2)
{
  # x <- fit0@individuals[[1]]
  # p1 <- 1/3
  # p2 <- 1/3
  # cut <- .25
  # verbose <- T
  # digits <- 2
  # stage1
  # cut_scale <- .25
    if (is.finite(cut_scale)) {
        stage1 <- isflat_location(x, p1, p2, cut, verbose, digits)
        tmp <- isflat_scale(x, p1, p2, cut_scale, verbose, digits, stage1)
        out <- !tmp
    } else {
        stage1 <- isflat_location(x, p1, p2, cut, verbose, digits)
        out <- !stage1[[1]]
    }
    return(out)
})

##' @docType methods
##' @rdname PickStuck-methods
setMethod("isflat", "list", function(x, p1 = 1/3, p2 = 1/3, cut = 0.25,
                                     cut_scale = Inf, verbose = FALSE,
                                     digits = 2)
{
  # x <- fit0
  # p1 <- 1/3
  # p2 <- 1/3
  # cut <- .25
  # cut_scale <- .25
  # verbose <- T
  # digits <- 2
  nsub <- length(x)

    if (is.finite(cut_scale)) {
        stage1 <- lapply(x, isflat_location, p1, p2, cut, FALSE, digits)
        out <- lapply(1:nsub, function(k){
          if(verbose) cat("Participant", k, "\n")
          isflat_scale(x[[k]], p1, p2, cut_scale, verbose, digits, stage1[[k]])
        })

    } else {
      out <- lapply(1:nsub, function(k){
        if(verbose) cat("Participant", k, "\n")
        isflat_location(x[[k]], p1, p2, cut, verbose, digits)
      })
    }


    return(sapply(out, function(xx) { !xx[[1]] }))
})


##' @importFrom stats IQR
##' @importFrom stats median
isflat_location <- function(x, p1, p2, cut, verbose, digits)
{
  # cut <- .25
  # cut_scale <- Inf
  theta  <- slot(x, "theta")
  npar   <- slot(x, "npar")
  pnames <- slot(x, "pnames")
  nchain <- slot(x, "nchain")
  nmc    <- slot(x, "nmc")
  nsam   <- nchain * nmc
  mat    <- matrix(aperm(theta, c(2, 3, 1)), ncol = npar) ## nsam x npar

  head_len <- round(nsam * p1)
  tail_len <- round(nsam * p2)

  ## Change in median relative to IQR (ie robust SD)
  zs <- apply (mat, 2, function(xx)
  {
    head_sam <- median(xx[1:head_len])
    tail_sam <- median(xx[(nsam - tail_len):nsam])
    abs(head_sam - tail_sam) / IQR(xx)
  })

  names(zs) <- paste0("m_", pnames)
  fail <- any(zs > cut)

  if (fail) {
    znames <- names(zs)
    parameter_value <- round(max(zs), digits)
    which_parameter <- znames[which.max(zs)]
    info <- paste(which_parameter, "=", parameter_value)
  } else {
    info <- ""
  }

  if (verbose)
  {
    if (!fail) {
      message("location parameters passed")
    } else {
      print(round(zs, digits))
      cat(paste("Parameter(s) with maximum location zs:", info, "\n\n"))
    }
  }

  return(list(fail, info))
}

##' @importFrom stats IQR
isflat_scale <- function(x, p1 = 1/3, p2 = 1/3, cut = 0.25,
                         verbose = FALSE, digits=2, location_info)
{
  # location_info <- stage1
  theta  <- slot(x, "theta")
  npar   <- slot(x, "npar")
  pnames <- slot(x, "pnames")
  nchain <- slot(x, "nchain")
  nmc    <- slot(x, "nmc")
  nsam   <- nchain * nmc
  mat    <- matrix(aperm(theta, c(2, 3, 1)), ncol = npar) ## nsam x npar
  head_len <- round(nsam * p1)
  tail_len <- round(nsam * p2)

  ## Change to IQR reltiave to overall IQR
  zs <- apply(mat, 2, function(xx)
  {
    head_sam <- IQR(xx[1:head_len])
    tail_sam <- IQR(xx[(nsam - tail_len):nsam])
    abs(head_sam - tail_sam)/IQR(xx)
  })

  names(zs) <- paste0("s_", pnames)

  if ( any(zs > cut) ) {
    znames <- names(zs)
    parameter_value <- round(max(zs), digits)
    which_parameter <- znames[which.max(zs)]
    scale_info <- paste(which_parameter, "=", parameter_value)

    fail <- location_info[[1]] | any(zs > cut)

  } else {
    fail <- location_info[[1]]
  }

  if (verbose)
  {
    if (fail) {
      print(round(zs, digits))
      cat(paste("Parameter(s) with maximum scale zs:", scale_info, "\n\n"))
    } else {
      message("scale parameters passed")
    }
  }

  return(fail)
}

### Brook-Gelman PSRF ------------------------
##' Potential scale reduction factor
##'
##' \code{gelman} function calls the function, \code{gelman.diag} in the
##' \pkg{coda} package to calculates PSRF.
##'
##' @param x posterior samples
##' @param hyper a Boolean switch, indicating posterior samples are from
##' hierarchical modeling
##' @param start start iteration
##' @param end end iteration
##' @param conf confident inteval
##' @param multivariate multivariate Boolean switch
##' @param subchain whether only calculate a subset of chains
##' @param verbose print more information
##' @param digits print out how many digits
##' @param ... other additional arguments
##' @export
##' @examples
##' \dontrun{
##' rhat1 <- gelman(hsam);
##' rhat2 <- gelman(hsam, end = 51);
##' rhat3 <- gelman(hsam, conf = .90);
##' rhat7 <- gelman(hsam, subchain = TRUE);
##' rhat8 <- gelman(hsam, subchain = 1:4);
##' rhat9 <- gelman(hsam, subchain = 5:7, digits = 1, verbose = TRUE);
##' }
##' @docType methods
##' @rdname gelman-methods
##' @export
setGeneric("gelman", function(x, ... ) {
  warning("Class ", class(x), " not defined for gelman")
  return(NULL)
})

##' @rdname gelman-methods
setMethod("gelman", "posterior", function (x, start = 1,
                                           end = NA, conf = 0.95,
                                           multivariate = TRUE, subchain = NA,
                                           digits = 2, verbose = FALSE)
{
  # x <- fit[[1]]
  # start <- 1
  # end <- NA
  # subchain <- T
  # nsubchain <- 3
  # multivariate <- T
  # conf <- .95
  # digits <- 2
  # subchain <- 1:3
  # subchain <- NA

  nchain <- slot(x, "nchain")
  if (is.na(end)) end <- slot(x, "nmc")
  if (x@nchain < 2) stop("Minimum two chains needed")

  if (!anyNA(subchain)) x <- sample_subchains(x, subchain)

  psrf <- gelman_diag(x, iter = start:end, conf, multivariate)

  if (verbose) {
    cat("Multivariate psrf: \n")
    print(round(psrf$psrf, digits))
  }
  return(psrf)
})

##' @rdname gelman-methods
setMethod("gelman", "list", function (x, start = 1,
                                      end = NA, conf = 0.95,
                                      multivariate = TRUE, subchain = NA,
                                      digits = 2, verbose = FALSE)
{

  out <- lapply(x, gelman, start, end, conf, multivariate, subchain,
                digits, FALSE) ## turn off verbose in gelman_one

  if (verbose && multivariate) {
    message("Diagnosing convergence of multi-participant separately")
    tmp   <- sapply(out, function(x){x$mpsrf})
    mpsrf <- c(mean(tmp), sort(tmp))
    names(mpsrf) <- c("Mean", names(sort(tmp)))
    print( round(mpsrf, digits) )
  }
  return(out)
})

##' @rdname gelman-methods
setMethod("gelman", "hyper", function(x, hyper = TRUE, start = 1, end = NA,
                                      conf = 0.95, multivariate = TRUE,
                                      subchain = NA, digits = 2,
                                      verbose = FALSE) {

  snames <- x@snames
  nsub   <- length(snames)

  if (hyper) {

    if (is.na(end)) end <- x@nmc
    psrf_hyper <- gelman_hyper(x, start, end, conf, multivariate, subchain)

    psrf <- lapply(x@individuals, gelman, start, end, conf,
                   multivariate, subchain, digits, FALSE)

    tmp   <- sapply(psrf, function(x){x$mpsrf})
    mpsrf <- c(psrf_hyper$mpsrf, mean(tmp), sort(tmp))

    if (verbose && multivariate) {
      message("Convergence checks for hyper and participant-level parameters")
      names(mpsrf) <- c("Hyper", "Mean", x@snames[order(tmp)] )
      print( round(mpsrf, digits) )
    }

    out <- vector(mode="list", length=1+nsub)
    for(i in 1:nsub) out[[i]] <- psrf[[i]]

    out[[1+nsub]] <- psrf_hyper
    names(out) <- c(snames, "Hyper")

  } else {
    out <- lapply(x@individuals, gelman, start, end, conf, multivariate,
                  subchain, digits, FALSE)

    if (verbose && multivariate) {
      message("Diagnosing convergence of multi-participant separately")
      tmp   <- sapply(out, function(x){x$mpsrf})
      mpsrf <- c(mean(tmp), sort(tmp))
      names(mpsrf) <- c( "Mean", x@snames[order(tmp)] )
      print( round(mpsrf, digits) )
    }
  }

  return(out)

})

##' @importFrom stats qf
gelman_hyper <- function(x, start, end, conf, multivariate, subchain) {

  phi <- extractPhi(x, start, end, subchain)
  npar   <- attr(phi, "npar")
  nchain <- attr(phi, "nchain")
  iter   <- attr(phi, "iter")
  niter  <- length(iter)
  pnames <- attr(phi, "pnames")

  tmp1 <- sapply(phi, var)
  tmp2 <- matrixStats::rowMeans2(tmp1)
  W  <- matrix(tmp2, nrow = npar)
  S2 <- array(tmp1, dim = c(npar, npar, nchain))

  vbar  <- lapply(1:nchain, function(k){ matrixStats::colMeans2(phi[[k]]) })

  xbar  <- matrix( unlist(vbar), nrow=npar)
  B     <- niter * var(t(xbar))
  mpsrf <- get_mpsrf(npar, niter, W, B, multivariate)

  w <- diag(W)
  b <- diag(B)
  s2 <- matrix(apply(S2, 3, diag), nrow = npar, ncol = nchain)
  muhat <- apply(xbar, 1, mean)
  var.w <- apply(s2, 1, var)/nchain
  var.b <- (2 * b^2)/(nchain - 1)
  cov.wb <- (niter/nchain) * diag(var(t(s2), t(xbar^2)) - 2 *
                                    muhat * var(t(s2), t(xbar)))

  ## Part 2
  V <- (niter - 1) * w/niter + (1 + 1/nchain) * b/niter
  var.V <- ((niter - 1)^2 * var.w + (1 + 1/nchain)^2 * var.b +
              2 * (niter - 1) * (1 + 1/nchain) * cov.wb)/niter^2
  df.V <- (2 * V^2)/var.V
  df.adj <- (df.V + 3)/(df.V + 1)
  B.df <- nchain - 1
  W.df <- (2 * w^2)/var.w
  R2.fixed <- (niter - 1)/niter
  R2.random <- (1 + 1/nchain) * (1/niter) * (b/w)
  R2.estimate <- R2.fixed + R2.random
  R2.upper <- R2.fixed + stats::qf( .5*(1 + conf), B.df, W.df ) * R2.random
  psrf <- cbind(sqrt(df.adj * R2.estimate), sqrt(df.adj * R2.upper))

  dimnames(psrf) <- list(pnames, c("Point est.", "Upper C.I."))

  ## Part 3
  return(list(psrf = psrf, mpsrf = mpsrf))
}

##' @importFrom stats qf
gelman_diag <- function(x, iter, conf, multivariate) {

  # iter <- start:end
  niter  <- length(iter)
  npar   <- slot(x, "npar")
  pnames <- slot(x, "pnames")
  nchain <- slot(x, "nchain")

  ## Part 1
  WS2 <- getW(x)
  W   <- WS2[[1]]
  S2  <- WS2[[2]]
  xbar <- get_xbar(x)
  B    <- niter * var(t(xbar))
  mpsrf <- get_mpsrf(npar, niter, W, B, multivariate)

  w <- diag(W)
  b <- diag(B)
  s2 <- matrix(apply(S2, 3, diag), nrow = npar, ncol = nchain)
  muhat <- apply(xbar, 1, mean)
  var.w <- apply(s2, 1, var)/nchain
  var.b <- (2 * b^2)/(nchain - 1)
  cov.wb <- (niter/nchain) * diag(var(t(s2), t(xbar^2)) - 2 *
                                    muhat * var(t(s2), t(xbar)))

  ## Part 2
  V <- (niter - 1) * w/niter + (1 + 1/nchain) * b/niter
  var.V <- ((niter - 1)^2 * var.w + (1 + 1/nchain)^2 * var.b +
              2 * (niter - 1) * (1 + 1/nchain) * cov.wb)/niter^2
  df.V <- (2 * V^2)/var.V
  df.adj <- (df.V + 3)/(df.V + 1)
  B.df <- nchain - 1
  W.df <- (2 * w^2)/var.w
  R2.fixed <- (niter - 1)/niter
  R2.random <- (1 + 1/nchain) * (1/niter) * (b/w)
  R2.estimate <- R2.fixed + R2.random
  R2.upper <- R2.fixed + stats::qf( .5*(1 + conf), B.df, W.df ) * R2.random
  psrf <- cbind(sqrt(df.adj * R2.estimate), sqrt(df.adj * R2.upper))


  ## Part 3
  dimnames(psrf) <- list(pnames, c("Point est.", "Upper C.I."))
  return(list(psrf = psrf, mpsrf = mpsrf))
}

sample_subchains <- function(x, subchain) {
  if ( any(x@nchain < subchain) ) stop("Chains not in the range.")

  chains <- subchain
  cat("Calculate chains:", chains, "\n")

  slot(x, "nchain") <- length(subchain)
  theta <- slot(x, "theta")
  lp    <- slot(x, "summed_log_prior")
  ll    <- slot(x, "log_likelihoods")
  slot(x, "theta") <- theta[, chains,]
  slot(x, "summed_log_prior") <- lp[chains,]
  slot(x, "log_likelihoods") <- ll[chains,]
  return(x)
}


### ismixed ------------------
##' @docType methods
##' @rdname PickStuck-methods
##' @export
setGeneric("ismixed", function(x, ...){
  warning("Class ", class(x), " not defined for ismixed")
  return(NULL)
})

##' @rdname PickStuck-methods
setMethod("ismixed", "posterior",  function(x, cut = 1.10, verbose = FALSE)
{

  pnames <- slot(x, "pnames")
  gnames <- c("mpsrf", pnames)

  tmp <- gelman(x)
  gds <- c(tmp$mpsrf, tmp$psrf[,1])
  names(gds) <- gnames
  out <- max(gds) > cut

  if (verbose) {

    if (!out) { cat("passed\n") } else {
      cat("Mixing check:\n")
      print(round(gds, 2))
    }
  }
  return(!out)
})

### Effective Size  -----------------------------
##' Effective Sample Size
##'
##' Posterior sample size adjusted for autocorrelation. The function is based
##' on the effectiveSize function in \code{coda} package.
##'
##' \code{hyper} argument does not work for list class (i.e., posterior
##' samples from a fixed-effect model fit).
##'
##' @param x posterior samples
##' @param hyper a Boolean switch to calculate phi
##' @param start start from iteration
##' @param end end at which iteraton
##' @param subchain calculate a subset of chains. This must be an integer vector
##' @param digits printing how many digits
##' @param verbose printing more information
##' @param ... other additional arguments
##' @references
##' \enumerate{
##' Plummer, M. Best, N., Cowles, K., Vines, K., Sarkar, D., Bates, D., Almond, R., & Magnusson, A. (2019). R package 'coda' \url{https://cran.r-project.org/web/packages/coda/}
##' }
##' @examples
##' #################################40
##' ## effectiveSize example
##' #################################40
##' \dontrun{
##' cat("Class:", class(fit), "\n")
##' es1 <- effectiveSize(fit, hyper=TRUE, verbose=FALSE)
##' es1 <- effectiveSize(fit, hyper=TRUE, verbose=TRUE)
##'
##' es1 <- effectiveSize(fit, hyper=TRUE, verbose=FALSE, subchain=7:9)
##' es1 <- effectiveSize(fit, hyper=TRUE, verbose=TRUE, subchain=7:9)
##'
##' es1 <- effectiveSize(fit, hyper=FALSE, verbose=FALSE)
##' es1 <- effectiveSize(fit, hyper=FALSE, verbose=TRUE)
##'
##' es1 <- effectiveSize(fit, hyper=FALSE, verbose=FALSE, subchain=4:6)
##' es1 <- effectiveSize(fit, hyper=FALSE, verbose=TRUE, subchain=4:6)
##'
##' cat("Starting a new fixed-effect model fit: \n")
##' fit0 <- StartNewsamples(dmi, p.prior, ncore=4)
##' fit  <- run(fit0, ncore=4)
##'
##' cat("Class:", class(fit), "\n")
##' es1 <- effectiveSize(fit, verbose=FALSE)
##' es1 <- effectiveSize(fit, verbose=TRUE)
##' es1 <- effectiveSize(fit, verbose=FALSE, subchain=4:6)
##' es1 <- effectiveSize(fit, verbose=TRUE, subchain=4:6)
##'
##' }
##'
##' @docType methods
##' @rdname effectiveSize-methods
##' @export
setGeneric("effectiveSize", function(x, ...)
{
  warning("Class ", class(x), " not defined for effectiveSize")
  return(NULL)
})

##' @importFrom matrixStats colVars
##' @importFrom matrixStats colSums2
##' @rdname effectiveSize-methods
setMethod("effectiveSize", "hyper",
          function(x, hyper = TRUE, start=1, end=NA, subchain=NA, digits = 2,
                   verbose = FALSE) {

  if (hyper) {
    if (is.na(end)) end <- x@nmc
    phi <- extractPhi(x, start, end, subchain)
    npar   <- attr(phi, "npar")
    nchain <- attr(phi, "nchain")
    iter   <- attr(phi, "iter")
    niter  <- length(iter)
    pnames <- attr(phi, "pnames")

    v <- lapply(1:nchain, function(k){
      spec <- spectrum0_ar( phi[[k]][iter,] )$spec     ## nmc x npar
      ans  <- niter * matrixStats::colVars( phi[[k]][iter,] ) / spec
      ans
    })

    out <- matrixStats::colSums2( do.call("rbind", v) )
    names(out) <- pnames
    if (verbose) print(round(out, digits))
  } else {
    out <- effectiveSize(x@individuals, start, end, subchain, digits, verbose)
  }
  return(out)
})


##' @importFrom matrixStats rowMeans2
##' @importFrom matrixStats rowSds
##' @importFrom matrixStats rowMaxs
##' @importFrom matrixStats rowMins
##' @rdname effectiveSize-methods
setMethod("effectiveSize", "list",
          function(x, start=1, end=NA, subchain=NA, digits=2, verbose=FALSE)
{
  # x <- fit
  # end <- NA
  # start <- 1
  # digits <- 2
  if (is.na(end)) end <- slot(x[[1]], "nmc")
  pnames <- slot(x[[1]], "pnames")

  out <- sapply(x, function(xx) {
    effectiveSize(xx, start, end, subchain, digits, FALSE)
  })

  if (verbose) {
    p1 <- matrixStats::rowMeans2(out)
    p2 <- matrixStats::rowSds(out)
    p3 <- matrixStats::rowMaxs(out)
    p4 <- matrixStats::rowMins(out)
    printthis <- rbind(p1, p2, p3, p4)
    rownames(printthis) <- c("MEAN", "SD", "MAX", "MIN")
    colnames(printthis) <- pnames
    print(round(printthis, digits))
  }
  return(out)
})

##' @importFrom matrixStats colVars
##' @importFrom matrixStats colSums2
##' @rdname effectiveSize-methods
setMethod("effectiveSize", "posterior",
          function(x, start=1, end=NA, subchain=NA, digits=2, verbose=FALSE)
{
  if (is.na(end)) end <- slot(x, "nmc")
  iter   <- start:end
  niter  <- length(iter)

  if (!anyNA(subchain)) x <- sample_subchains(x, subchain)

  theta  <- slot(x, "theta") ## npar x nchain x nmc
  nchain <- slot(x, "nchain")
  pnames <- slot(x, "pnames")

  v <- lapply(1:nchain, function(k){
    spec <- spectrum0_ar( t(theta[,k,iter]) )$spec
    ans  <- niter * matrixStats::colVars( t(theta[,k,iter]) ) / spec
    ans ## This mini function returns this
  })

  out <- matrixStats::colSums2( do.call("rbind", v) )
  names(out) <- pnames
  if (verbose) print(round(out, digits))
  return(out)
})

### DIC ---------------
##' Deviance Information Criteria
##'
##' Calculate DIC and BPIC.
##'
##' This function implements three different definitions of the "effective
##' number of parameters of the model". First is from Spiegelhalter et al
##' (2002, p. 587), "... that pD can be considered as a 'mean deviance minus
##' the deviance of the means'". Second is from Gelman et al (2014, p. 173,
##' equation 7.10), and third subtracts the minimal value of the deviance from
##' the mean of the deviance.
##'
##' @param object posterior samples from one participant
##' @param start start from which iteration.
##' @param end end at which iteration. For example, set
##' \code{start = 101} and \code{end = 1000}, instructs the function to
##' calculate from 101 to 1000 iteration.
##' @param BPIC a Boolean switch to calculate BPIC, instead of DIC
##' @param ... other plotting arguments passing through dot dot dot.
##' @examples
##' ## Calculate DIC from data of one participant
##' \dontrun{
##' model <- BuildModel(
##'   p.map     = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1",
##'                   st0 = "1"),
##'   match.map = list(M = list(s1 = 1, s2 = 2)),
##'   factors   = list(S = c("s1", "s2")),
##'   constants = c(st0 = 0, sd_v = 1),
##'   responses = c("r1", "r2"),
##'   type      = "norm")
##'
##' p.vector <- c(A = .75, B = 1.25, t0 = .15, mean_v.true = 2.5,
##'               mean_v.false = 1.5)
##' ntrial <- 50
##' dat <- simulate(model, nsim = ntrial, ps = p.vector)
##' dmi <- BuildDMI(dat, model)
##'
##' p.prior <- BuildPrior(
##'   dists = c("tnorm", "tnorm", "beta", "tnorm", "tnorm"),
##'   p1    = c(A = 1, B = 1, t0 = 1, mean_v.true = 1, mean_v.false = 1),
##'   p2    = c(1,  1,  1, 1, 1),
##'   lower = c(rep(0, 3),  rep(NA, 2)),
##'   upper = c(rep(NA, 2), 1, rep(NA, 2)))
##'
##' ## Sampling
##' fit0 <- StartNewsamples(dmi, p.prior)
##' fit  <- run(fit0, thin = 8)
##'
##' DIC(fit)
##' DIC(fit)
##' DIC(fit, start=100, end=200)
##' DIC(fit, BPIC=TRUE)
##' DIC(fit, BPIC=TRUE, start=201, end=400)
##' }
##'
##' ## Calculate DICs from data of 8 participant
##' \dontrun{
##' model <- BuildModel(
##' p.map     = list(a = "1", v = "F", z = "1", d = "1", sz = "1", sv = "1",
##'                  t0 = "1", st0 = "1"),
##' match.map = list(M = list(s1 = "r1", s2 = "r2")),
##' factors   = list(S = c("s1", "s2"), F = c("f1", "f2")),
##' constants = c(st0 = 0, d = 0),
##' responses = c("r1", "r2"),
##' type      = "rd")
##' npar <- length(Get_pnames(model))
##'
##' ## Population distribution
##' pop.mean  <- c(a=2,   v.f1=4,  v.f2=3,  z=0.5, sz=0.3, sv=1,  t0=0.3)
##' pop.scale <- c(a=0.5, v.f1=.5, v.f2=.5, z=0.1, sz=0.1, sv=.3, t0=0.05)
##' pop.prior <- BuildPrior(
##'   dists = rep("tnorm", npar),
##'   p1    = pop.mean,
##'   p2    = pop.scale,
##'   lower = c(0,-5, -5, 0, 0, 0, 0),
##'   upper = c(5, 7,  7, 1, 2, 1, 1))
##'
##' ## Simulate some data
##' dat <- simulate(model, nsub = 8, nsim = 10, prior = pop.prior)
##' dmi <- BuildDMI(dat, model)
##' ps <- attr(dat, "parameters")
##'
##' p.prior <- BuildPrior(
##'   dists = rep("tnorm", npar),
##'   p1    = pop.mean,
##'   p2    = pop.scale*5,
##'   lower = c(0,-5, -5, 0, 0, 0, 0),
##'   upper = c(5, 7,  7, 1, 2, 2, 1))
##'
##' ## Sampling
##' fit0 <- StartNewsamples(dmi, p.prior, ncore=1)
##' fit  <- run(fit0, ncore=4)  ## No printing when running in RStudio
##'
##' ## Calculate DIC for participant 1
##' DIC(fit[[1]])
##'
##' ## Calculate all participants
##' res <- DIC(fit)
##'
##' ## BPIC
##' res <- DIC(fit, BPIC = TRUE)
##' }
##'
##' @docType methods
##' @rdname DIC-methods
##' @export
setGeneric("DIC", function(object, ... ) {
  warning("Class ", class(object), " not defined for DIC")
  return(NULL)
})

##' @rdname DIC-methods
setMethod("DIC", "posterior", function(object, start = 1, end = NA,
                                       BPIC = FALSE)
{
  if (is.na(end)) end <- object@nmc
  if ( end <= start ) { stop("End must be greater than start") }

  ds <- deviance_model(object, start, end)
  ## Three different definitions of "the effective number of parameters of the
  ## model" (1) Spiegelhalter et al (2002, p. 587); (2) ???; (3) Gelman et al
  ## (2014, p. 173, equation 7.10)
  pds <- list(Pmean = ds$meanD-ds$Dmean, ## (1)
              Pmin  = ds$meanD-ds$minD,
              Pvar  = 2*ds$varD)         ## (3) Never been used?

  if (ds$minD < ds$Dmean) { pd <- pds$Pmin }
  else                    { pd <- pds$Pmean }

  ## Both return the same number
  ## pd + ds$meanD    ## Spiegelhalter et al's method
  ## ds$Dmean + 2*pd  ## Gelman et al (2014, p. 173, the equation after 7.10)
  ## ds$meanD + 2*pd  ## Ando (2007)
  if (BPIC) { out <- ds$meanD+2*pd }
  else      { out <- ds$meanD+pd }

  return(out)
})

##' @rdname DIC-methods
setMethod("DIC", "list", function(object, start = 1, end = NA, BPIC = FALSE)
{
  cat("in list class")
  out <- sapply(object, DIC, start, end, BPIC)
  if (BPIC) { message("Summed BPIC: ", appendLF = FALSE) }
  else      { message("Summed DIC: ", appendLF = FALSE) }
  cat(sum(out), "\n")
  out
})

##' @rdname DIC-methods
setMethod("DIC", "hyper", function(object, start = 1, end = NA, BPIC = FALSE)
{
  names(object@individuals) <- object@snames
  ## This will go to DIC posterior
  out <- sapply(object@individuals, DIC, start, end, BPIC)
  if (BPIC) { message("Summed BPIC: ", appendLF = FALSE) }
  else      { message("Summed DIC: ", appendLF = FALSE) }
  cat(sum(out), "\n")
  out
})



### loglik ----------
##' Extract Posterior Log-Likelihood
##'
##' This function is to extract posterior log-likelihood in a "model" object.
##'
##' @param object posterior samples
##' @param hyper whether to summarise hyper parameters
##' @param start start from which iteration.
##' @param end end at which iteration. For example, set
##' \code{start = 101} and \code{end = 1000}, instructs the function to
##' calculate from 101 to 1000 iteration.
##' @param ... other arguments passing through dot dot dot.
##' @docType methods
##' @rdname logLik-methods
##' @export
setGeneric("logLik", function(object, ... ) {
  warning("Class ", class(object), " not defined for logLik")
  return(NULL)
})

##' @rdname logLik-methods
setMethod("logLik", "posterior", function(object, start = 1, end = NA) {

  if (is.na(end)) end <- object@nmc
  if ( end <= start ) { stop("End must be greater than start") }

  tmp <- dim(object@summed_log_prior)
  if (tmp[2] < tmp[1]) stop("Did you use DMC or earlier version? Remember to transpose the matrices.")

  lp <- object@summed_log_prior[,start:end]
  ll <- object@log_likelihoods[,start:end]
  return(lp + ll)
})

##' @rdname logLik-methods
setMethod("logLik", "list", function(object, start = 1, end = NA)
{
  if (is.na(end)) end <- object@nmc
  if ( end <= start ) stop("End must be greater than start")

  out <- lapply(object, function(x)
  { x@summed_log_prior[, start:end] + x@log_likelihoods[, start:end] })
  return(out)
})

##' @rdname logLik-methods
setMethod("logLik", "hyper", function(object, start = 1, end = NA)
{
  hyper <- attr(object, "hyper")
  if (is.null(hyper)) stop("Samples are not from a hierarhcial model fit")
  if (is.na(end)) end <- object$nmc
  if (end <= start) stop("End must be greater than start")

  lp <- hyper$h_summed_log_prior[,start:end]
  ll <- hyper$h_log_likelihoods[,start:end]
  out <-  lp + ll
  return(out)
})

### Print ------------------------------------------------------------
##' ggdmc Printing Methods
##'
##' The function is an extension of the print function in \code{base} pacakge.
##' It prints a model object set up by \code{BuildModel} and a prior object
##' set up by \code{BuildPrior}.
##'
##' The print method for a prior object merely rearranges a prior object
##' as a data frame for the inspection convenience.
##'
##' @param x a model object.
##' @param ps a parameter vector
##' @param ... Additional argument passing via dot dot dot.
##'
##' @return The original model object, a list of parameter matrices or a prior
##' matrix
##'
##' @export
##' @docType methods
##' @rdname print-methods
##' @examples
##' model <- BuildModel(
##'           p.map     = list(A = "1", B = "1", t0 = "1", mean_v = "M",
##'                       sd_v = "1", st0 = "1"),
##'           match.map = list(M = list(s1 = 1, s2 = 2)),
##'           factors   = list(S = c("s1", "s2")),
##'           constants = c(st0 = 0, sd_v = 1),
##'           responses = c("r1", "r2"),
##'           type      = "norm")
##'
##' p.vector <- c(A = .75, B = 1.25, t0 = .15, mean_v.true = 2.5,
##'               mean_v.false = 1.5)
##'
##' print(model)
##' print(model, ps=p.vector)
##'
##' dat <- simulate(model, nsim = 10, ps = p.vector);
##' dmi <- BuildDMI(dat, model)
##' p.prior <- BuildPrior(
##'   dists = c("tnorm", "tnorm", "beta", "tnorm", "tnorm"),
##'   p1    = c(A = 1, B = 1, t0 = 1, mean_v.true = 1, mean_v.false = 1),
##'   p2    = c(1,  1,  1, 1, 1),
##'   lower = c(rep(0, 3),  rep(NA, 2)),
##'   upper = c(rep(NA, 2), 1, rep(NA, 2)))
##'
##' print(p.prior)
##'
##' ## A different example printing a prior object
##' pop.mean  <- c(a=1,  v.f1=1,  v.f2=.2, z=.5, sz=.3,  sv.f1=.25, sv.f2=.23,
##'                t0=.3)
##' pop.scale <- c(a=.2, v.f1=.2, v.f2=.2, z=.1, sz=.05, sv.f1=.05, sv.f2=.05,
##'                t0=.05)
##'
##' p.prior <- BuildPrior(
##'   dists = rep("tnorm", 8),
##'   p1    = pop.mean,
##'   p2    = pop.scale,
##'   lower = c(0, -5, -5, 0, 0, 0, 0, 0),
##'   upper = c(2,  5,  5, 1, 2, 2, 1, 1))
##'
##' print(p.prior)
setGeneric("print", function(x, ...) { standardGeneric("print") })


##' @rdname print-methods
setMethod("print", "model", function (x, ps = NULL, ...) {

  if (!is.array( x@model )) stop("model is not an array. Is input a posterior object?")

  if (is.null(ps)) {

    nr <- length(x@responses)

    for (i in 1:nr) {
      dim3 <- x@dimnames[[3]]
      cat(dim3[i], "\n")
      print(x@model[,, i])
    }

    message("Attributes: ")
    print(names(attributes(x)))
    return(invisible(x))

  } else {
    dim1 <- x@dimnames[[1]]

    out <- lapply(dim1, function(xx) {
      print(xx)
      print(TableParameters(ps, xx, x, TRUE))
    })
    return(invisible(out))
  }

} )

##' @rdname print-methods
setMethod("print", "prior", function (x, ...) {

  obj  <- x@priors
  npar <- x@npar

  ncol <- 7
  bucket  <- matrix(numeric(npar*ncol), npar);

  for(i in 1:npar) {
    add1 <- attr(obj[[i]], "dist");
    add1 <- ifelse(add1 == 1, "tnorm",
    ifelse(add1 == 2, "beta",
    ifelse(add1 == 3, "gamma",
    ifelse(add1 == 4, "lnorm",
    ifelse(add1 == 5, "unif",
    ifelse(add1 == 6, "constant",
    ifelse(add1 == 7, "tnorm2",
    ifelse(add1 == 8, "cauchy", NA))))))))


    add2    <- attr(obj[[i]], "untrans");
    if (is.na(add1)) {
      tmp <- unlist(obj[[i]])
      rowObj  <- c(NA, NA, NA, NA, tmp[5], NA, add2)
    } else if (add1 == "constant") {
      tmp <- unlist(obj[[i]])
      rowObj  <- c(c(tmp[1], tmp[2], NA, NA, tmp[5]), add1, add2);
    } else if (add1 == "gamma_l") {
      tmp <- unlist(obj[[i]])
      rowObj <- c(tmp[1:5], "gamma_l", add2)
    } else {
      rowObj  <- c(unlist(obj[[i]]), add1, add2);
    }
    bucket[i,] <- rowObj
  }


  out <- data.frame(bucket)
  names(out) <- c("p1", "p2", "lower", "upper", "lg", "dist", "untrans")
  rownames(out) <- slot(x, "pnames")
  return(out)

})


### Summary ------------------------------------------------------------
##' ggdmc Summary Methods
##'
##' Summarise posterior samples. Note when recovery = TRUE, the prob vector
##' will be fixed at the default values.
##'
##' @param object an object storing posterior samples.
##' @param hyper a Boolean switch to plot hyper parameters
##' @param start start from which iteration.
##' @param end end at which iteration. For example, set
##' \code{start = 101} and \code{end = 1000}, instructs the function to
##' calculate from 101st to 1000th iteration.
##' @param prob a numeric vector, indicating the quantiles to calculate
##' @param recovery a Boolean switch indicating if samples are from a recovery
##' study.
##' @param ps true parameter values.  This is only for recovery studies
##' @param type calculate type 1 = location or type 2 = scale hyper parameters
##' @param verbose print more information
##' @param digits printing digits
##' @param ... Additional argument passing via dot dot dot.
##'
##' @return NULL
##'
##' @export
##' @docType methods
##' @rdname summary-methods
##' @examples
##' \dontrun{
##' model <- BuildModel(
##'      p.map    = list(a = "1", v = "F", z = "1", d = "1", sz = "1", sv = "1",
##'                 t0 = "1", st0 = "1"),
##'     match.map = list(M = list(s1 = "r1", s2 = "r2")),
##'     factors   = list(S = c("s1", "s2"), F = c("f1", "f2")),
##'     constants = c(st0 = 0, d = 0),
##'     responses = c("r1", "r2"),
##'     type      = "rd")
##' npar <- model@npar
##'
##' ## Population distribution
##' pop.mean  <- c(a=2,   v.f1=4,  v.f2=3,  z=0.5, sz=0.3, sv=1,  t0=0.3)
##' pop.scale <- c(a=0.5, v.f1=.5, v.f2=.5, z=0.1, sz=0.1, sv=.3, t0=0.05)
##' pop.prior <- BuildPrior(
##'   dists = rep("tnorm", npar),
##'   p1    = pop.mean,
##'   p2    = pop.scale,
##'   lower = c(0,-5, -5, 0, 0, 0, 0),
##'   upper = c(5, 7,  7, 1, 2, 1, 1))
##'
##' ## Simulate some data
##' dat <- simulate(model, nsub = 30, nsim = 30, prior = pop.prior)
##' dmi <- BuildDMI(dat, model)
##' ps <- attr(dat, "parameters")
##'
##' p.prior <- BuildPrior(
##'   dists = rep("tnorm", npar),
##'   p1    = pop.mean,
##'   p2    = pop.scale*5,
##'   lower = c(0,-5, -5, 0, 0, 0, 0),
##'   upper = c(5, 7,  7, 1, 2, 1, 1))
##'
##' mu.prior <- ggdmc::BuildPrior(
##'   dists = rep("tnorm", npar),
##'   p1    = pop.mean,
##'   p2    = pop.scale*5,
##'   lower = c(0,-5, -5, 0, 0, 0, 0),
##'   upper = c(5, 7,  7, 1, 2, 1, 1)
##' )
##' sigma.prior <- BuildPrior(
##'   dists = rep("beta", npar),
##'   p1    = c(a=1, v.f1=1,v.f2 = 1, z=1, sz=1, sv=1, t0=1),
##'   p2    = rep(1, npar),
##'   upper = rep(2, npar))
##'
##' priors <- list(pprior=p.prior, location=mu.prior, scale=sigma.prior)
##'
##' ## Sampling
##' ## Processing time: 394.37 secs.
##' fit0 <- StartNewsamples(dmi, priors, thin = 2)
##' fit  <- run(fit0)
##' fit  <- run(fit, 1e2, add=TRUE)
##'
##' ## By default the type = 1 for location parameters
##' ## When recovery = TRUE, one must enter the true parameter to ps
##' est0 <- summary(fit, recovery = TRUE, ps = pop.mean, verbose = TRUE)
##' ## Explicitly enter type = 1
##' est0 <- summary(fit, recovery = TRUE, ps = pop.mean,  type=1, verbose = TRUE)
##' est0 <- summary(fit, recovery = TRUE, ps = pop.scale, type=2, verbose = TRUE)
##'
##' ## When recovery = FALSE (default), the function return parameter estimates
##' est0 <- summary(fit, verbose = TRUE, type=1)
##' est0 <- summary(fit, verbose = TRUE, type=2)
##'
##' ## To estimate individual participants, one must enter hyper = FALSE for a
##' ## hierarchical model fit
##' est0 <- summary(fit, hyper=FALSE, verbose = TRUE)
##' }
setGeneric("summary", function(object, ...) {
  standardGeneric("summary")
})

##' @rdname summary-methods
setMethod("summary", "posterior", function(object, start = 1,
                      end = NA, prob = c(0.025, 0.25, 0.5, 0.75, 0.975),
                      recovery = FALSE, ps = NA, verbose = FALSE,
                      digits = max(3, getOption("digits")-3)) {

  if (recovery)
  {
    message("Recovery summarises only default quantiles: ",
            appendLF = FALSE)
    prob <- c(0.025, 0.5, 0.975)
    cat(paste0(prob*100, "%"), "\n")

    if (any(is.na(ps))) stop("Some true values are NAs.")

    out <- summary_recoverone(object, start, end, ps, digits, prob, verbose)

  } else {
    if(verbose) {
        message("Summarises the following quantiles: ", appendLF = FALSE)
        cat(paste0(prob*100, "%"), "\n")
    }

    out <- summary_one(object, start, end, prob)

  }

  return(out)

})

##' @rdname summary-methods
setMethod("summary", "list", function(object, start = 1, end = NA,
                                      prob = c(0.025, 0.25, 0.5, 0.75, 0.975),
                                      recovery = FALSE, ps = NA,
                                      verbose = FALSE,
                                      digits = max(3, getOption("digits")-3))
{
  if (recovery) {
    message("Recovery shows only default quantiles: ", appendLF = FALSE)
    prob <- c(0.025, 0.25, 0.5, 0.75, 0.975)
    cat(paste0(prob*100, "%"), "\n")

    if (any(is.na(ps))) stop("Some true values are NAs.")
    if (!is.matrix(ps)) stop("Must provide true parameter matrix")
    out <- summary_recovermany(object, start, end, ps, digits, prob, verbose)
  } else {
    out <- summary_many(object, start, end, prob, digits, verbose)
  }
  return(out)
})

##' @rdname summary-methods
setMethod("summary", "hyper", function(object, hyper = TRUE, start = 1, end = NA,
                                       prob = c(0.025, 0.25, 0.5, 0.75, 0.975),
                                       recovery = FALSE, ps = NA, type = 1,
                                       verbose = FALSE,
                                       digits = max(3, getOption("digits")-3))
{
  if (recovery) {
    message("Recovery shows only default quantiles: ", appendLF = FALSE)
    prob <- c(0.025, 0.25, 0.5, 0.75, 0.975)
    cat(paste0(prob*100, "%"), "\n")

    if (hyper) {
      out <- summary_recoverhyper(object, start, end, ps, type, digits, prob,
                                  verbose)
    } else {
      if (any(is.na(ps))) stop("Some true values are NAs.")
      if (!is.matrix(ps)) stop("Must provide true parameter matrix")

      names(object@individuals) <- object@snames
      out <- summary_recovermany(object@individuals, start, end, ps, digits,
                                 prob, verbose)
    }
  } else {

    if (hyper) {
      out <- summary_hyper(object, start, end, type, prob, digits, verbose)
      if (verbose) {
        type_label <- ifelse(type==1, "location", "scale")
        message("Summarise ", type_label, " parameters:")
        print(round( out$quantiles, digits) )
      }

    } else {
      names(object@individuals) <- object@snames
      out <- summary_many(object@individuals, start, end, prob, digits, verbose)
    }
  }
  return(out)
})

### Simulate method ----------------
##' Get a n-cell matrix
##'
##' Constructs a matrix, showing how many responses to in each
##' cell. The function checks whether the format of \code{n} and \code{ns}
##' conform.
##'
##' \code{n} can be:
##' \enumerate{
##' \item an integer for a balanced design,
##' \item a matrix for an unbalanced design, where rows are subjects and
##' columns are cells. If the matrix is a row vector, all subjects
##' have the same \code{n} in each cell. If it is a column vector, all
##' cells have the same \code{n}. Otherwise each entry specifies the \code{n}
##' for a particular subject x cell combination. See below for concrete
##' examples.}
##'
##' @param ncell number of cells.
##' @param n number of trials.
##' @param ns number of subjects.
##' @examples
##' model <- BuildModel(
##'   p.map     = list(A = "1", B = "R", t0 = "1", mean_v = "M", sd_v = "M",
##'                   st0 = "1"),
##'   match.map = list(M = list(s1 = 1, s2 = 2)),
##'   constants = c(sd_v.false = 1, st0 = 0),
##'   factors   = list(S = c("s1","s2")),
##'   responses = c("r1", "r2"),
##'   type      = "norm")
##'
##' #######################30
##' ## Example 1
##' #######################30
##' cells <- as.numeric(sapply(model@factors, length))
##' ncell <- prod(cells)
##' GetNsim(ncell, ns = 2, n = 1)
##' #      [,1] [,2]
##' # [1,]    1    1
##' # [2,]    1    1
##'
##' #######################30
##' ## Example 2
##' #######################30
##' n <- matrix(c(1:2), ncol = 1)
##' #      [,1]
##' # [1,]    1  ## subject 1 has 1 response for each cell
##' # [2,]    2  ## subject 2 has 2 responses for each cell
##'
##' GetNsim(ncell, ns = 2, n = n)
##' #      [,1] [,2]
##' # [1,]    1    1
##' # [2,]    2    2
##'
##' #######################30
##' ## Example 3
##' #######################30
##' n <- matrix(c(1:2), nrow = 1)
##' #      [,1] [,2]
##' # [1,]    1    2
##' GetNsim(ncell, ns = 2, n = n)
##' #     [,1] [,2]
##' # [1,]   1    2 ## subject 1 has 1 response for cell 1 and 2 responses for cell 2
##' # [2,]   1    2 ## subject 2 has 1 response for cell 1 and 2 responses for cell 2
##'
##' #######################30
##' ## Example 4
##' #######################30
##' n <- matrix(c(1:4), nrow=2)
##' #      [,1] [,2]
##' # [1,]    1    3
##' # [2,]    2    4
##' GetNsim(ncell, ns = 2, n = n)
##' #      [,1] [,2]
##' # [1,]    1    3 ## subject 1 has 1 response for cell 1 and 3 responses for cell 2
##' # [2,]    2    4 ## subject 2 has 2 responses for cell 1 and 4 responses for cell 2
##
##' @export
GetNsim <- function (ncell, n, ns)
  ## Use only in simulate_many R function
{
  if (ns <= 1) stop("ns must be greater than or equal to 1.")
  if (is.vector(n) & (length(n) != 1))
  {
    if (!is.matrix(n)) stop("n must be a scalar, a vector, or a matrix")
  }

  ## Untested
  if (is.matrix(n))
  {
    dim1 <- dim(n)[1]
    dim2 <- dim(n)[2]
    if (dim1 == 1) {
      nmat <- matrix(rep(n, each = ns), ns) ## only cells differ
    } else if (dim2 == 1) {
      nmat <- matrix(rep.int(n, ncell), ns) ## only subjects differ
    } else {
      nmat <- n
    }
  } else {
    nmat <- matrix(rep.int(n, ns*ncell), ns)
  }

  dim1 <- dim(nmat)[1]
  dim2 <- dim(nmat)[2]

  if ( ns != dim1 )    stop(paste0("The n matrix must have ", ns, " rows"))
  if ( ncell != dim2 ) stop(paste0("The n matrix must have ", ncell, " columns"))
  return(nmat)
}

##' Constructs a ns x npar matrix,
##'
##' The matrix is used to simulate data. Each row represents one set of
##' parameters for a participant.
##'
##' One must enter either a vector or a matrix as true parameters
##' to the argument, \code{ps}, when presuming to simulate data based on a
##' fixed-effect model. When the assumption is to simulate data based on a
##' random-effect model, one must enter a prior object to the argument,
##' \code{prior} to first randomly generate a true parameter matrix.
##'
##' @param object a model object
##' @param nsub number of subjects.
##' @param prior a prior object
##' @param ps a vector or a matirx.
##' @param seed an integer specifying a random seed.
##' @return a ns x npar matrix
##' @examples
##' model <- BuildModel(
##' p.map     = list(a ="1", v = "1",z = "1", d = "1", sz = "1", sv = "1",
##'             t0 = "1", st0 = "1"),
##' match.map = list(M = list(s1 = "r1", s2 ="r2")),
##' factors   = list(S = c("s1", "s2")),
##' constants = c(st0 = 0, d = 0),
##' responses = c("r1", "r2"),
##' type      = "rd")
##'
##' p.prior <- BuildPrior(
##'   dists = c("tnorm", "tnorm", "beta", "beta", "tnorm", "beta"),
##'   p1    = c(a = 1, v = 0, z = 1, sz = 1, sv = 1, t0 = 1),
##'   p2    = c(a = 1, v = 2, z = 1, sz = 1, sv = 1, t0 = 1),
##'   lower = c(0, -5, NA, NA, 0, NA),
##'   upper = c(2,  5, NA, NA, 2, NA))
##'
##' ## Example 1: Randomly generate 2 sets of true parameters from
##' ## parameter priors (p.prior)
##' GetParameterMatrix(model, nsub=2, p.prior)
##' ##            a         v         z        sz       sv        t0
##' ## [1,] 1.963067  1.472940 0.9509158 0.5145047 1.344705 0.0850591
##' ## [2,] 1.512276 -1.995631 0.6981290 0.2626882 1.867853 0.1552828
##'
##' ## Example 2: Use a user-selected true parameters
##' true.vector  <- c(a=1, v=1, z=0.5, sz=0.2, sv=1, t0=.15)
##' GetParameterMatrix(model, nsub=2, ps = true.vector)
##' ##   a v   z  sz sv   t0
##' ## 1 1 1 0.5 0.2  1 0.15
##' ## 2 1 1 0.5 0.2  1 0.15
##'
##' ## Example 3: When a user enter arbritary sequence of parameters.
##' ## Note sv is before sz. It should be sz before sv
##' ## See correct sequence, by entering "model@pnames"
##' ## GetParameterMatrix will rearrange the sequence.
##' true.vector  <- c(t0=15, a=1, v=1, z=0.5, sv=1, sz = .2)
##' GetParameterMatrix(model, nsub=2, ps= true.vector)
##' ##   a v   z  sz sv   t0
##' ## 1 1 1 0.5 0.2  1 0.15
##' ## 2 1 1 0.5 0.2  1 0.15
##'
##' @export
GetParameterMatrix <- function(object, nsub, prior, ps, seed = NULL)
{
  ## Used in simulate_one and simulate_many
  pnames <- object@pnames
  if (missingArg(prior) && missingArg(ps)) stop("Neither prior nor p.vector was found")

  ## use ps; fixed-effect model
  if ( missingArg(prior) )
  {
    if (is.vector(ps)) {
      if (check_pvec(ps, object)) stop("ps is not incompatible with model")
      ps    <- ps[pnames]
      psmat <- matrix(rep(ps, each = nsub), nrow=nsub, dimnames = list(NULL, pnames))
    } else if (is.matrix(ps)) {
      ps    <- ps[, pnames]
      psmat <- matrix(ps, nrow=nsub, dimnames = list(NULL, pnames))
    } else {
      stop("ps must be a vector or a matrix")
    }
    rownames(psmat) <- 1:nsub
  }
  else  ## use prior; random-effect model
  {
    if ( !all( pnames %in% prior@pnames) ) stop("priors must match the model")
    set.seed(seed)
    psmat <- rprior(prior, nsub)
  }

  return(psmat)
}


##' Simulate Choice Responses
##'
##' The function is an extension of the simulate function in \code{stats}
##' pacakge. It simulates the data from either two-alternative force choice
##' tasks, multiple-alternative force choice task, or continuous report tasks.
##'
##' The function simulates data either for one participant or multiple
##' participants. The simulation process is based on the model object, entering
##' via \code{object} argument. For simulating one participant, one must supply
##' a true parameter vector to the \code{ps} argument.
##'
##' For simulating multiple participants, one can enter a matrix or a row
##' vector as true parameters. Each row is used to generate the data for a
##' participant. This process is usually dubbed the fixed-effect modelling.
##' To generate data via the random-effect modelling, one must supply a set of
##' prior distributions. In this case, \code{ps} argument is unused. Note in
##' some cases, a random-effect modelling may fail to draw data from the model,
##' because true parameters are randomly drawn
##' from prior distributions. This would happen sometimes for example in the
##' diffusion decision model, because certain parameter combinations are
##' considered invalid (e.g., t0 < 0, zr > a) for obvious reasons.
##'
##' \code{ps} can be a row vector, in which case each participant has one set
##' of identical parameters. It can also be a matrix with one row per
##' participant, in which case it must have \code{ns} rows. The true values will
##' be saved as \code{parameters} attribute in the output.
##'
##' @param object a model object.
##' @param nsim number of observations. \code{nsim} can be a single number
##' for a balanced design or a matrix for an unbalanced design, where rows
##' are participants and columns are design cells. If the matrix has one row
##' than all participants have the same \code{nsim} in each cell, if it has one
##' column then all cells have the same \code{nsim}; Otherwise each entry
##' specifies the \code{nsim} for a particular participant x design cell
##' combination.
##' @param nsub number of participants
##' @param prior a prior object
##' @param ps a true parameter vector or matrix.
##' @param seed a user specified random seed.
##' @param ... additional optional arguments.
##' @return a data frame
##' @examples
##' model <- BuildModel(
##'   p.map     = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1",
##'                    st0 = "1"),
##'   match.map = list(M = list(s1 = 1, s2 = 2)),
##'   factors   = list(S = c("s1", "s2")),
##'   constants = c(st0 = 0, sd_v = 1),
##'   responses = c("r1", "r2"),
##'   type      = "norm")
##'
##' p.vector <- c(A = .75, B = 1.25, t0 = .15, mean_v.true = 2.5,
##'               mean_v.false = 1.5)
##' ntrial <- 100
##' dat <- simulate(model, nsim = ntrial, ps = p.vector)
##' @export
##' @docType methods
##' @rdname simulate-methods
setMethod("simulate", "model", function(object, nsim = 1, seed = NULL,
                                        nsub, prior = NA, ps = NA)
{
  if ( missingArg(nsub) ) {

    if (anyNA(ps)) { ps <- GetParameterMatrix(object, 1, prior, seed) }
    out <- simulate_one(object, nsim, ps, seed)
    attr(out, "parameters") <- ps

  } else {
    if (anyNA(prior@priors) & anyNA(ps)) {
      stop("Must supply either sampling distribution or a true vector.")
    }
    out <- simulate_many(object, nsim, nsub, prior, ps, seed)
  }

  return(out)
})

simulate_one <- function(model, n, ps, seed)
{
  if (check_pvec(ps, model)) stop("p.vector and model incompatible")
  resp <- model@responses
  type <- model@type
  levs <- model@factors
  facs <- createfacsdf(model)
  nvec <- check_n(n, facs)
  dat  <- nadf(nvec, facs, levs, type)
  row1 <- 1

  dfnames <- names(dat)
  reserved_names <- c("RT", "R")
  if (type == "cddm") reserved_names <- c("R", "RT", "A")
  for (i in 1:nrow(facs))
  {
    ## simulate use n1.order == FALSE for LBA
    pmat <- TableParameters(ps, i, model, FALSE)
    if (nvec[i] == 0) next
    rown <- row1 + nvec[i] - 1
    dat[row1:rown, reserved_names] <- random(type, pmat, nvec[i], seed)
    row1 <- rown+1
  }

  if (type %in% c("norm", "norm_pda", "norm_pda_gpu", "rd")) {
    dat$R <- factor(dat$R, levels = 1:length(resp), labels = resp)
    if (type == "rd") dat <- FlipResponse_rd(model, dat, facs)
  } else if (type == "cddm") {
    dat$R <- LabelTheta(dat, resp)
  }

  return(dat)
}

simulate_many <- function(object, n, ns, prior, ps, seed)
{
  # model
  # n = 30
  # ns = 8
  # prior = pop.prior

  ## sapply(facs, length) return the numbers of level in each factor in a named
  ## vector. Then prod multiplies them to get the number of cell
  cells <- as.numeric(sapply(object@factors, length))
  ncell <- prod(cells)

  n  <- GetNsim(ncell, n, ns)
  ps <- GetParameterMatrix(object, ns, prior, ps, seed)

  multimodel <- is_multimodel(object, ns)
  if (multimodel) modeli <- object[[1]] else modeli <- object

  ndatai <- cumsum(c(0, matrixStats::rowSums2(n))); ## index boundaries
  datr <- (ndatai[1] + 1):(ndatai[2]); ## First subject's trial index

  ## Simulate first subject; modeli must be 'model' class
  dat <- cbind(s = rep.int(1, length(datr)),
               simulate_one(modeli, n[1,], ps[1,], seed))

  # i <- 2
  if (ns > 1)
  {
    for (i in 2:ns)
    {
      if (multimodel) modeli <- object[[i]]
      datr <- (ndatai[i] + 1):(ndatai[i + 1]) ## continue to index trials
      dat  <- rbind(dat,
                    cbind(s = rep.int(i, length(datr)),
                          simulate_one(modeli, n[i,], ps[i,], seed)))
    }
  }

  dat$s <- factor(dat$s)
  attr(dat, "parameters") <- ps
  ## if ps is not auto-created by p.prior, save the user's p.prior in 'attribute'

  if(!anyNA(prior@priors)) attr(dat, "priors") <- prior@priors
  return(dat)
}

