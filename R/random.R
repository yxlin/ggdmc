######### r-functions -----------------------------------
rdiffusion <- function (n,
                        a, v, z = 0.5*a, d = 0, sz = 0, sv = 0, t0 = 0, st0 = 0,
                        s = 1, precision = 3, stop_on_error = TRUE)
{
  ## @author Underlying C code by Jochen Voss and Andreas Voss. Porting and R
  ## wrapping by Matthew Gretton, Andrew Heathcote, Scott Brown, and Henrik
  ## Singmann.
  ## \code{qdiffusion} by Henrik Singmann. This function is extracted from
  ## rtdists written by the above authors.

  if(any(missing(a), missing(v), missing(t0)))
    stop("a, v, and/or t0 must be supplied")

  s <- rep(s, length.out = n)
  a <- rep(a, length.out = n)
  v <- rep(v, length.out = n)
  z <- rep(z, length.out = n)
  z <- z/a  # transform z from absolute to relative scale (which is currently required by the C code)
  d <- rep(d, length.out = n)
  sz <- rep(sz, length.out = n)
  sz <- sz/a # transform sz from absolute to relative scale (which is currently required by the C code)
  sv <- rep(sv, length.out = n)
  t0 <- rep(t0, length.out = n)
  st0 <- rep(st0, length.out = n)
  t0 <- recalc_t0 (t0, st0)

  # Build parameter matrix (and divide a, v, and sv, by s)
  params <- cbind (a/s, v/s, z, d, sz, sv/s, t0, st0)

  # Check for illegal parameter values
  if(ncol(params)<8)
    stop("Not enough parameters supplied: probable attempt to pass NULL values?")
  if(!is.numeric(params))
    stop("Parameters need to be numeric.")
  if (any(is.na(params)) || !all(is.finite(params)))
    stop("Parameters need to be numeric and finite.")

  randRTs    <- vector("numeric",length=n)
  randBounds <- vector("numeric",length=n)

  #uniques <- unique(params)
  parameter_char <- apply(params, 1, paste0, collapse = "\t")
  parameter_factor <- factor(parameter_char, levels = unique(parameter_char))
  parameter_indices <- split(seq_len(n), f = parameter_factor)

  for (i in seq_len(length(parameter_indices)))
  {
    ok_rows <- parameter_indices[[i]]

    # Calculate n for this row
    current_n <- length(ok_rows)

    out <- r_fastdm (current_n,
                     params[ok_rows[1],1:8],
                     precision,
                     stop_on_error=stop_on_error)
    #current_n, uniques[i,1:8], precision, stop_on_error=stop_on_error)

    randRTs[ok_rows]    <- out$rt
    randBounds[ok_rows] <- out$boundary
  }
  response <- factor(randBounds, levels = 0:1, labels = c("lower", "upper"))
  data.frame(rt = randRTs, response)
}

######### Generic functions -----------------------------------
##' Generate random numbers
##'
##' A wrapper function for generating random numbers of either
##' the model type, \code{rd}, or \code{norm}.
##'
##' @param type a character string of the model type
##' @param pmat a matrix of response x parameter
##' @param n number of observations
##' @param seed an integer specifying a random seed
##' @export
random <- function(type, pmat, n, seed = NULL)
{

  set.seed(seed)

  if (type == "rd") {

    out <- rdiffusion(n, a = pmat$a[1], v = pmat$v[1],
      t0 = pmat$t0[1],
      z  = pmat$z[1]*pmat$a[1], # convert to absolute
      d  = pmat$d[1],
      sz = pmat$sz[1]*pmat$a[1],
      sv = pmat$sv[1], st0 = pmat$st0[1], stop_on_error = TRUE)

  } else if (type %in% c("norm", "norm_pda", "norm_pda_gpu")) {
    ## pmat: A b t0 mean_v sd_v st0

    out <- rlba_norm(n, pmat[, 1], pmat[, 2], pmat[, 4], pmat[, 5],
      pmat[,3], pmat[1,6], TRUE) ## posdrift is RUE

  } else {
    stop("Model type yet created")
  }

  attr(out, "seed") <- seed
  return(out)
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
##' @param x a model object
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
##' GetParameterMatrix(model, 2, p.prior)
##' ##            a         v         z        sz       sv        t0
##' ## [1,] 1.963067  1.472940 0.9509158 0.5145047 1.344705 0.0850591
##' ## [2,] 1.512276 -1.995631 0.6981290 0.2626882 1.867853 0.1552828
##'
##' ## Example 2: Use a user-selected true parameters
##' true.vector  <- c(a=1, v=1, z=0.5, sz=0.2, sv=1, t0=.15)
##' GetParameterMatrix(model, 2, NA, true.vector)
##' ##      a v   z  sz sv   t0
##' ## [1,] 1 1 0.5 0.2  1 0.15
##' ## [2,] 1 1 0.5 0.2  1 0.15
##' GetParameterMatrix(model, 2, ps = true.vector)
##'
##' ## Example 3: When a user enter arbritary sequence of parameters.
##' ## Note sv is before sz. It should be sz before sv
##' ## See correct sequence, by entering "attr(model, 'p.vector')"
##' ## GetParameterMatrix will rearrange the sequence.
##' true.vector  <- c(a=1, v=1, z=0.5, sv=1, sz = .2, t0=.15)
##' GetParameterMatrix(model, 2, NA, true.vector)
##' ##      a v   z  sz sv   t0
##' ## [1,] 1 1 0.5 0.2  1 0.15
##' ## [2,] 1 1 0.5 0.2  1 0.15
##'
##' @export
GetParameterMatrix <- function(x, nsub, prior = NA, ps = NA, seed = NULL)
## Used in simulate.model and simulate_many
{
  message1 <- "Parameters are incompatible with model"
  pnames <- names(attr(x, "p.vector"))

  ## use ps; fixed-effect model
  if (anyNA(prior))
  {
    if (is.vector(ps)) {
      if (check_pvec(ps, x)) stop(message1)
      ps    <- ps[pnames]
      pss   <- rep(ps, each = nsub)
      psmat <- matrix(pss, nsub, dimnames = list(NULL, pnames))
    } else if (is.matrix(ps)) {
      ps    <- ps[, pnames]
      psmat <- matrix(ps, nsub, dimnames = list(NULL, pnames))
    } else {
      if ((nsub != dim(ps)[1])) stop("ps matrix must have nsub rows")
      if (check_pvec(ps[1,], x)) stop(message1)
    }

    rownames(psmat) <- 1:nsub
  }
  else  ## use prior; random-effect model
  {
    if (!all( pnames %in% names(prior))) stop(message1)
    set.seed(seed)
    psmat <- rprior(prior[pnames], nsub)

  }

  return(psmat)
}


simulate_one <- function(model, n, ps, seed)
{
  if (check_pvec(ps, model)) stop("p.vector and model incompatible")
  resp <- attr(model, "responses")
  type <- attr(model, "type")
  levs <- attr(model, "factors")
  facs <- createfacsdf(model)
  nvec <- check_n(n, facs)
  dat  <- nadf(nvec, facs, levs, type)
  row1 <- 1

  dfnames <- names(dat)
  Xnames   <- dfnames[!dfnames %in% c("R", "N", "Y")]

  # i <- 2
  for (i in 1:nrow(facs)) {
    pmat <- TableParameters(ps, i, model, FALSE) ## simulate use n1.order == FALSE
    rown <- row1 + nvec[i] - 1

    dat[row1:rown, c("RT", "R")] <- random(type, pmat, nvec[i], seed)
    row1 <- rown+1
  }

  dat$R <- factor(dat$R, levels = 1:length(resp), labels = resp)
  if (type == "rd") dat <- FlipResponse_rd(model, dat, facs)
  return(dat)
}

simulate_many <- function(model, n, ns, prior, ps, seed)
{

  n  <- GetNsim(model, n, ns)
  ps <- GetParameterMatrix(model, ns, prior, ps, seed)

  ismanysub <- ismanymodels(model, ns)
  if(ismanysub) modeli <- model[[1]] else modeli <- model

  ndatai <- cumsum(c(0, matrixStats::rowSums2(n))); ## index boundaries
  datr <- (ndatai[1] + 1):(ndatai[2]); ## First subject's trial index

  ## Simulate first subject; modeli must be 'model' class
  dat <- cbind(s = rep.int(1, length(datr)),
    simulate_one(modeli, n[1,], ps[1,], seed))

  if (ns > 1) {
    for (i in 2:ns) {
      if (ismanysub) modeli <- model[[i]] else modeli <- model
      datr <- (ndatai[i] + 1):(ndatai[i + 1]) ## continue to index trials
      dat  <- rbind(dat,
        cbind(s = rep.int(i, length(datr)),
          simulate_one(modeli, n[i,], ps[i,], seed)))
    }
  }

  dat$s <- factor(dat$s)
  attr(dat, "parameters") <- ps
  ## if ps is not auto-created by p.prior, save the user's p.prior in 'attribute'
  if(!anyNA(prior)) attr(dat, "p.prior") <- prior
  return(dat)
}

##' Simulate response time data
##'
##' Simulate response time data either for one subject or multiple subjects.
##' The simulation is based on a model object. For one subject, one must supply
##' a true parameter vector to the \code{ps} argument.
##'
##' For multiple subjects, one can enter a matrix (or a row vector) as true
##' parameters. Each row is to generate data separately for a subject.  This is
##' the fixed-effect model. To generate data based on a random-effect
##' model, one must supply a prior object.  In this case, \code{ps} argument
##' is unused. Note in some cases, a random-effect model may fail to draw data
##' from the model, because true parameters are randomly drawn from
##' a prior object.  This would happen sometimes in diffusion model, because
##' certain parameter combinations are considered invalid.
##'
##' \code{ps} can be a row vector, in which case each subject has identical
##' parameters.  It can also be a matrix with one row per subject, in which
##' case it must have \code{ns} rows. The true values will be saved as
##' \code{parameters} attribute in the output object.
##'
##' @param object a model object.
##' @param nsim number of trials / responses. \code{n} can be a single number
##' for a balanced design or a matrix for an unbalanced design, where rows
##' are subjects and columns are design cells. If the matrix has one row then
##' all subjects have the same \code{n} in each cell, if it has one column then
##' all cells have the same \code{n}; Otherwise each entry specifies the
##' \code{n} for a particular subject x design cell combination.
##' @param nsub number of subjects
##' @param prior a prior object
##' @param ps a true parameter vector or matrix.
##' @param seed a user specified random seed.
##' @param ... additional optional arguments.
##' @return a data frame
##' @importFrom stats simulate
##' @export
simulate.model <- function(object, nsim = NA, seed = NULL, nsub = NA,
  prior = NA, ps = NA, ...)
{
  if (is.na(nsub)) {
    if (is.na(nsim)) stop("How many response you want to generate? Must supply n")
    if (anyNA(ps)) {
      ps <- GetParameterMatrix(object, 1, prior, seed)
      ## stop("Some true parameters missing")
    }
    out <- simulate_one(object, nsim, ps, seed)
    attr(out, "parameters") <- ps
  } else {
    message1 <- "Must supply either sampling distribution or a true vector."
    if (anyNA(prior) & anyNA(ps)) stop(message1)
    out <- simulate_many(object, nsim, nsub, prior, ps, seed)
  }
  return(out)
}

#### Post-predictive functions --------------
# predict_one <- function(object, npost = 100, rand = TRUE, factors = NA,
#                         xlim = NA, seed = NULL)
# {
#   # object <- fit
#   # factors = NA
#   model <- attributes(object$data)$model
#   facs <- names(attr(model, "factors"))
#   class(object$data) <- c("data.frame", "list")
#
#   if (!is.null(factors))
#   {
#     if (any(is.na(factors))) factors <- facs
#     if (!all(factors %in% facs))
#       stop(paste("Factors argument must contain one or more of:",
#                  paste(facs, collapse=",")))
#   }
#
#   resp <- names(attr(model, "responses"))
#   ns   <- table(object$data[,facs], dnn = facs)
#   npar   <- object$n.pars
#   nchain <- object$n.chains
#   nmc    <- object$nmc
#   ntsample <- nchain * nmc
#   pnames   <- object$p.names
#   # str(object$theta) ## npar x nchain x nmc
#   # str(thetas)       ## (nchain x nmc) x npar
#   thetas <- matrix(aperm(object$theta, c(3,2,1)), ncol = npar)
#
#   # head(thetas)
#   # head(object$theta[,,1:2])
#
#   colnames(thetas) <- pnames
#
#   if (is.na(npost)) {
#     use <- 1:ntsample
#   } else {
#     if (rand) {
#       use <- sample(1:ntsample, npost, replace = F)
#     } else {
#       use <- round(seq(1, ntsample, length.out = npost))
#     }
#   }
#
#   npost  <- length(use)
#   posts   <- thetas[use, ]
#   nttrial <- sum(ns) ## number of total trials
#
#   ## should replace with parallel
#   v <- lapply(1:npost, function(i) {
#     simulate_one(model, n = ns, ps = posts[i,], seed = seed)
#   })
#   out <- data.table::rbindlist(v)
#   # names(out) <- names(object$data)
#   reps <- rep(1:npost, each = nttrial)
#   out <- cbind(reps, out)
#
#   if (!any(is.na(xlim)))
#   {
#     out <- out[RT > xlim[1] & RT < xlim[2]]
#   }
#
#   attr(out, "data") <- object$data
#   return(out)
# }

######### Utility functions -----------------------------------
## [MG 20150616]
## In line with LBA, adjust t0 to be the lower bound of the non-decision time distribution
## rather than the average
## Called from prd, drd, rrd (Extracted from rtdists)
recalc_t0 <- function (t0, st0) { t0 <- t0 + st0/2 }

"%w/o%" <- function(x, y) x[!x %in% y] #--  x without y


