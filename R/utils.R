### Hierarchical tools ------------------------------------------------
##' Get n-cell matrix
##'
##' Constructs a matrix, showing how many responses to in each
##' cell. It also checks if the format of \code{n} and \code{ns}
##' conform to the standard.
##'
##' \code{n} can be:
##' \enumerate{
##' \item one integer for a balanced design,
##' \item a matrix for an unbalanced design, where rows are subjects and
##' columns are cells. If the matrix is a row vector, all subjects
##' have the same \code{n} in each cell. If it is a column vector, all
##' cells have the same \code{n}. Otherwise each entry specifies the \code{n}
##' for a particular subject x design cell combination. See below
##' for concrete examples.
##' }
##'
##' @param model a model object
##' @param n number of trials / responses.
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
##' GetNsim(model, ns = 2, n = 1)
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
##' GetNsim(model, ns = 2, n = n)
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
##' GetNsim(model, ns = 2, n = n)
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
##' ggdmc::GetNsim(model, ns = 2, n = n)
##' #      [,1] [,2]
##' # [1,]    1    3 ## subject 1 has 1 response for cell 1 and 3 responses for cell 2
##' # [2,]    2    4 ## subject 2 has 2 responses for cell 1 and 4 responses for cell 2
##
##' @export
GetNsim <- function(model, n, ns) {
  if (ns <= 1) stop("Use simulate.model instead to simulate one participant.")
  if (is.vector(n) & (length(n) != 1)) {
    if (!is.matrix(n)) stop("n must be a scalar, a vector, or a matrix")
  }

  facs  <- attr(model, "factors")
  ## lapply(facs, length) return the numbers of level in each factor
  ## unlist and prod then multiply them to get the number of cell
  ncell <- prod(unlist(lapply(facs, length)))

  if (is.matrix(n)) {     ## Untested
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

  dim1 <- dim(nmat)[1]   ## n has been altered in ifelse
  dim2 <- dim(nmat)[2]

  if ( ns != dim1 ) stop(paste0("The n matrix must have ", ns, " rows"))
  if ( ncell != dim2 ) stop(paste0("The n matrix must have ", ncell, " columns"))
  return(nmat)
}

##' Test whether `model` object has many models
##'
##' @param model a model object
##' @param ns number of subjects.
##'
##' @export
ismanymodels <- function(model, ns = NA) {
  if (is.list(model)) {
    if (length(model) != ns) stop("number of participants not equal to number of models")
    out <- TRUE
  } else {
    if (is.na(ns)) stop("Please indicate the number of participants")
    out <- FALSE
  }
  return(out)
}

### PDA ------------------------------------------------
checkforce <- function(force = FALSE, nsamp = NA) {
  ## This is a tentative function
  if( is.na(nsamp) ) stop("number of MC samples not found")
  if (is.numeric(force)) force <- c(rep(FALSE, force - 1), TRUE)
  force <- c(TRUE, rep(force, length.out = nsamp - 1))

  if (!is.logical(force) || (length(force) != nsamp))
    stop(paste("force argument must be a logical vector of length", nsamp - 1))
  return(force)
}

##' Create a vector for re-calculation
##'
##' The is a PDA function. It creates an index vector, indicating which
##' iteration to recalculate likelihoods.
##'
##' @param x posterior samples
##' @param nth every nth step to recalculate
##' @examples
##' \dontrun{
##'   MakeForce(hsam[[1]], 3)
##' }
##'
##' @export
MakeForce <- function(x, nth) {
  nsamp <- 1 + (x$nmc - x$start) * x$thin
  if (!nth) {
    out <- rep.int(0, nsamp)
  } else {
    out <- rep.int(0, nsamp)
    if (!is.numeric(nth)) stop("nth must be an integer.")
    for(i in 1:nsamp) { if (i %% nth == 0) {out[i] <- 1} }
  }
  return(out)
}


### Sampling ------------------------------------------------
CheckHyperDMI <- function(data = NULL, nchain = NULL) {
  if (is.null(data)) stop("No data")
  if (!is.list(data)) stop ("data must be a list")
  if (is.data.frame(data)) stop("data is a list with each if its elements is data.frame")
  model1 <- attr(data[[1]], "model")
  pnames <- GetPNames(model1)
  if (is.null(nchain)) {
    nchain <- 3*length(pnames)
    message("nchain is not supplied. Use default ", nchain, " chains")
    ## message("If you are using DGMC, this nchain may be inappropriate.")
  }
  return(nchain)
}

##' Check Data Model Instance
##'
##' Return a model object extracted either from a data model instance or
##' an object storing posterior samples. The function checks also
##' \enumerate{
##' \item When x stores DMI of one participant, if DMI is a \code{data.frame},
##' \item When x is an object of posterior samples, if x is a list of many subjects,
##' \item whether model is successfully created,
##' \item whether prior is suppled or we can extract it from 'samples'
##' }
##'
##' @param x a data-model instance
##' @param prior a parameter prior list
##' @param theta1 a user-supplied theta cube
##' @param nchain number of MCMC chains
##' @export
CheckDMI <- function(x = NULL, prior = NULL, theta1 = NULL, nchain = NULL) {
  ## x == data == DMI

  if (!is.null(x) && !is.data.frame(x)) stop("Data must be a data frame")
  if (is.null(x)) {
    stop("No data model instance")
  } else {
    model <- attr(x, "model")
  }

  npar <- length(GetPNames(model))
  if (is.null(nchain)) nchain <- 3*npar
  if (is.null(model)) stop("Must specify a model")
  if (is.null(prior)) stop("Must specify a prior argument")
  if (!is.null(theta1) && !is.matrix(theta1) || (!all(dim(theta1)==c(nchain, npar))))
    stop("theta1 must be a nchain x npar matrix")
  return(model)
}

CheckSamples <- function(samples = NULL, prior = NULL, theta1 = NULL) {
  ## x == samples

  if (is.null(samples)) stop("Must supply samples") else data <- samples$data
  if (!is.null(data) && !is.data.frame(data)) stop("Data must be a data frame")
  if (!is.null(samples) && is.null(samples$theta)) stop("Use StartHypersamples")
  if (is.null(data)) model <- attr(samples$data, "model") else model <- attr(data, "model")
  npar <- length(GetPNames(model))
  nchain <- samples$n.chain

  if (is.null(model)) stop("Must specify a model")
  if (is.null(prior) && is.null(samples)) stop("Must specify a p.prior argument")
  if (!is.null(theta1) && !is.matrix(theta1) || (!all(dim(theta1)==c(nchain, npar))))
    stop("theta1 must be a nchain x npar matrix")
  return(model)
}

checkblocks <- function(blocks, hyper) {
  ## This function is tentative
  if ( any(is.na(blocks)) ) {
    blocks <- as.list(1:hyper$n.pars)
  } else { # check
    if (any(unlist(lapply(blocks, function(x) {
      length(x)==1 || all(hyper$has.sigma[x][1] == hyper$has.sigma[x][-1])
    }))))
      stop("Cant mix hyper-paramaters with and without sigma in a block")
  }
  return(blocks)
}

##' Check rejection rate
##'
##' Check rejection rate for posterior samples
##'
##' @param object posterior samples
##' @param verbose print more information
##' @export
CheckRJ <- function(object, verbose = TRUE) {
  nchain <- object$n.chains
  # nsamp <- 1 + (object$nmc - object$start) * object$thin;
  nsamp <- object$nmc

  for(i in 1:nchain) {
    mr <- sum(object$rejection_rate[i, 2:nsamp] == 2 |
        object$rejection_rate[i, 2:nsamp] == 1) / nsamp
    rj <- sum(object$rejection_rate[i, 2:nsamp] == 2 |
        object$rejection_rate[i, 2:nsamp] == 4) / nsamp
    if (verbose) {
      cat("Chain ", i)
      cat(": migration and rejection rates: ", round(mr, 2), " ",
        round(rj, 2), "\n")
      # cat(": rejection rates: ", round(rj, 2), "\n")
    }
  }
}

### MCMC ------------------------------------------------
##' Create a MCMC list
##'
##' @param x posterior samples
##' @param start start from which iteration
##' @param end end at which iteration
##' @param pll a Boolean switch for calculating posterior log-likelihood
##' @export
mcmc_list.model <- function(x, start = 1, end = NA, pll = TRUE) {

  if (is.null(x$theta)) stop("Use hyper mcmc_list")
  if ( is.na(end) ) end <- x$nmc
  if ( end <= start ) stop("End must be greater than start")

  if (pll) {
    message("Convert posterior log-likelihood")
    lp <- x$summed_log_prior[start:end, ] + x$log_likelihoods[start:end, ]
    colnames(lp) <- 1:ncol(lp)
    out <- coda::mcmc.list(lapply(data.frame(lp), function(x) mcmc(x)))
  } else {
    message("Convert marginal log-likelihood")
    lp <- x$log_likelihoods[start:end, ]
    out <- theta2mcmclist(x, start = start, end = end, thin = 1)

  }
  return(out)
}

##' Prepare posterior samples for plotting version 1
##'
##' Convert MCMC chains to a data frame for plotting
##'
##' @param x posterior samples
##' @param start which iteration to start
##' @param end end at which iteration
##' @param pll a Boolean switch to make posterior log likelihood
##' @export
ConvertChains <- function(x, start = 1, end = NA, pll = TRUE) {

  if (x$n.chains == 1) stop ("MCMC needs multiple chains to check convergence")
  if (is.null(x$theta)) stop("Use hyper mcmc_list")
  if ( is.na(end) ) end <- x$nmc
  if ( end <= start ) stop("End must be greater than start")

  nchain <- x$n.chain
  npar   <- x$n.pars
  pnames <- x$p.names
  iter   <- start:end
  mcmclist <- mcmc_list.model(x, start, end, pll)

  v <- lapply(seq_along(mcmclist), function(k) {
    tmp1 <- sapply(mcmclist[k][[1]], c)
    if (pll) {
      d <- data.frame(Iteration = iter,
        Parameter = "lp", value = tmp1)
    } else {
      d <- data.frame(Iteration = rep(iter, npar),
        Parameter = rep(pnames, each = length(iter)), value = tmp1)
    }

    d$Chain <- k
    d[, c("Iteration", "Chain", "Parameter", "value")]
  })


  D <- data.table::rbindlist(v)
  D$Parameter <- factor(D$Parameter)
  D$Chain <- factor(D$Chain)
  attr(D, "nThin")       <- x$thin
  attr(D, "nIterations") <- end
  attr(D, "nChains")     <- nchain
  attr(D, "nParameters") <- npar

  return(D)
}


##' Prepare posterior samples for plotting version 2
##'
##' Convert MCMC chains to a data frame for plotting
##'
##' @param x posterior samples
##' @param pll a Boolean switch to make posterior log likelihood
##' @export
ConvertChains2 <- function(x, pll) {

  nchain <- attr(x, "nchain")
  npar <- attr(x, "npar")
  pnames <- attr(x, "pnames")
  start <- attr(x, "start")
  end <- attr(x, "end")

  if (nchain == 1) stop ("MCMC needs multiple chains to check convergence")
  # if (is.null(x$theta)) stop("Use hyper mcmc_list")
  if ( is.na(end) ) end <- attr(x, "nmc")
  if ( end <= start ) stop("End must be greater than start")
  iter   <- start:end
  # mcmclist <- mcmc_list.model(x, start, end, pll)

  v <- lapply(seq_along(x), function(k) {
    tmp1 <- sapply(x[k][[1]], c)
    if (pll) {
      d <- data.frame(Iteration = iter,
        Parameter = "lp", value = tmp1)
    } else {
      d <- data.frame(Iteration = rep(iter, npar),
        Parameter = rep(pnames, each = length(iter)), value = tmp1)
    }

    d$Chain <- k
    d[, c("Iteration", "Chain", "Parameter", "value")]
  })


  D <- data.table::rbindlist(v)
  D$Parameter <- factor(D$Parameter)
  D$Chain <- factor(D$Chain)
  attr(D, "nThin")       <- attr(x, "thin")
  attr(D, "nIterations") <- end
  attr(D, "nChains")     <- nchain
  attr(D, "nParameters") <- npar

  return(D)
}



### Bayes --------------------------------------------------------------------
##' Add log-likelihoods across subjects at the hyper level
##'
##' An external sum log likelihoods for hyper parameters.
##'
##' @param ps a nsubject x npar matrix
##' @param pp a temporary pp.prior values extracted from phi. The temporary
##' pp.prior has no attached parameter names. phi is two-element list. First
##' element is a location array of nchain x npar x nmc; second element is a
##' scale array of nchain x npar x nmc.
##' @param prior prior distributions at data level
##'
##' @export
hsumloglike <- function(ps, pp, prior) {
  term1 <- apply(ps, 1, sumlogpriorNV, assign_pp(pp, prior))
  return(apply(term1, 2, sum))
}

##' Slot pp values into p.prior
##'
##' Assign hyper-prior distributions to prior distributions
##'
##' @param pp pp.prior
##' @param prior p.prior
##' @export
assign_pp <- function(pp, prior) {
  for (i in 1:length(prior)) prior[[i]][1:2] <- c(pp[[1]][i], pp[[2]][i])
  return(prior)
}


### Old DMC tools ------------------------------------------------
##' Censor missing values and RT outliers
##'
##' \code{censor} requests a data frame minimally with three columns, indicating
##' stimulus (S factor), response (R factor) and response time
##' (RT).
##'
##' @param x a data frame
##' @param xlim the lower and upper censoring boundaries
##' @examples
##' model <- BuildModel(
##'         p.map     = list(a = "1", v = "1", z = "1", d = "1", t0 = "1",
##'                     sv = "1", sz = "1", st0 = "1"),
##'         constants = c(st0 = 0, d = 0, sz = 0, sv = 0),
##'         match.map = list(M = list(s1 = "r1", s2 = "r2")),
##'         factors   = list(S = c("s1", "s2")),
##'         responses = c("r1", "r2"),
##'         type      = "rd")
##' p.vector <- c(a = 1, v = 1, z = 0.5, t0 = .15)
##' dat <- simulate(model, 100, ps = p.vector)
##' dmi <- BuildDMI(dat, model)
##'
##' ## Trim off RTs below .05 and above 10 s
##' censored_data <- censor(dat, xlim = c(.05, 10))
##' @export
censor <- function(x, xlim = c(0, Inf)) {
  if (sum( names(x) %in% c("S","R","RT") ) < 3)
    stop("The data.frame must has, at least, three columns, S, R and RT.")
  if (!is.factor(x$R)) x$R <- factor(x$R) ## Make R column as.factor
  ## p.na <- mean(is.na(x$RT))   ## proportion of na response
  is.in <- !is.na(x$RT)       ## TRUE/FALSE index which row is not NA
  is.in[is.in] <- x$RT[is.in] > xlim[1] & x$RT[is.in] < xlim[2]
  return(x[is.in,])
}


##' Convert factor levels to a data frame
##'
##' \code{fac2df} takes a model object created by BuildModel and returns a
##' data frame with all combination of factor levels.
##'
##' @param x a model object
##' @return a data frame
##' @examples
##' model <- BuildModel(
##'  p.map     = list(a = "1", v = "1", z = "1", d = "1", sz = "1", sv = "1",
##'                   t0 = "1", st0 = "1"),
##'  constants = c(st0 = 0, d = 0),
##'  match.map = list(M = list(s1 = "r1",s2 = "r2")),
##'  factors   = list(S = c("s1","s2")),
##'  responses = c("r1", "r2"),
##'  type      = "rd")
##'
##' df1 <- fac2df(model)
##'
##'
##' model <- BuildModel(
##'             p.map     = list(A = "1", B = "1", v = "M", sv = "M", t0 = "1",
##'                         st0 = "1"),
##'             constants = c(st0 = 0, sv.false = 1),
##'             match.map = list(M = list(s1 = 1, s2 = 2)),
##'             factors   = list(S = c("s1", "s2")),
##'             responses = c("r1", "r2"),
##'             type      = "lnorm")
##'
##' fac2df(model)
##'
##' model <- BuildModel(
##'   p.map     = list(A = "1", B = "R", t0 = "1", mean_v = c("F", "M"),
##'               sd_v = "M", st0 = "1"),
##'   match.map = list(M = list(s1 = 1, s2 = 2)),
##'   factors   = list(S = c("s1", "s2"), F = c("f1", "f2")),
##'   constants = c(sd_v.false = 1,st0 = 0),
##'   responses = c("r1", "r2"),
##'   type="norm")
##'
##' fac2df(model)
##' ##    S  F
##' ## 1 s1 f1
##' ## 2 s2 f1
##' ## 3 s1 f2
##' ## 4 s2 f2
##' @export
fac2df  <- function(x) {
  facs  <- lapply(strsplit(dimnames(x)[[1]], "\\."), function(xx){xx[-length(xx)]})
  facs  <- facs[ 1: ( length(facs) / length(attr(x, "responses")) ) ]
  fnams <- names(attr(x, "factors"))
  facs  <- data.frame(t(matrix(unlist(facs), nrow = length(fnams))))
  names(facs) <- fnams
  return(facs)
}

### Generic  ---------------------------------------------
##' Retrieve OS information
##'
##' A convenient wrapper to extract system information from \code{Sys.info}
##' and \code{.Platform}
##'
##' @examples
##' get_os()
##' ## sysname
##' ## "linux"
##' @export
get_os <- function() {
  sysinf <- Sys.info()
  ostype <- .Platform$OS.type

  ## Probe using Sys.info: Windows, Linux or Darwin
  if (!is.null(sysinf)) {
    os <- sysinf["sysname"]
    if (os == "Darwin") os <- "osx"
  } else {
    ## If something gets wrong with Sys.info, probe using .Platform
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))   os <- "osx"
    if (grepl("unix", R.version$os))      os <- "osx"
    if (grepl("linux-gnu", R.version$os)) os <- "linux"
    if (grepl("mingw32", R.version$os))   os <- "windows"
  }
  tolower(os)
}

getboundsR <- function(data) {
  model <- attr(data, "model")
  bound <- rep("lower", dim(model)[1])
  bound[as.vector(sapply(
    paste("", names(attr(model, "match.map")$M),
      attr(model, "match.map")$M, sep="*"),
    function(x){ grep(glob2rx(x), row.names(model)) }))] <- "upper"
  names(bound) <- row.names(model)
  bound
}

##' Calculate autocorrelation
##'
##' This function calculate autocorrealtion for a Markov Chain
##'
##' @param x likelihood vector
##' @param nlag number of autocorrelation lags
##' @importFrom stats cor
##' @export
ac <- function (x, nlag) {
  out <- data.frame(Lag = 1:nlag, Autocorrelation = cor(ac_(x, nlag),
    use = "pairwise.complete.obs")[, 1])
  return(out)
}


##' Extract parameter names from a model object
##'
##' @param model a model object
##'
##' @export
GetPNames <- function(model) { return(names(attr(model, "p.vector"))) }

checklba <- function(x) {
  m <- attr(x$data, "model")

  parnames <- attr(m, "par.names")
  if ( (which(parnames == "A") != 1) |
       (which(parnames == "B") != 2) |
       (which(parnames == "t0") != 3) |
       (which(parnames == "mean_v") != 4) |
       (which(parnames == "sd_v") != 5) |
       (which(parnames == "st0") != 6) ) {
    cat(parnames, "\n")
    message("parnames / p.vector must be ordered as: A, B, t0, mean_v, sd_v, & st0.")
    stop("Check p.map")
  }
}


