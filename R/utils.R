### Hierarchical tools ------------------------------------------------
is_multimodel <- function(model, ns = NA)
## Used in simulate_many
{
  if (is.list(model))
  {
    if (length(model) != ns)
      stop("number of participants not equal to number of models")
    out <- TRUE
  }
  else
  {
    if (is.na(ns)) stop("Must indicate the number of participants")
    out <- FALSE
  }
  return(out)
}

##' Extract trial log likelihoods
##'
##' This function simply run trial_loglik to loop through one subject after
##' another to extracts trial_log_likes from a list of subject fits and
##' concatanates the result into an array.
##'
##' @param samples posterior samples
##' @param thin thinnng length
##' @param verbose whether print information
##' @export
trial_loglik_hier <- function(samples, thin = 1, verbose=FALSE)
{
  check <- function(x) {
    nmc <- min(x[1,])

    if ( !all(x[3,1]==x[3,-1]) )
      warning(paste("Subjects do not all have the same number of interations, using first",
                    nmc,"for all."))

    if ( !all(x[2,1]==x[2,-1]) )
      stop("Subjects must have the same number of chains")
  }

  PrintSize <- function(x, thin) {

    tll <- trial_loglik(x@individuals[[1]], thin)

    nsub <- length(x@snames)
    sdim <- dim(tll); ## sdim ##  ntrial x nchain  x nnmc_thin
    size <- sum(nsub * prod(sdim))

    ## nchain and nmc per subject
    nmc_nchain_mat  <- sapply(x, function(xx) { dim(xx@log_likelihoods) } )

    dimnames(nmc_nchain_mat) <- list(c("Chains", "Trials"), names(x))
    message("Log-likelihood dimension")
    print(nmc_nchain_mat)

    names(sdim) <- c("Trials","Chains","Iterations")
    message("\nSubject 1")
    print(sdim)
    message("Total log-likelihoods: ", appendLF=FALSE)
    cat(round(size/1e6,2), "millions)\n")
    return(NULL)
  }

  if(verbose) PrintSize(samples, thin)

  tlls <- lapply(samples@individuals, trial_loglik, thin)

  sdims <- sapply(tlls, dim)  ## DMC nmc_thin, nchain, ntrial
  nmc <- min(sdims[3,])       ## nnmc_thin
  check(sdims)

  ### nmc_thin x nchain x (ntrial x nsub)
  out <- array(dim=c(nmc, dim(tlls[[1]])[2], sum(sdims[1,])))
  nsub <- length(samples@snames)

  start <- 1;
  end <- sdims[1,1]

  for (i in 1:nsub)
  {
    out[,,start:end] <- aperm(tlls[[i]], c(3, 2, 1))
    if (i < nsub)
    {
      start <- end+1
      end   <- start - 1 + sdims[1, i+1]
    }
  }
  return(out)
}


##' Calculate the autocorrelation of a vector
##'
##' Calculate the autocorrelation of a vector.
##'
##' @param x a vector storing parameter values
##' @param nLags the maximum number of lags
##' @return A data.frame
##' @export
##' @examples
##' res <- ac(1:100)
##' ## List of 2
##' ## $ Lag            : int [1:50] 1 2 3 4 5 6 7 8 9 10 ...
##' ## $ Autocorrelation: num [1:50] 1 1 1 1 1 1 1 1 1 1 ...
##'
##' res <- ac(rnorm(100))
##' str(res)
##' ## List of 2
##' ## $ Lag            : int [1:50] 1 2 3 4 5 6 7 8 9 10 ...
##' ## $ Autocorrelation: num [1:50] 1 -0.0485 0.0265 -0.1496 0.0437 ...
ac <- function(x, nLags = 50)
{
  tmp <- ac_(x, nLags)
  return(  list(Lag = 1:nLags, Autocorrelation=tmp[,1]))
}


### Generic  ---------------------------------------------
##' Retrieve information of operating system
##'
##' A wrapper function to extract system information from \code{Sys.info}
##' and \code{.Platform}
##'
##' @examples
##' get_os()
##' ## sysname
##' ## "linux"
##' @export
get_os <- function()
{
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



##' @importFrom methods slot
checklba <- function(x)
{
  model <- slot(x, "model")
  if (slot(model, "type") == "norm" )
  {
    parnames <- slot(model, "par.names")
    if ( (which(parnames == "A")      != 1) |
         (which(parnames == "B")      != 2) |
         (which(parnames == "t0")     != 3) |
         (which(parnames == "mean_v") != 4) |
         (which(parnames == "sd_v")   != 5) |
         (which(parnames == "st0")    != 6) )
    {
      cat("Your p.vector is order as: ", parnames, "\n")
      message("It must be in the order of: A, B, t0, mean_v, sd_v, & st0.")
      stop("Check p.map")
    }
  }

}

##' Extract parameter names from a model object
##'
##' GetPNames will be deprecated. Please extract pnames directly via S4 slot 'model@pnames'
##' @param x a model object
##'
##' @export
GetPNames <- function(x) {
  warning("GetPNames will be deprecated. Please extract pnames directly via S4 slot 'model@pnames'")
  return(x@pnames)
}

