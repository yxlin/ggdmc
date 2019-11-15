






##' Unstick posterios samples (One subject)
##'
##' @param x posterior samples
##' @param bad a numeric vector, indicating which chains to remove
##' @export
unstick_one <- function(x, bad) {
  # cat("unstick_one")

  nchain <- x$n.chains
  if (length(bad) > 0)
  {
    if (!all(bad %in% 1:nchain))
      stop(paste("Index of bad chains must be in 1 to ", nchain))

    x$theta            <- x$theta[,-bad,]
    x$summed_log_prior <- x$summed_log_prior[-bad,]
    x$log_likelihoods  <- x$log_likelihoods[-bad,]
    x$n.chains         <- x$n.chains - length(bad)
  }

  return(x)
}


##' Model checking functions
##'
##' The function tests whether we have drawn enough samples.
##'
##' @param x posterior samples
##' @param minN specify the size of minimal effective samples
##' @param nfun specify to use the \code{mean} or \code{median} function to
##' calculate effective samples
##' @param verbose print more information
##' @export
iseffective <- function(x, minN, nfun, verbose = FALSE) {
  n <- do.call(nfun, list(effectiveSize(x, verbose = verbose)))
  fail <- n < minN
  if (verbose) {
    cat("Length check")
    if (!fail) cat(": OK\n") else cat(paste(":",n,"\n"))
  }
  fail
}

##' @rdname PickStuck-methods
CheckConverged <- function(x)
{
  stuck <- isstuck(x, verbose = FALSE, cut = 10)
  flat  <- isflat(x, p1 = 1/3, p2 = 1/3,
                  cut_location = 0.25, cut_scale = Inf, verbose = FALSE)
  mix  <- ismixed(x, cut = 1.05, verbose = FALSE)
  size <- iseffective(x, minN = 500, nfun = "mean", FALSE)
  isstuck <- TRUE
  if (stuck == 0) isstuck <- FALSE

  out <- c(isstuck, flat, mix, size)
  names(out) <- c("Stuck", "Flat", "Mix", "ES")
  return(out)
}


