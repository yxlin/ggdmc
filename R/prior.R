### Prior density functions (plot_prior) ------
##' A modified dbeta function
##'
##' @param x quantile
##' @param p1 shape1 parameter
##' @param p2 shape2 parameter
##' @param lower lower bound
##' @param upper upper bound
##' @param lg logical; if TRUE, return log density.
##' @export
##' @importFrom stats dbeta
dbeta_lu <- function(x, p1, p2, lower, upper, lg = FALSE)
## Used with beta prior
{
  if (!lg) {
    out <- dbeta((x-lower)/(upper-lower), p1, p2, log = FALSE) / (upper-lower)
  } else {
    out <- dbeta((x-lower)/(upper-lower), p1, p2, log = TRUE) - log(upper-lower)
  }
  return(out)
}

##' A modified dgamma function
##'
##' @param x quantile
##' @param p1 shape parameter
##' @param p2 scale parameter
##' @param lower lower bound
##' @param upper upper bound
##' @param lg log density?
##' @importFrom stats dgamma
##' @export
dgamma_l <- function(x, p1, p2, lower, upper, lg = FALSE) {
  dgamma(x - lower, shape = p1, scale = p2, log = lg)
}

##' A modified dlnorm functions
##'
##' @param x quantile
##' @param p1 meanlog parameter
##' @param p2 sdlog parameter
##' @param lower lower bound
##' @param upper upper bound
##' @param lg log density?
##' @importFrom stats dlnorm
##' @export
dlnorm_l <- function(x, p1, p2, lower, upper, lg = FALSE) {
  dlnorm(x - lower, p1, p2, log = lg)
}


##' A modified dcauchy functions
##'
##' @param x quantile
##' @param p1 location parameter
##' @param p2 scale parameter
##' @param lg log density?
##' @importFrom stats dcauchy
##' @export
dcauchy_l <- function(x, p1, p2, lg = FALSE) {
  dcauchy(x, p1, p2, log = lg)
}


##' A pseudo constant function to get constant densities
##'
##' Used with constant prior
##'
##' @param x quantile
##' @param p1 constant value
##' @param p2 unused argument
##' @param lower dummy varlable
##' @param upper dummy varlable
##' @param lg log density?
##' @export
dconstant <- function(x, p1, p2, lower, upper, lg = FALSE) {
  den <- as.numeric((x == p1))
  if (lg) out <- base::log(den) else out <- den
  return(out)
}

##' @importFrom stats dunif
dunif_ <- function(x, p1, p2, lower, upper, lg = FALSE) {
  dunif(x, min = p1, max = p2, log = lg)
}


### Main functions --------------------
##' Specifying Prior Distributions
##'
##' \code{BuildPrior} sets up prior distributions for each model
##' parameter. \code{p1} and \code{p2} refer to the first and second parameters
##' a prior distribution. \code{p1} must comes with parameter names.
##'
##' Four distribution types are implemented:
##' \enumerate{
##' \item Normal and truncated normal distribution, where: p1 = mean, p2 = sd.
##'       When the lower and upper are not provided, they are set to -Inf and
##'       Inf, rendering a normal distribution. Type name is "tnorm".
##' \item Beta distribution, where: p1 = shape1 and p2 = shape2 (see \link{pbeta}).
##'       Note the uniform distribution is a special case of the beta with p1 = 1and
##'       and p2 = 1. Type name is "beta".
##' \item Gamma distribution, where p1 = shape and p2 = scale (see \link{pgamma}).
##'       Note p2 is scale, not rate. Type name is "gamma".
##' \item Log-normal, where p1 = meanlog and p2 = sdlog (see \link{plnorm}).
##' \item Uniform distribution. The bounds are not c(0, 1). The option comes handy.
##'       Type name is "unif".
##' }
##'
##' @param p1 the first parameter of a distribution
##' @param p2 the second parameter of a distribution
##' @param lower lower support (boundary)
##' @param upper upper support (boundary)
##' @param dists a vector of character string specifying a distribution.
##' @param untrans whether to do log transformation. Default is not
##' @param types available distribution types
##' @return a list of list
##' @examples
##' ## Show using dbeta to visualise a uniform distribution with bound (0, 1)
##' x <- seq(-.1, 1.1, .001)
##' plot(x, dbeta(x, 1, 1), type="l", ylab="Density", xlab="x", lwd=2)
##'
##' ## BuildPrior
##' pop.mean  <- c(a=2,   v=4, z=0.5, t0=0.3)
##' pop.scale <- c(a=0.5, v=.5, z=0.1, t0=0.05)
##'
##' pop.prior <- BuildPrior(
##'   dists = rep("tnorm", 4),
##'   p1    = pop.mean,
##'   p2    = pop.scale,
##'   lower = c(0,-5,  0, 0),
##'   upper = c(5, 7,  1, 1))
##'
##' p.prior <- BuildPrior(
##'   dists = rep("tnorm", 4),
##'   p1    = pop.mean,
##'   p2    = pop.scale*5,
##'   lower = c(0,-5, 0, 0),
##'   upper = c(5, 7, 1, 1))
##'
##' mu.prior <- BuildPrior(
##'   dists = rep("tnorm", 4),
##'   p1    = pop.mean,
##'   p2    = pop.scale*5,
##'   lower = c(0,-5,  0, 0),
##'   upper = c(5, 7,  1, 1))
##'
##' sigma.prior <- BuildPrior(
##'   dists = rep("beta", 4),
##'   p1    = c(a=1, v=1, z=1, t0=1),
##'   p2    = rep(1, 4),
##'   upper = rep(1, 4))
##'
##' ## Bind three priors together for hierarchical modelling
##' priors <- list(pprior=p.prior, location=mu.prior, scale=sigma.prior)
##'
##' @export
##' @importFrom methods new
BuildPrior <- function(p1, p2,
                       lower   = rep(NA, length(p1)),
                       upper   = rep(NA, length(p1)),
                       dists   = rep("tnorm", length(p1)),
                       untrans = rep("identity", length(p1)),
                       types   = c("tnorm", "beta", "gamma", "lnorm", "unif",
                                   "constant", "tnorm2", NA))
{
  npar <- length(p1)
  if (length(p2) == 1) { p2 <- rep(p2, npar) }
  name.untrans <- check_BuildPrior(p1, p2, lower, upper, dists, untrans, types)
  out <- vector("list", npar)
  pnames <- names(p1)
  names(out) <- pnames
  if ( is.null(pnames) ) stop("p1 must comes with parameter names")

  for (i in 1:npar) {
    out[[i]] <- switch(dists[i],
                       tnorm = {
                         if (is.na(lower[i])) lower[i] <- -Inf
                         if (is.na(upper[i])) upper[i] <- Inf
                         p <- c(p1[i], p2[i], lower[i], upper[i])
                         names(p) <- c("p1", "p2", "lower", "upper")
                         p <- as.list(p)
                         attr(p, "dist") <- 1  ## "tnorm"
                         p
                       },
                       beta  = {
                         if (is.na(lower[i])) lower[i] <- 0
                         if (is.na(upper[i])) upper[i] <- 1
                         p <- c(p1[i], p2[i], lower[i], upper[i])
                         names(p) <- c("p1", "p2", "lower","upper")
                         p <- as.list(p)
                         attr(p, "dist") <- 2 ##"beta_lu"
                         p
                       },
                       gamma = {
                         if (is.na(lower[i])) lower[i] <- 0
                         if (is.na(upper[i])) upper[i] <- Inf
                         p <- c(p1[i], p2[i], lower[i], upper[i])
                         names(p) <- c("p1", "p2", "lower", "upper")
                         p <- as.list(p)
                         attr(p, "dist") <- 3 ## "gamma_l"
                         p
                       },
                       lnorm = {
                         if (is.na(lower[i])) lower[i] <- 0
                         if (is.na(upper[i])) upper[i] <- Inf
                         p <- c(p1[i], p2[i], lower[i], upper[i])
                         names(p) <- c("p1", "p2", "lower", "upper")
                         p <- as.list(p)
                         attr(p, "dist") <- 4 ## "lnorm_l"
                         p
                       },
                       unif  = {
                         p <- c(p1[i], p2[i], -Inf, Inf)
                         names(p) <- c("p1", "p2", "lower", "upper")
                         p <- as.list(p)
                         attr(p, "dist") <- 5 ## "unif_"
                         p
                       },
                       constant = {
                         p <- c(p1[i], p2[i], -Inf, Inf)
                         names(p) <- c("p1", "p2", "lower", "upper")
                         p <- as.list(p)
                         attr(p, "dist") <- 6 ## "constant"
                         p
                       },
                       tnorm2 = {
                         if (is.na(lower[i])) lower[i] <- -Inf
                         if (is.na(upper[i])) upper[i] <- Inf
                         p <- c(p1[i], p2[i], lower[i], upper[i])
                         names(p) <- c("p1", "p2", "lower", "upper")
                         p <- as.list(p)
                         attr(p, "dist") <- 7  ## "tnorm2"
                         p
                       },

                       {
                         p <- c(p1[i], p2[i], -Inf, Inf)
                         names(p) <- c("p1", "p2", "lower", "upper")
                         p <- as.list(p)
                         attr(p, "dist") <- NA
                         p
                       }
    )
    out[[i]]$lg <- TRUE

    if (!name.untrans) {
      attr(out[[i]], "untrans") <- untrans[i]
    } else {
      if (is.na(untrans[pnames[i]])) {
        attr(out[[i]], "untrans") <- "identity"
      } else {
        attr(out[[i]], "untrans") <- untrans[pnames[i]]
      }
    }
  }


  out <- new("prior",
             npar = npar,
             pnames = pnames,
             priors = out)
  return(out)
}

BuildPrior_string <- function(p1, p2,
                           lower   = rep(NA, length(p1)),
                           upper   = rep(NA, length(p1)),
                           dists   = rep("tnorm", length(p1)),
                           untrans = rep("identity", length(p1)),
                           types   = c("tnorm", "beta", "gamma", "lnorm", "unif", "constant", "empty"))
{
  np1 <- length(p1)
  if (length(p2) == 1) p2 <- rep(p2, np1)
  name.untrans <- check_BuildPrior(p1, p2, lower, upper, dists, untrans, types)
  out <- vector(mode = "list", length = np1)
  pnames <- names(p1)
  names(out) <- pnames

  for (i in 1:np1) {
    out[[i]] <- switch(dists[i],
                       tnorm = {
                         if (is.na(lower[i])) lower[i] <- -Inf
                         if (is.na(upper[i])) upper[i] <- Inf
                         p <- c(p1[i], p2[i], lower[i], upper[i])
                         names(p) <- c("mean", "sd", "lower", "upper")
                         p <- as.list(p)
                         attr(p, "dist") <- "tnorm"
                         p
                       },
                       beta  = {
                         if (is.na(lower[i])) lower[i] <- 0
                         if (is.na(upper[i])) upper[i] <- 1
                         p <- c(p1[i], p2[i], lower[i], upper[i])
                         names(p) <- c("shape1", "shape2", "lower","upper")
                         p <- as.list(p)
                         attr(p, "dist") <- "beta_lu"
                         p
                       },
                       gamma = {
                         if (is.na(lower[i])) lower[i] <- 0
                         if (is.na(upper[i])) upper[i] <- Inf
                         p <- c(p1[i], p2[i], lower[i], upper[i])
                         names(p) <- c("shape", "scale", "lower", "upper")
                         p <- as.list(p)
                         attr(p, "dist") <- "gamma_l"
                         p
                       },
                       lnorm = {
                         if (is.na(lower[i])) lower[i] <- 0
                         if (is.na(upper[i])) upper[i] <- Inf
                         p <- c(p1[i], p2[i], lower[i], upper[i])
                         names(p) <- c("meanlog","sdlog", "lower", "upper")
                         p <- as.list(p)
                         attr(p, "dist") <- "lnorm_l"
                         p
                       },
                       unif  = {
                         if (is.na(lower[i])) lower[i] <- p1[i]
                         if (is.na(upper[i])) upper[i] <- p2[i]
                         p <- c(p1[i], p2[i], lower[i], upper[i])
                         names(p) <- c("min", "max", "lower", "upper")
                         p <- as.list(p)
                         attr(p, "dist") <- "unif_"
                         p
                       },
                       constant = {
                         p <- c(p1[i], p2[i], -Inf, Inf)
                         names(p) <- c("p1", "p2", "lower", "upper")
                         p <- as.list(p)
                         attr(p, "dist") <- "constant"
                         p
                       },
                       tnorm2 = {
                         p <- c(p1[i], p2[i], -Inf, Inf)
                         names(p) <- c("p1", "p2", "lower", "upper")
                         p <- as.list(p)
                         attr(p, "dist") <- "tnorm2"
                         p
                       },
                       {
                         p <- c(p1[i], p2[i], -Inf, Inf)
                         names(p) <- c("p1", "p2", "lower", "upper")
                         p <- as.list(p)
                         attr(p, "dist") <- NA
                         p
                       }
    )
    out[[i]]$lg <- TRUE

    if (!name.untrans) {
      attr(out[[i]], "untrans") <- untrans[i]
    } else {
      if (is.na(untrans[pnames[i]])) {
        attr(out[[i]], "untrans") <- "identity"
      } else {
        attr(out[[i]], "untrans") <- untrans[pnames[i]]
      }
    }
  }

  class(out) <- c("prior", "list")
  return(out)
}



check_BuildPrior <- function(p1, p2, lower, upper, dists, untrans,types) {
  np1 <- length(p1)  ## number of parameter
  np2 <- length(p2)

  if (np2 == 1) { p2 <- rep(p2, np1); np2 <- np1 }
  if (np1 != np2) stop("p1 and p2 must be equal length")
  if (np1 != length(lower) ) stop("p1 and lower must be equal length")
  if (np1 != length(upper) ) stop("p1 and upper must be equal length")
  if (np1 != length(dists) ) stop("p1 and dists must be equal length")

  both.not.na <- !is.na(upper) & !is.na(lower)

  if ( any(upper[both.not.na] <= lower[both.not.na]) )
    stop("All elements of upper must be greater than lower")
  if ( !all(dists %in% types) )
    stop(paste("Unsupported distribution, allowable types are:",
               paste(types, collapse = ", ")))

  name.untrans <- length(untrans) != np1
  if (name.untrans & (is.null(names(untrans)) | is.null(names(p1))))
    stop("If untrans vector is not the same length as p1 it must have p1 names")
  if (!(all(names(untrans) %in% names(p1)))) {
    stop("untrans vector has names not in p1 names")
  }
  return(name.untrans)
}

