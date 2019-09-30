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
##' Specifying Parameter Prior Distributions
##'
##' \code{BuildPrior} sets up parameter prior distributions for each model
##' parameter. \code{p1} and \code{p2} refer to the first and second parameters
##' a prior distribution.
##'
##' Four distribution types are implemented:
##' \enumerate{
##' \item Normal and truncated normal, where: p1 = mean, p2 = sd. It specifies
##' a normal distribution when bounds are set -Inf and Inf,
##' \item Beta, where: p1 = shape1 and p2 = shape2 (see \link{pbeta}). Note the
##'       uniform distribution is a special case of the beta with p1 and
##'       p2 = 1),
##' \item Gamma, where p1 = shape and p2 = scale (see \link{pgamma}). Note p2 is
##'       scale, not rate,
##' \item Lognormal, where p1 = meanlog and p2 = sdlog (see \link{plnorm}).
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
##' @export
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

  class(out) <- c("list", "prior")
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

##' Parameter Prior Distributions
##'
##' Probability density functions and random generation for parameter prior
##' distributions.
##'
##' @param prior a list of list usually created by BuildPrior to store the
##' information about parameter prior distributions.
##' @param n number of observations/random draws
##'
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
##' pvec <- c(a=1, v=1, z=0.5, sz=0.25, sv=0.2,t0=.15)
##' p.prior  <- BuildPrior(
##'   dists = rep("tnorm", 6),
##'   p1    = c(a=2,   v=2.5, z=0.5, sz=0.3, sv=1,  t0=0.3),
##'   p2    = c(a=0.5, v=.5,  z=0.1, sz=0.1, sv=.3, t0=0.05) * 5,
##'   lower = c(0,-5, 0, 0, 0, 0),
##'   upper = c(5, 7, 2, 2, 2, 2))
##'
##' @export
rprior <- function(prior, n = 1)
{
  rprior_mat(prior, n)
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

### Prior tools ------------------------------------------------------------
##' Print Prior Distribution
##'
##' a convenient function to rearrange \code{p.prior} or an element in a
##' \code{pp.prior} as a data frame for inspection.
##'
##' @param x a list of prior distributions list, usually created by
##' \code{BuildPrior}
##' @param ... other arguments
##' @return a data frame listing prior distributions and their settings
##' @examples
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
##' @export
print.prior <- function(x, ...) {

  ncol <- 7
  npar <- length(x);
  bucket  <- matrix(numeric(npar*ncol), npar);

  for(i in 1:npar) {
    add1    <- attr(x[[i]], "dist");
    add2    <- attr(x[[i]], "untrans");
    if (is.na(add1)) {
      tmp <- unlist(x[[i]])
      rowObj  <- c(NA, NA, NA, NA, tmp[5], NA, add2)
    } else if (add1 == "constant") {
      tmp <- unlist(x[[i]])
      rowObj  <- c(c(tmp[1], tmp[2], NA, NA, tmp[5]), add1, add2);
    } else if (add1 == "gamma_l") {
      tmp <- unlist(x[[i]])
      rowObj <- c(tmp[1:5], "gamma_l", add2)
    } else {
      rowObj  <- c(unlist(x[[i]]), add1, add2);
    }
    bucket[i,] <- rowObj
  }

  out <- data.frame(bucket)
  names(out) <- c("p1", "p2", "lower", "upper", "lg", "dist", "untrans")
  rownames(out) <- names(x)
  return(out)
}

### Back up (to be removed) ------------------------------------------------------
