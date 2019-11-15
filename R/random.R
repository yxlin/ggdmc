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

  # Build parameter matrix (and divide a, v, and sv, by s)
  params <- cbind (a, v, z, d, sz, sv, t0, st0, s)

  # Check for illegal parameter values
  if(ncol(params)<8)
    stop("Not enough parameters supplied: probable attempt to pass NULL values?")
  if(!is.numeric(params))
    stop("Parameters need to be numeric.")
  if (any(is.na(params)) || !all(is.finite(params)))
    stop("Parameters need to be numeric and finite.")

  randRTs    <- vector("numeric",length=n)
  randBounds <- vector("numeric",length=n)

  parameter_char <- apply(params, 1, paste0, collapse = "\t")
  parameter_factor <- factor(parameter_char, levels = unique(parameter_char))
  parameter_indices <- split(seq_len(n), f = parameter_factor)

  for (i in seq_len(length(parameter_indices)))
  {
    ok_rows <- parameter_indices[[i]]

    # Calculate n for this row
    current_n <- length(ok_rows)

    out <- r_fastdm (current_n,
                     params[ok_rows[1],1:9],
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
##' Random number generation
##'
##' A wrapper function for generating random numbers of either
##' the model type, \code{rd}, \code{norm}, \code{norm_pda},
##' \code{norm_pda_gpu}, or \code{cddm}. \code{pmat} is generated usually by
##' \code{TableParameter}.
##'
##' Note PM model uses \code{norm} type.
##'
##' @param type a character string of the model type
##' @param pmat a matrix of response x parameter
##' @param n number of observations. This can be a scalar or a integer vector.
##' @param seed an integer specifying a random seed
##' @param ... other arguments
##' @export
random <- function(type, pmat, n, seed = NULL, ...)
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
    ## posdrift is always TRUE in simulation
    ## pmat: A b t0 mean_v sd_v st0 (nacc)
    if (ncol(pmat) == 7) ## PM model
    {
      ## The only difference is to selectively pick 2 or 3 accumulators.
      ## This has been done by BuildModel and p_df.
      nacc <- pmat[1,7]
      A <- pmat[1:nacc,1]
      b <- pmat[1:nacc,2]
      t0 <- pmat[1:nacc,3]
      mean_v <- pmat[1:nacc,4]
      sd_v <- pmat[1:nacc,5]
      st0 <- pmat[1:nacc,6]
      posdrift <- TRUE

      out <- rlba_norm(n, A, b, mean_v, sd_v, t0, st0, posdrift)
    }
    else
    {
      out <- rlba_norm(n, pmat[, 1], pmat[, 2], pmat[, 4], pmat[, 5],
                       pmat[,3], pmat[1,6], TRUE)
    }

  } else if (type=="cddm") {
    # nw <- nrow(pmat)
    pvec <- as.vector(t(pmat[1,]))
    out  <- rcircle(n=n, P=pvec[1:8], tmax=pvec[9], h=pvec[10],
                    nw = nrow(pmat))
  } else {
    stop("Model type yet created")
  }

  attr(out, "seed") <- seed
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
"%w/o%" <- function(x, y) x[!x %in% y] #--  x without y


LabelTheta <- function(dat, response)
{
  ## This is a new function for cddm only
  nw <- length(response)
  w <- 2*pi/nw
  theta <- seq(-pi, pi-w, w)
  R <- rep(NA, nrow(dat))
  for(i in 1:nw)
  {
    idx <- which(round(dat$R, 2) == round(theta[i], 2))
    R[idx] <- response[i]
  }

  out <- factor(R, levels=response)
  return(out)
}