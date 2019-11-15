### MCMC -------------------------------------------------------
##' @importFrom methods slot
extractPhi <- function(x, start, end, subchain) {

  thin   <- x@thin
  nchain <- x@nchain
  pnames <- x@pnames
  iter <- start:end
  niter <- length(iter)
  chains <- 1:nchain

  location <- x@prior_loc
  scale    <- x@prior_sca
  not_na0 <- sapply(slot(location, "priors"), function(y) !is.na(attr(y, "dist")))
  not_na1 <- sapply(slot(scale,    "priors"), function(y) !is.na(attr(y, "dist")))

  ok1 <- paste(pnames[not_na0], "h1", sep = ".")
  ok2 <- paste(pnames[not_na1], "h2", sep = ".")

  location_names <- pnames[not_na0]
  scale_names    <- pnames[not_na1]
  new_loc_names <- paste(location_names, "h1", sep=".")
  new_sca_names <- paste(scale_names,    "h2", sep=".")

  if (!anyNA(subchain)) {
    if ( any(x@nchain < subchain) ) stop("Chains not in the range.")
    chains <- subchain
    nchain <- length(chains)
    cat("Calculate chains:", chains, "\n")
  }

  v <- lapply(1:nchain, function(k){

    tmp1 <- t( x@phi_loc[not_na0, chains[k], iter] )
    tmp2 <- t( x@phi_sca[not_na1, chains[k], iter] )
    dimnames(tmp1)[[2]] <- new_loc_names
    dimnames(tmp2)[[2]] <- new_sca_names

    cbind(tmp1, tmp2)
  })

  if (length(ok1) != x@npar | length(ok2) != x@npar) {
    npar  <- sum(not_na0, not_na1)
    newpnames <- c( new_loc_names[not_na0], new_sca_names[not_na1] )
  } else {
    npar <- x@npar * 2
    newpnames <- c(new_loc_names, new_sca_names)
  }

  attr(v, "npar")   <- npar
  attr(v, "pnames") <- newpnames
  attr(v, "thin")   <- thin
  attr(v, "nchain") <- nchain
  attr(v, "iter")   <- iter
  attr(v, "chains") <- chains
  return(v)
}

##' @importFrom stats var
##' @importFrom matrixStats rowMeans2
##' @importFrom methods slot
getW <- function (x) {
  theta  <- slot(x, "theta")
  nchain <- slot(x, "nchain")
  npar   <- slot(x, "npar")
  v <- lapply (1:nchain, function(k) { t(theta[,k,]) } )
  tmp1 <- sapply(v, var)
  tmp2 <- matrixStats::rowMeans2(tmp1)
  W  <- matrix(tmp2, nrow=npar)
  S2 <- array(tmp1, dim = c(npar, npar, nchain))
  return(list(W=W, S2=S2))
}

##' @importFrom matrixStats rowMeans2
##' @importFrom methods slot
get_xbar <- function (x) {
  theta  <- slot(x, "theta")
  nchain <- slot(x, "nchain")
  npar   <- slot(x, "npar")
  v   <- lapply(1:nchain, function(k){ matrixStats::rowMeans2(theta[,k,]) })
  out <- matrix( unlist(v), nrow=npar)
  return(out)
}

get_mpsrf <- function(npar, niter, W, B, multivariate) {

  if (npar > 1 && multivariate) {

    if (is.R()) {
      CW   <- chol(W)
      emax <- eigen(backsolve(CW, t(backsolve(CW, B, transpose = TRUE)),
                    transpose = TRUE), symmetric = TRUE,
                    only.values = TRUE)$values[1]
    }
    else {
      emax <- eigen(qr.solve(W, B), symmetric = FALSE, only.values = TRUE)$values
    }

    mpsrf <- sqrt((1 - 1/niter) + (1 + 1/npar) * emax/niter)
  } else {
    mpsrf <- NULL
  }

  return(mpsrf)

}


### Summary ------------------------------------------------------
##' @importFrom stats ar
##' @importFrom stats residuals
##' @importFrom stats sd
##' @importFrom stats lm
spectrum0_ar <- function (x)
{
  x <- as.matrix(x)
  v0 <- order <- numeric(ncol(x))
  names(v0) <- names(order) <- colnames(x)
  z <- 1:nrow(x)

  for (i in 1:ncol(x)) {
    lm.out <- stats::lm(x[, i] ~ z)
    if (identical(all.equal(sd(residuals(lm.out)), 0), TRUE)) {
      v0[i] <- 0
      order[i] <- 0
    }
    else {
      ar.out <- ar(x[, i], aic = TRUE)
      v0[i] <- ar.out$var.pred/(1 - sum(ar.out$ar))^2
      order[i] <- ar.out$order
    }
  }
  return(list(spec = v0, order = order))
}

safespec0 <- function(x) {
  result <- try(spectrum0_ar(x)$spec)
  if (class(result) == "try-error") result <- NA
  if (class(result) == "try") result <- NA
  result
}

##' @importFrom methods new
summary_hyper <- function(object, start, end, type, prob, digits, verbose)
{
  # x <- fit
  # end <- NA
  # start <- 1
  # subchain <- NA
  # type <- 1
  # prob <- c(.025, .5, .975)

  if (is.na(end)) end <- object@nmc
  res <- check_nonna(object, type)

  samples <- new("posterior",
                 theta  = res[[1]],
                 ## Not used dummy ##
                 summed_log_prior = object@individuals[[1]]@summed_log_prior,
                 log_likelihoods  = object@individuals[[1]]@log_likelihoods,
                 dmi              = object@individuals[[1]]@dmi,
                 prior            = object@individuals[[1]]@prior,
                 ## Not used dummy ##
                 start  = object@start,
                 npar   = res[[2]],
                 pnames = res[[3]],
                 nmc    = object@nmc,
                 thin   = object@thin,
                 nchain = object@nchain)
  out <- summary_one(samples, start, end, prob)
  return(out)
}

##' @importFrom methods slot
summary_one <- function(x, start, end, prob) {
  # x <- fit0

  if ( is.na(end) ) end <- slot(x, "nmc")
  if ( length(prob) == 1) stop("prob must be more than one element")

  nchain <- slot(x, "nchain")
  iter   <- start:end
  niter  <- length(iter)

  statnames <- c("Mean", "SD", "Naive SE", "Time-series SE")
  varstats  <- matrix(nrow = x@npar, ncol = length(statnames),
                     dimnames = list(x@pnames, statnames))
  xtsvar    <- matrix(nrow = nchain, ncol = x@npar)
  x0        <- NULL

  for (k in 1:nchain)
  {
    for (i in 1:x@npar)
    {
      xtsvar[k, i] <- safespec0(x@theta[i,k,]) ## npar x nchain x nmc
    }
    x0 <- cbind(x0, x@theta[,k,])
  }

  xmean    <- matrixStats::rowMeans2(x0)
  xvar     <- matrixStats::rowVars(x0)
  xtsvar2  <- matrixStats::colMeans2(xtsvar)
  varquant <- matrixStats::rowQuantiles(x0, probs = prob)

  varstats[, 1] <- xmean
  varstats[, 2] <- sqrt(xvar)
  varstats[, 3] <- sqrt(xvar/niter * nchain)
  varstats[, 4] <- sqrt(xtsvar2/niter * nchain)
  varquant <- drop(varquant)
  varstats <- drop(varstats)
  rownames(varquant) <- x@pnames

  out <- list(statistics = varstats, quantiles = varquant)
  attr(out, "start")  <- start
  attr(out, "end")    <- end
  attr(out, "thin")   <- x@thin
  attr(out, "nchain") <- nchain
  return(out)
}

##' @importFrom matrixStats colMeans2
summary_many <- function(x, start, end, prob, digits, verbose) {

  tmp0 <- lapply(x, summary_one, start, end, prob)

  if (verbose) {
    message("Show and keep the estimates of the mean of 50% quantile ")
    message("Use verbose = FALSE & recovery = FALSE to retrieve all information")

    v <- lapply(tmp0, function(xx) { xx[[1]][, 1] } )
    vmat <- do.call("rbind", v)
    out <- rbind(do.call("rbind", v), matrixStats::colMeans2(vmat))
    row.names(out) <- c(names(x), "Mean")
    print( round(out, digits))
  } else {
    out <- tmp0
  }

  return(out)
}

summary_recoverone <- function(object, start, end, ps, digits, prob, verbose)
{
#
#   object <- samples
#   start <- 1
#   end <- 500
#   prob <- c(0.025, 0.25, 0.5, 0.75, 0.975)
#   ps = pop.mean
#   digits <- 2
#   verbose <- T

  qs <- summary_one(object, start, end, prob)$quantiles
  parnames <- dimnames(qs)[[1]]

  if (!is.null(ps) && ( !all(parnames %in% names(ps)) ) )
    stop("Names of p.vector do not match parameter names in samples")

  est  <- qs[names(ps), "50%"]
  op.vector <- ps[order(names(ps))]
  oest <- est[order(names(est))]
  bias <- oest- op.vector

  lo  <- qs[names(ps), "2.5%"]
  hi  <- qs[names(ps), "97.5%"]
  olo <- lo[order(names(lo))]
  ohi <- hi[order(names(hi))]

  out  <- rbind(
    'True'          = op.vector,
    '2.5% Estimate' = olo,
    '50% Estimate'  = oest,
    '97.5% Estimate'= ohi,
    'Median-True'   = bias)

  if (verbose) print(round(out, digits))
  return(out)
}

summary_recovermany <- function(object, start, end, ps, digits, prob, verbose)
{
  # object <- fit
  # start <- 1
  # end <- NA
  # digits <- 2
  est <- summary_many(object, start, end, prob, digits, TRUE)

  mean.est <- matrixStats::colMeans2(est)
  mean.ps <- matrixStats::colMeans2(ps)
  sd.est <- matrixStats::colSds(est)
  sd.ps <- matrixStats::colSds(ps)

  loc <- rbind(mean.est, mean.ps, mean.ps - mean.est)
  sca <- rbind(sd.est, sd.ps, sd.ps - sd.est)
  out <- rbind(loc, sca)

  rownames(out) <- c("Mean", "True", "Diff", "Sd", "True", "Diff")
  colnames(out) <- object[[1]]@pnames

  if (verbose) print(round(out, digits))
  return(list(summary=out, estimate = est))
}

check_nonna <- function(x, type) {
  # x <- fit
  # type <- 1

  if (type==1) {
    not_na_idx <- sapply(x@prior_loc@priors, function(y) !is.na( attr(y, "dist")))
  } else {
    not_na_idx <- sapply(x@prior_loc@priors, function(y) !is.na( attr(y, "dist")))
  }

  ## npar x nchain x nmc
  if (type==1) theta <- x@phi_loc[not_na_idx,,] else theta <- x@phi_sca[not_na_idx,,]
  pnames  <- x@pnames[not_na_idx]
  newnpar <- sum(not_na_idx)

  return(list(theta, newnpar, pnames))
}

##' @importFrom methods new
summary_recoverhyper <- function(object, start, end, ps, type, digits, prob,
  verbose) {
  # object <- fit
  # type <- 1

  res   <- check_nonna(object, type)

  ## The selected phi (either location or scale renamed as theta)
  samples <- new("posterior",
                 theta  = res[[1]],
                 ## Not used ##
                 summed_log_prior = object@individuals[[1]]@summed_log_prior,
                 log_likelihoods  = object@individuals[[1]]@log_likelihoods,
                 dmi              = object@individuals[[1]]@dmi,
                 prior            = object@individuals[[1]]@prior,
                 ## Not used ##
                 start  = object@start,
                 npar   = res[[2]],
                 pnames = res[[3]],
                 nmc    = object@nmc,
                 thin   = object@thin,
                 nchain = object@nchain)
  return( summary_recoverone(samples, start, end, ps, digits, prob, verbose) )
}

### Model Selection-------------------------------------------------
##' Calculate the statistics of model complexity
##'
##' Calculate deviance for a model object for which a
##' log-likelihood value can be obtained, according to the formula
##' -2*log-likelihood.
##'
##' @param object posterior samples
##' @param start start iteration
##' @param end end iteration
##' @param ... other plotting arguments passing through dot dot dot.
##' @references
##' Spiegelhalter, D. J., Best, N. G., Carlin, B. P., & van der Linde, A.
##' (2002). Bayesian Measures of Model Complexity and Fit. Journal of the Royal
##' Statistical Society, Series B (Statistical Methodology), 64(4), 583--639.
##' doi:10.1111/1467-9868.00353\cr
##'
##' Ando, T. (2007). Bayesian predictive information criterion for the
##' evaluation of hierarchical Bayesian and empirical Bayes models.
##' Biometrika. 94(2), 443–458. doi:10.1093/biomet/asm017.
##' @importFrom stats var
##' @export
deviance_model <- function(object, start, end, ...)
{
  ## The deviance; constant term drop when comparing models
  ## See equation (10) "where" (p. 587) and
  ## https://en.wikipedia.org/wiki/Deviance_information_criterion
  D <- -2*object@log_likelihoods[,start:end]

  ## Average across chains and iterations
  ## the expectation of theta; ie the posterior mean of the parameters
  mtheta <- apply(object@theta[,, start:end], 1, mean)
  ## The deviance based on the expectation of theta
  Dmean <- -2*sum(log(likelihood(mtheta, object@dmi)))

  # A list with mean (meanD), variance (varD) and min (minD)
  # of Deviance, deviance of mean theta (Dmean)
  out <- list(np=object@npar, meanD = mean(D), varD = var(as.vector(D)),
              minD = min(D), Dmean = Dmean, D = D)
  return(out)
}

pd <- function(ds)
{
  # ds is an object returned by deviance.model
  # effective number of parameters calculated by mean, min and var methods
  list(Pmean=ds$meanD-ds$Dmean,Pmin=ds$meanD-ds$minD,Pvar=ds$varD/2)
}

##' @importFrom stats density
##' @importFrom graphics hist
##' @importFrom graphics abline
posterior_logLik <- function(D1, D2, main="", plot=FALSE, den=TRUE)
{
  # Aitkin, M., Boys, R. J., & Chadwick, T. (2005). Bayesian point null
  # hypothesis testing via the posterior likelihood ratio. Statistics and
  # Computing, 15(3), 217-230.
  if (is.list(D1)) D1 <- D1$D
  if (is.list(D2)) D2 <- D2$D
  n <- base::pmin(length(D1), length(D2))
  dD <- D1[1:n] - D2[1:n]

  if (plot)
  {
    if (den) {
      plot(stats::density(dD),xlab="D1-D2",main=main)
    } else {
      graphics::hist(dD,breaks="fd",xlab="D1-D2",main=main)
    }

    if (min(dD) < 0 & max(dD) > 0) graphics::abline(v=0) ## Difference covers 0
  }

  ## Return the probability of difference < 0
  return( c(pD1 = mean(dD < 0)) )
}

weightIC <- function(ics, BPIC=FALSE)
{
  # Calculate weights for a set of models
  # ics is a vector of ICs
  d <- ics - min(ics)
  w <- exp(-d/2)/sum(exp(-d/2))

  if (!is.null(names(ics))) {
    mnams <- names(ics)
  } else {
    mnams <- paste0("M", 1:length(d))
  }

  out <- matrix(data = c(d, w),
                nrow = length(d),
                dimnames = list(mnams, c("IC-min","w")))

  print(round( t(out), 2))
  invisible(out)
}

##' @importFrom loo waic.matrix
WAIC_one <- function(object, mc_se = FALSE, ...)
{
  ## TODO: remove dependency on loo::waic.matrix, which is just some algebra via
  ## matrixStats
  ## ... for passing digits to print function
  ntrial <- dim(object)[1]
  nchain <- dim(object)[2]

  tmp0 <- t( matrix(object, nrow = ntrial) )
  out  <- loo::waic.matrix(tmp0)

  if (mc_se) {

    alist <- vector(mode="list", length=nchain)
    for (i in 1:nchain) alist[[i]] <- loo::waic.matrix(t(object[,i,]))
    # alist[[1]]$estimates
    #              Estimate          SE
    # elpd_waic  1086.04103 209.9144045
    # p_waic        6.37261   0.1980855
    # waic      -2172.08205 419.8288090
    tmp1 <- sapply(alist, function(x) {x$estimates[3,1]} )
    mc   <- sd(tmp1) / sqrt(nchain)

    print(out$estimates, ...)
    print(c(mc_se_waic = mc), ...)
  } else {
    print(out$estimates, ...)
  }
  invisible(out)
}

##' @importFrom loo loo.matrix
LOOIC_one <- function(object, mc_se = FALSE, ...)
{
  ## ... for passing digits to print function
  # object <- res0
  ntrial <- dim(object)[1]
  nchain <- dim(object)[2]

  tmp0 <- t( matrix(object, nrow = ntrial) )
  out  <- loo::loo.matrix(tmp0)
  pareto_k <- out$diagnostics$pareto_k

  ######Diagnosis####################40
  if (all(pareto_k < .5)) {
    message("All Pareto k estimates OK (k < 0.5)")
  } else {
    msg <- "See PSIS-LOO description (?'loo-package') for more information"

    if (any(pareto_k > 1)) {
      tmp1 <- sum(pareto_k > 1)
      tmp2 <- mean(pareto_k > 1)
      msg1 <- paste0(tmp1, " (", round(100*tmp2, 2),
                     "%) Pareto k estimates greater than 1\n")
    } else {
      tmp1 <- sum(pareto_k > .5)
      tmp2 <- mean(pareto_k > 5)
      msg1 <- paste0(tmp1," (", round(100*tmp2, 2),
                     "%) Pareto k estimates between 0.5 and 1\n")
    }
    warning(msg1, msg)
  }

  ############Not tested (could break)#####40
  if (mc_se) {
    alist <- vector(mode="list", length=nchain)
    for (i in 1:nchain) alist[[i]] <- loo::loo.matrix(t(object[,i,]))
    # alist[[1]]$estimates
    #              Estimate          SE
    # elpd_waic  1086.04103 209.9144045
    # p_waic        6.37261   0.1980855
    # waic      -2172.08205 419.8288090
    tmp1 <- sapply(alist, function(x) {x$estimates[3,1]} )
    mc   <- sd(tmp1) / sqrt(nchain)

    print(out$estimates, ...)
    print(c(mc_se_waic = mc), ...)
  } else {
    print(out$estimates, ...)
  }
  invisible(out)
}

##' @importFrom loo compare
compare <- function(loo1, loo2=NULL,...)
{
  if ( !is.null(loo2) ) {
    tmp <- loo::compare(loo1, loo2)
    out <- c(waic_diff = -2*tmp[[1]], se=2*tmp[[2]])
    print(out,...)
  } else {
    if ( !class(loo1)[1] =="list" ) {
      stop("Must bind multiple loo objects as a list")
    }

    ics <- sapply(loo1, function(x) {  x$estimates[3,1] })

    d <- ics - min(ics)
    w <- exp(-d/2) / sum(exp(-d/2))
    if (!is.null(names(ics))) {
      mnams <- names(ics)
    } else {
      mnams <- paste0("M", 1:length(d))
    }

    out <- t(matrix(c(d, w), nrow=length(d),
                    dimnames=list(mnams, c("IC-min","w"))))

    print( out, ...)
  }

  invisible(out)
}

posterior_predictive_test <- function(x, fun, hyper=FALSE, ptype=1, pnams=NA)
{
  #################################################################70
  ## This is a tmp tmp tmp function; not check yet for its accurayc
  #################################################################70
  ## AH's DMC explantion
  # Another way of testing if it is necessary to have a difference between v for
  # a1 and a2 is to do posterior predictive tests to see if the credible interval
  # for their sampled difference contains zero. Recall that in model 1 v.a1=.75
  # and v.a2=1.25, so we would expect v.a2-v.a1 = 0.5, whereas in model 0 they are
  # the same so we would expect no difference.

  # vdiff <- function(p) {p["v.a2"]-p["v.a1"]}
  #
  # vdiffs1 <- p.fun.dmc(samples=samples1.true, fun=vdiff, pnams=c("v.a2", "v.a1"))
  # round(c(mean(vdiffs1), quantile(vdiffs1,probs=c(.025,.975))),2)
  #        2.5% 97.5%
  #  0.47  0.42  0.52
  # vdiffs0 <- p.fun.dmc(samples=samples0.false,fun=vdiff, pnams=c("v.a2", "v.a1"))
  # round(c(mean(vdiffs0),quantile(vdiffs0,probs=c(.025,.975))),2)
  #        2.5% 97.5%
  #  0.01 -0.04  0.06

  if (!hyper) {

    if (any(is.na(pnams))) pnams <- dimnames(x$phi)[[2]]
    tmp0 <- apply(x$theta[,pnams,], c(1,3), fun)
    out <- as.vector(tmp0)

  } else {
    if (any(is.na(pnams))) pnams <- dimnames(attr(x, "hyper")$phi[[ptype]])[[2]]
    tmp0 <- apply(attr(x, "hyper")$phi[[ptype]][,pnams,],c(1,3),fun)
    out <- as.vector(tmp0)
  }
  return(out)
}
