### Sampling -------------------------------------------------------------------

##' @importFrom parallel parLapply
##' @importFrom parallel makeCluster
##' @importFrom parallel mclapply
##' @importFrom parallel stopCluster
run_many <- function(dmi, prior, nchain, nmc, thin, report, rp, gammamult,
                     pm0, pm1, block, ncore)
{

  if (get_os() == "windows" & ncore > 1)
  {
    cl  <- parallel::makeCluster(ncore)
    message("fits multi-participant in parallel")
    out <- parallel::parLapply(cl = cl, X = dmi,
                               init_new, prior, nchain, nmc, thin, report,
                               rp, gammamult, pm0, pm1, block)
    parallel::stopCluster(cl)
  }
  else if (ncore > 1)
  {
    message("fits multi-participant in parallel")
    out <- parallel::mclapply(dmi, init_new, prior, nchain, nmc, thin,
                                    report, rp, gammamult, pm0, pm1, block,
                              mc.cores=getOption("mc.cores", ncore))
  }
  else
  {
    message("fits multi-participant with lapply")
    out <- lapply(dmi, init_new, prior, nchain, nmc, thin, report, rp,
                  gammamult, pm0, pm1, block)
  }

  return(out)
}

##' @importFrom parallel parLapply
##' @importFrom parallel makeCluster
##' @importFrom parallel mclapply
##' @importFrom parallel stopCluster
rerun_many <- function(samples, nmc, thin, report, rp, gammamult, pm0, pm1,
                       block, add, ncore)
{
  message("Fit multi-participant")

  if (get_os() == "windows" & ncore > 1)
  {
    cl  <- parallel::makeCluster(ncore)
    out <- parallel::parLapply(cl = cl, X = samples,
                               init_old, nmc, thin, report, rp,
                               gammamult, pm0, pm1, block, add)
    parallel::stopCluster(cl)
  }
  else if (ncore > 1)
  {
    out <- parallel::mclapply(samples, init_old, nmc, thin, report, rp,
                              gammamult, pm0, pm1, block, add,
                              mc.cores=getOption("mc.cores", ncore))
  }
  else
  {
    out <- lapply(samples, init_old, nmc, thin, report, rp,
                  gammamult, pm0, pm1, block, add)
  }

  return(out)
}

##' Start new model fits
##'
##' Fit a hierarchical or a fixed-effect model, using Bayeisan
##' optimisation.  We use a specific type of pMCMC algorithm, the DE-MCMC. This
##' particular sampling method includes crossover and two different migration
##' operators. The migration operators are similar to random-walk algorithm.
##' They wouold be less efficient to find the target parameter space, if been
##' used alone.
##'
##' @param dmi a data model instance or a list of data model instances
##' @param samples posterior samples.
##' @param prior prior objects.  For hierarchical model, this must be a
##' list with three sets of prior distributions. Each is respectively named,
##' "pprior", "location", and "scale".
##' @param nmc number of Monte Carlo samples
##' @param thin thinning length
##' @param nchain number of chains
##' @param report progress report interval
##' @param rp tuning parameter 1
##' @param gammamult tuning parameter 2. This is the step size.
##' @param pm0 probability of migration type 0 (Hu & Tsui, 2010)
##' @param pm1 probability of migration type 1 (Turner et al., 2013)
##' @param block Only for hierarchical modeling. A Boolean switch for update one
##' parameter at a time
##' @param ncore Only for non-hierarchical, fixed-effect models with many
##' subjects.
##' @param add Boolean whether to add new samples
##'
##' @export
StartNewsamples <- function(dmi, prior, nmc=2e2, thin=1, nchain=NULL,
                            report=1e2, rp=.001, gammamult=2.38, pm0=.05,
                            pm1=.05, block=TRUE, ncore=1)
{
  if ( !missingArg(prior) &&
       all(c("pprior", "location", "scale") %in% names(prior)) &&
       length(prior) == 3 )
  {
    nchain <- CheckHyperDMI(dmi, prior, nchain)
    tmp0   <- sapply(dmi, checklba)

    message("Hierarchical model: ", appendLF = FALSE)
    t0 <- Sys.time()
    out <- init_newhier(prior[[1]], prior[[2]], prior[[3]], dmi, nchain, nmc,
                        thin, report, rp, gammamult, pm0, pm1, block)
    t1 <- Sys.time()
  } else if ( class(dmi) == "list" ) {
    nchain <- sapply(dmi, CheckDMI, prior=prior, nchain=nchain)
    nchain <- unique(nchain)
    tmp0   <- sapply(dmi, checklba)

    message("Fixed-effect model ", appendLF=FALSE)
    t0 <- Sys.time()
    out <- run_many(dmi, prior, nchain, nmc, thin, report, rp, gammamult, pm0,
                    pm1, block, ncore)
    t1 <- Sys.time()

  } else if (class(dmi) == "dmi") {
    nchain <- CheckDMI(dmi, prior, nchain)
    if(is.null(nchain)) nchain <- 3*slot(prior, "npar")
    checklba(dmi)

    message("Fixed-effect model (ncore has no effect): ", appendLF = FALSE)
    t0 <- Sys.time()
    out <- init_new(dmi, prior, nchain, nmc, thin, report, rp, gammamult, pm0,
                    pm1, block)
    t1 <- Sys.time()
  } else {
    stop("Class undefined")
  }

  proc_time <- difftime(t1, t0, units = "secs")[[1]]
  message("Processing time: ", round(proc_time, 2), " secs.")

  return(out)
}

##' @rdname StartNewsamples
##' @export
run <-  function(samples, nmc=5e2, thin=1, report=1e2, rp=.001,
                 gammamult=2.38, pm0=0, pm1=0, block=TRUE, ncore=1,
                 add=FALSE)
{

  if ( class(samples) == "hyper" )
  {
    t0 <- Sys.time()
    out <- init_oldhier(samples, nmc, thin, report, rp, gammamult, pm0, pm1,
                        block, add)
    t1 <- Sys.time()

  }
  else if ( is.list(samples) )
  {
    t0 <- Sys.time()
    out <- rerun_many(samples, nmc, thin, report, rp, gammamult, pm0, pm1,
                      block, add, ncore)
    t1 <- Sys.time()
    names(out) <- names(samples)

  }
  else
  {
    t0 <- Sys.time()
    out <- init_old(samples, nmc, thin, report, rp, gammamult, pm0, pm1,
                    block, add)
    t1 <- Sys.time()
  }

  proc_time <- difftime(t1, t0, units = "secs")[[1]]
  message("Processing time: ", round(proc_time, 2), " secs.")

  return(out)
}

## Utilities ------------------------------------------------------------------
CheckDMI <- function(d = NULL, prior = NULL, nchain=NULL)
{
  # data <- dmi
  # prior <- p.prior
  if (is.null(d)) stop("No data model instance")
  if (!is.data.frame(d@data)) stop("Data must be a data frame")

  model  <- d@model
  priors <- prior@priors
  npar   <- model@npar

  if (is.null(model))  stop("Must specify a model")
  if (is.null(priors)) stop("Must specify prior distributions")
  if (prior@npar != npar)
  {
    stop("data prior incompatible with model")
  }

  if (is.null(nchain)) nchain <- 3*npar
  return(nchain)
}

CheckHyperDMI <- function(data, prior = NULL, nchain = NULL)
{
  ## data
  if (missingArg(data)) stop("No data-model instance")
  if (is.data.frame(data)) stop("data must not be a data.frame")

  ## Prior
  for (i in 1:3) if (is.null(prior[[i]])) stop("No prior distribution")
  if ( prior[[1]]@npar < prior[[2]]@npar )
    stop("Location and scale priors differ in the numbers of parameters")

  model1 <- data[[1]]@model
  pnames <- model1@pnames

  isppriorok <- pnames %in% prior[[1]]@pnames
  islpriorok <- pnames %in% prior[[2]]@pnames
  isspriorok <- pnames %in% prior[[3]]@pnames

  if (!all(isppriorok))  {
    cat("Here is the parameter in your model\n")
    print(pnames)
    cat("Here is the parameter in your data prior\n")
    print(slot(prior[[1]], "pnames"))
    stop("data prior incompatible with model")
  }
  if (!all(islpriorok)) stop("location prior incompatible with model")
  if (!all(isppriorok)) stop("scale prior incompatible with model")

  ## nchain
  if (is.null(nchain)) {
    nchain <- 3*length(pnames)
    message("Use ", nchain, " chains")
  }

  return(nchain)
}
