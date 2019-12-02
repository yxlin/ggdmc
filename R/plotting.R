### One Subject  ---------------------------------------------
##' @import ggplot2
preplot_one <- function(x, start, end, pll) {

  # pll <- F
  # start <- 499

  if ( is.na(end) ) end <- slot(x, "nmc") ## This is neceaary, bcuz of plot_many
  nchain <- slot(x, "nchain")
  if (nchain == 1) stop ("MCMC needs multiple chains to check convergence")
  iter <- start:end

  if (pll) {
    lp <- x@summed_log_prior[,start:end] + x@log_likelihoods[,start:end]
    v <- lapply(1:nchain, function(k) {
      dd <- data.frame(Iteration = iter,
                       Chain     = k,
                       Parameter = "log-posterior",
                       value     = sapply(lp[k, ], c))
    })

  } else {
    theta  <- slot(x, "theta") ## npar x nchain x nmc
    pnames <- slot(x, "pnames")
    npar   <- slot(x, "npar")

    v <- lapply(1:nchain, function(k) {

      dd <- data.frame(Iteration = rep(iter, npar),
                       Chain     = k,
                       Parameter = rep(pnames, each = length(iter)),
                       value     = sapply(t (theta[,k,start:end]), c))
    })
  }

  DT <- data.table::rbindlist(v)

  DT$Parameter <- factor(DT$Parameter)
  DT$Chain     <- factor(DT$Chain)
  attr(DT, "nThin")       <- slot(x, "thin")
  attr(DT, "nIterations") <- end
  attr(DT, "nChains")     <- nchain
  attr(DT, "nParameters") <- slot(x, "npar")
  attr(DT, "pnames") <- slot(x, "pnames")
  return(DT)
}

##' @import ggplot2
plot_one <- function(x, start, end, pll, save, den, subchain, nsubchain,
  chains, ...) {

  if ( is.na(end) ) end <- slot(x, "nmc")
  if ( end <= start ) stop("End must be greater than start")
  DT <- preplot_one(x, start, end, pll)

  if (subchain) {
    if (missing(nsubchain)) stop("Please supply nsubchain")
    if (any(is.na(chains))) chains <- base::sample( unique(DT$Chain), nsubchain)
    cat("Plot chains: ", chains, "\n")
    DT <- DT[ DT$Chain %in% chains, ]
  }

  if (den)
  {
    f1 <- ggplot(DT, aes_string(x = "value", colour = "Chain", fill = "Chain")) +
      geom_density(alpha = 0.3) +
      scale_fill_discrete(name = "Chain") +
      scale_colour_discrete(name = "Chain") +
      facet_wrap(~Parameter, scales = "free") +
      xlab("") + ylab("Density") + theme(legend.position="none") +
      geom_rug(alpha = 0.1)
    print(f1)
  }
  else if (pll)
  {
    f1 <- ggplot(DT, aes_string(x = "Iteration", y = "value", color = "Chain")) +
      geom_line() +
      ylab("Posterior log-likelihood") +
      theme(legend.position = "none")
    print(f1)
  }
  else
  {
    f1 <- ggplot(DT, aes_string(x = "Iteration", y = "value", colour = "Chain")) +
      geom_line(alpha = 0.7) +
      scale_colour_discrete(name = "Chain") +
      facet_wrap(~Parameter, scales = "free") +
      ylab("") + theme(legend.position = "none")
    print(f1)
  }
  if (save) { return(DT) } else { return(f1) }
}

# autocor <- function(x, start = 1, end = NA, nLags = 50, nsubchain = NULL,
#                     base_size = 14)
# {
#
#   if (x$n.chains == 1) stop ("MCMC needs multiple chains to check convergence")
#   if (is.null(x$theta)) stop("Use hyper mcmc_list")
#   if ( is.na(end) ) end <- x$nmc
#   if ( end <= start ) stop("End must be greater than start")
#   if (!is.null(nsubchain)) idx <- sample(1:x$n.chains, nsubchain) else idx <- 1:x$n.chains
#
#   d <- ConvertChains(x, start, end, FALSE)
#   DT <- d[, .SD[, .(Lag = 1:nLags,
#                     Autocorrelation = ac_(value, nLags))],
#           .(Parameter, Chain)]
#
#   p0 <- ggplot(DT[Chain %in% idx],
#                aes(x = Lag, y = Autocorrelation, colour = Chain, fill = Chain)) +
#     geom_bar(stat = "identity", position = "identity") +
#     ylim(-1, 1) +
#     scale_fill_discrete(name = "Chain") +
#     scale_colour_discrete(name = "Chain") +
#     facet_grid(Parameter ~ Chain) +
#     theme_minimal(base_size = base_size) +
#     theme(legend.position = "none") +
#     theme(axis.text.x  = element_blank())
#   print(p0)
#   return(invisible(p0))
# }


### Many Subjects  ---------------------------------------------
##' @import ggplot2
plot_many <- function(x, start, end, pll, save, den, subchain, nsubchain,
  chains, ...) {

  nsub <- length(x)
  DT   <- lapply(x, preplot_one, start, end, pll)

  if (subchain) {
    if (missing(nsubchain)) stop("Please supply nsubchain")

    nchain <- slot(x[[1]], "nchain")

    if (any(is.na(chains))) chains <- base::sample(nchain, nsubchain)
    if (nsubchain != length(chains))
      message("nsubchain and chains are inconsistent. Ignore nsubchain")
    cat("Plot chains: ", chains, "\n")
  }

  if(is.null(names(DT))) {
    subjectnames <- 1:nsub
  } else {
    subjectnames <- names(DT)
  }

  x0 <- NULL
  for(i in 1:nsub) {
    tmp <- DT[[i]]
    if (subchain) tmp <- tmp[ tmp$Chain %in% chains, ]
    tmp$s <- factor(subjectnames[i])
    x0 <- rbind(x0, tmp)
  }

  if (pll) {

    f1 <- ggplot(x0, aes_string(x = "Iteration", y = "value", color = "Chain")) +
      geom_line() +
      facet_wrap(~s, scales = "free") +
      ylab("Log-posterior likelihood") +
      theme(legend.position = "none")
    print(f1)

  } else if (den) {

    f1 <- ggplot(x0, aes_string(x = "value", colour = "Chain", fill = "Chain")) +
      geom_density(alpha = 0.3) +
      scale_fill_discrete(name = "Chain") +
      scale_colour_discrete(name = "Chain") +
      facet_wrap(~s+Parameter, scales = "free") +
      xlab("") + ylab("Density") + theme(legend.position="none") +
      geom_rug(alpha = 0.1)
    print(f1)

  } else {
    f1 <- ggplot(x0, aes_string(x = "Iteration", y = "value", colour = "Chain")) +
      geom_line(alpha = 0.7) +
      scale_colour_discrete(name = "Chain") +
      facet_wrap(~s+Parameter, scales = "free") +
      ylab("") + theme(legend.position = "none")
    print(f1)
  }

  if (save) { return(x0) } else { return(f1) }
}

### Phi  ---------------------------------------------
preplot_phi <- function(x, start = 1, end = NA, pll = TRUE, ...) {
  # x <- fit
  # start <- 1
  # pll <- F
  # end <- NA
  # hyper <- attr(x, "hyper")
  if ( end <= start ) stop("End must be greater than start")
  nchain <- x@nchain
  iter   <- start:end

  if (pll) {
    lp <- x@summed_log_prior[,start:end] + x@log_likelihoods[,start:end]
    rownames(lp) <- 1:nchain
    v <- lapply(1:nchain, function(k) {
      dd <- data.frame(Iteration = iter,
                       Chain     = k,
                       Parameter = "log-posterior",
                       value     = sapply(lp[k, ], c))
    })

    DT <- data.table::rbindlist(v)

  } else {


    pnames   <- x@pnames
    npar     <- x@npar

    ## npar x nchain x nmc
    v0 <- lapply(1:nchain, function(k) {

      dd <- data.frame(Iteration = rep(iter, npar),
                       Chain     = k,
                       Parameter = rep(pnames, each = length(iter)),
                       Type      = "h1",
                       value     = sapply(t (x@phi_loc[,k,start:end]), c)
                       )
    })

    v1 <- lapply(1:nchain, function(k) {

      dd <- data.frame(Iteration = rep(iter, npar),
                       Chain     = k,
                       Parameter = rep(pnames, each = length(iter)),
                       Type      = "h2",
                       value     = sapply(t (x@phi_sca[,k,start:end]), c)
      )
    })

    tmp0 <- data.table::rbindlist(v0)
    tmp1 <- data.table::rbindlist(v1)
    DT <- rbind(tmp0, tmp1)
    DT$Type <- factor(DT$Type)
  }

  DT$Parameter <- factor(DT$Parameter)
  DT$Chain     <- factor(DT$Chain)
  attr(DT, "nThin")       <- x@thin
  attr(DT, "nIterations") <- end
  attr(DT, "nChains")     <- nchain
  attr(DT, "nParameters") <- x@npar
  attr(DT, "pnames") <- x@pnames
  return(DT)
}


##' @import ggplot2
plot_phi <- function(x, start, end, pll, save, den, subchain, nsubchain,
  chains, ...) {

  DT <- preplot_phi(x, start, end, pll)

  if (subchain) {
    if (missing(nsubchain)) stop("Please supply nsubchain")
    if (any(is.na(chains))) chains <- sample(unique(DT$Chain), nsubchain)
    cat("Plot chains: ", chains, "\n")
    DT <- DT[ DT$Chain %in% chains, ]
  }


  if (pll) {
    # DT$Parameter <- "lp"
    ## Output 1
    f1 <- ggplot(DT) +
      geom_line(aes_string(x = "Iteration", y = "value", color = "Chain")) +
      ylab("Posterior log-likelihood") +
      theme_bw(base_size = 14) +
      theme(legend.position = "none")
    print(f1)

  } else if (den) {

    f1 <- ggplot(DT, aes_string(x = "value", colour = "Chain", fill = "Chain")) +
      geom_density(alpha = 0.3) +
      scale_fill_discrete(name = "Chain") +
      scale_colour_discrete(name = "Chain") +
      facet_wrap(~Parameter+Type, scales = "free") +
      xlab("") + ylab("Density") + theme(legend.position="none") +
      geom_rug(alpha = 0.1)
    print(f1)

  } else {
    f1 <- ggplot(DT, aes_string(x = "Iteration", y = "value", colour = "Chain")) +
      geom_line(alpha = 0.7) +
      scale_colour_discrete(name = "Chain") +
      facet_wrap(~Parameter+Type, scales = "free") +
      ylab("") + theme(legend.position = "none")
    print(f1)

  }
  if (save) { return(DT) } else { return(f1) }
}


### Prior  ---------------------------------------------
plot_prior <- function(i, prior, xlim = NA, natural = TRUE, npoint = 100,
                       trans = NA, save = FALSE, ... ) {

  pnames <- slot(prior, "pnames")
  priors <- slot(prior, "priors")

  ### Check 1 ###
  if (is.numeric(i)) i <- pnames[i]
  if (any(is.na(trans))) trans <- ifelse(natural, attr(priors[[i]], "untrans"), "identity")
  if ( !(i %in% pnames) ) stop("Parameter not in prior")

  ## Finish checks. Calculate the data frame
  ## Save all parameter setting to another object called p
  p    <- priors[[i]]
  dist <- attr(p, "dist")

  if ( is.na(dist) | dist == 0 ) {
    disttype <- NA
  } else {
    disttype <- switch(dist, "tnorm", "beta_lu", "gamma_l", "lnorm_l",
                       "unif_", "constant", "tnorm2")
  }

  ## Do an educated guess for xlim
  ## xlim can be a scalar or vector. test if any of its elements is NA. If it
  if ( any(is.na(xlim)) ) {
    ## does contain NAs, we set a xlim for the user
    ## check what distribution the user wants to as her prior distribution
    ## 4 options: tnorm (default), beta_lu, gamma_l, and lnorm_l
    ## Basically, we use parameter 1 and 2 to do some (random?) arithmetic
    ## and pick (1) the bigger and smaller number (tnorm); (2) use lower and
    ## upper directly (beta_lu); (3) use lower and a different arithmetic
    ## (gamma_l); (4) use lower and another different arithmetic (lnorm_l)
    x <- switch(disttype,
                ## Now we get xlim. Then we want to set all the x points for
                ##  plotting. By default we plot 100 point (n.point = 1e2)
                tnorm    = {
                     xlim <- c(pmax(p$lower, p[[1]] - 3*p[[2]]),
                       pmin(p$upper, p[[1]] + 3*p[[2]]))
                     seq(xlim[1], xlim[2], length.out = npoint)
                },
                beta_lu  = {
                     xlim <- c(p$lower, p$upper)
                     seq(xlim[1], xlim[2], length.out = npoint)
                },
                gamma_l  = {
                     xlim <- c(p$lower, p[[1]]*p[[2]] + 3*sqrt(p[[1]])*p[[2]])
                     seq(xlim[1], xlim[2], length.out = npoint)
                },
                lnorm_l  = {
                     xlim <- c(p$lower, exp(p[[1]]+2*p[[2]]))
                     seq(xlim[1], xlim[2], length.out = npoint)
                },
                unif_    = {
                     xlim <- c(p[[1]], p[[2]])
                     seq(xlim[1], xlim[2], length.out = npoint)
                },
                constant = {
                     xlim <- c(p[[1]] - p[[2]], p[[1]] + p[[2]])
                     tmp <- seq(xlim[1], xlim[2], length.out = npoint)
                     sort(c(p[[1]], tmp))
                },
                tnorm2  = {
                  sd <- 1/sqrt(p[[2]])
                  xlim <- c(pmax(p$lower, p[[1]] - 3 * sd),
                            pmin(p$upper, p[[1]] + 3 * sd))
                  seq(xlim[1], xlim[2], length.out = npoint)
                },
                {
                     c(-1, 1)
                }
    )
  } else {
    x <- c(-1, 1)
  }
  p$x <- x
  p$lg <- FALSE    ## Turn off logarithmic transformation

  ## 1. call a log/identity transform function to change x value
  ## 2. call a density function for a distribution to set y (density) value
  if (is.na(disttype) | disttype == 0) {
    message(paste0(i, " prior not available"))
    invisible(NULL)
  } else if (save) {
    xpos <- do.call(trans, list(x = x))
    ypos <- do.call(paste0("d", disttype), p)
    dat <- data.frame(xpos  = xpos, ypos  = ypos,
                      Parameter = rep(names(priors[i]), length(x)))
    return(dat)
  } else {
    xpos <- do.call(trans, list(x = x))
    ypos <- do.call(paste0("d", disttype), p)
    dat <- data.frame(xpos  = xpos, ypos  = ypos,
                      Parameter = rep(names(priors[i]), length(x)))

    p0 <- ggplot(dat, aes_string(x="xpos", y="ypos")) + geom_line() +
      xlab("")+ ylab("Density") + facet_grid(. ~Parameter) +
      theme_bw(base_size = 16)
    print(p0)
    invisible(p0)
  }

}

### Autocorrelation ---------------------------------------------
##' Autocorrelation Plot
##'
##' Plot the autocorrelation of posterior samples,
##'
##' @param x posterior samples
##' @param start start from which iteration.
##' @param end end at which iteration
##' @param nLags the maximum number of lags.
##' @param subchain a Boolean switch to plot a subset of chains.
##' @param pll a Boolean switch for plotting parameter values or posterior
##' log likelihoods
##' @examples
##' ## Model 1
##' ## 27 elements with 20 levels
##' FR <- list(S = c("n","w","p"), cond=c("C","F", "H"), R=c("N", "W", "P"))
##' lev <- c("CnN","CwN", "CnW","CwW",
##'          "FnN","FwN","FpN", "FnW","FwW","FpW", "fa","FpP",
##'          "HnN","HwN","HpN", "HnW","HwW","HpW", "HpP",
##'          "FAKERATE")
##' map_mean_v <- ggdmc:::MakeEmptyMap(FR, lev)
##' map_mean_v[1:27] <- c(
##'   "CnN","CwN","FAKERATE", "FnN","FwN","FpN", "HnN","HwN","HpN",
##'   "CnW","CwW","FAKERATE", "FnW","FwW","FpW", "HnW","HwW","HpW",
##'   "FAKERATE","FAKERATE","FAKERATE", "fa","fa","FpP", "fa","fa","HpP")
##'
##' model0 <- BuildModel(
##'   p.map     = list(A = "1", B = c("cond", "R"), t0 = "1", mean_v = c("MAPMV"),
##'                    sd_v = "1", st0 = "1", N = "cond"),
##'   match.map = list(M = list(n = "N", w = "W", p = "P"), MAPMV = map_mean_v),
##'   factors   = list(S = c("n","w","p"), cond = c("C","F", "H")),
##'   constants = c(N.C = 2, N.F = 3, N.H = 3, st0 = 0, B.C.P = Inf,
##'                 mean_v.FAKERATE = 1, sd_v = 1),
##'   responses = c("N", "W", "P"),
##'   type      = "norm")
##'
##' npar <- model0@npar
##'
##' p.vector <- c(A = .3, B.C.N = 1.3,  B.F.N = 1.3,  B.H.N = 1.3,
##'               B.C.W = 1.3,  B.F.W = 1.4,  B.H.W = 1.5,
##'               B.F.P = 1.1,  B.H.P = 1.3,
##'
##'               t0=.1,
##'
##'               mean_v.CnN = 2.8,  mean_v.CwN = -0.3, mean_v.CnW=-1,
##'               mean_v.CwW = 2.9,  mean_v.FnN = 2.8,  mean_v.FwN=-.3,
##'
##'               mean_v.FpN = -1.6, mean_v.FnW = -1,   mean_v.FwW = 2.9,
##'               mean_v.FpW = .5 ,  mean_v.fa = -2.4,  mean_v.FpP = 2.5,
##'
##'               mean_v.HnN = 2.8, mean_v.HwN = -.5,   mean_v.HpN = -.6,
##'               mean_v.HnW = -.7, mean_v.HwW = 3.0,   mean_v.HpW = 1.6,
##'               mean_v.HpP = 2.3)
##'
##' acc_tab0 <- TableParameters(p.vector, 1, model0, FALSE)
##' acc_tab1 <- TableParameters(p.vector, "w.C.N", model0, FALSE)
##' acc_tab2 <- TableParameters(p.vector, "w.F.P", model0, FALSE)
##' print(acc_tab0); print(acc_tab1); print(acc_tab2)
##'
##' \dontrun{
##' dat0 <- simulate(model0, nsim=50, ps=p.vector)
##' dmi0 <- BuildDMI(dat0, model0)
##' }
##' p1 <- rep(1, npar)
##' names(p1) <- model0@pnames
##'
##' p.prior0 <- BuildPrior(
##'   dists = c(rep("tnorm", 9), "beta", rep("tnorm", 19)),
##'   p1    = p1,
##'   p2    = c(rep(2, 9), 1, rep(2, 19)),
##'   lower = c(rep(0, 10),  rep(NA, 19)),
##'   upper = c(rep(NA, 9), 1, rep(NA, 19)))
##'
##' # plot(p.prior0, ps = p.vector)
##' ## Sampling
##' ## 18.4 & 36.17 s
##' \dontrun{
##' fit0 <- StartNewsamples(dmi0, p.prior0, block = FALSE, thin=4)
##' fit0_correct  <- run(fit0, thin=4, block = FALSE)
##'
##' hat <- gelman(fit0_correct, verbose=TRUE);
##' p0 <- autocorr(fit0_correct, subchain=1:3, pll=TRUE)
##' }
##' @export
autocorr <- function(x, start = 1, end = NA, nLags = 50, pll = TRUE,
                     subchain = FALSE)
{
  if (x@nchain == 1) stop ("MCMC needs multiple chains to check convergence")
  if (is.null(x@theta)) stop("autocorr plot only posterior class")
  if (is.na(end)) end <- x@nmc
  if ( end <= start ) stop("End must be greater than start")
  if (end < nLags)
  {
    warning(sprintf("nLags = %d, which is larger than number of iterations, computing until max possible lag %d", nLags, end))
    nLags <- end
  }

  DT <- preplot_one(x, start, end, pll)
  if (any(subchain)) {
    nsubchain <- length(subchain)
    # chains <- base::sample(unique(DT$Chain), nsubchain)
    cat("Plot chains: ", subchain, "\n")
    DT <- DT[ DT$Chain %in% subchain, ]
    DT$Chain <- droplevels(DT$Chain)
  }

  D_ac <- DT[, .(Lag = ac(value, nLags)[[1]],
                 Autocorrelation = ac(value, nLags)[[2]]), .(Parameter, Chain)]

  if ( length(unique(D_ac$Chain)) <= 1) {
    p0 <- ggplot(D_ac, aes(x = Lag, y = Autocorrelation)) +
      geom_bar(stat = "identity", position = "identity") +
      facet_wrap(~Parameter) +
      theme_bw(base_size = 14) +
      theme(legend.position="none")

  } else {
    p0 <- ggplot(D_ac, aes(x = Lag, y = Autocorrelation, colour = Chain,
                           fill = Chain)) +
      geom_bar(stat = "identity", position = "identity") +
      scale_fill_discrete(name = "Chain") +
      scale_colour_discrete(name = "Chain") +
      facet_grid(Parameter ~ Chain) +
      theme_bw(base_size = 14) +
      theme(legend.position="none")
  }

  print(p0)
}

### Class Methods ---------------------------------------------------------
##' ggdmc Plotting Methods
##'
##' The function plots prior distributions or posterior samples depending on
##' whether the first argument \code{x} is a prior object or an object
##' storing posterior samples.
##'
##' @param x a prior object or posterior samples.
##' @param y NULL
##' @param hyper a Boolean switch, indicating posterior samples are from
##' hierarchical modeling
##' @param ps a parameter vector
##' @param start start from iteration
##' @param end end at which iteraton
##' @param pll a Boolean switch whether to plot posterior log likelihoods
##' @param den a Boolean switch whether for density plots
##' @param subchain a Boolean switch whether to plot a subset of chains.
##' @param nsubchain number of subchain
##' @param chains indicate the subchains to plot. This must be an integer vector
##' @param save a Boolean switch whether to save plotting data
##' @param ... Additional argument passing via dot dot dot.
##'
##' @return NULL
##'
##' @export
##' @docType methods
##' @rdname plot-methods
##' @examples
##' p.prior <- BuildPrior(
##'            dists = rep("tnorm", 7),
##'            p1    = c(a = 2,   v.f1 = 4,  v.f2 = 3,  z = 0.5, sv = 1,
##'                      sz = 0.3, t0 = 0.3),
##'            p2    = c(a = 0.5, v.f1 = .5, v.f2 = .5, z = 0.1, sv = .3,
##'                      sz = 0.1, t0 = 0.05),
##'            lower = c(0,-5, -5, 0, 0, 0, 0),
##'            upper = c(5, 7,  7, 1, 2, 1, 1))
##' plot(p.prior)
setGeneric("plot", function(x, y = NULL, ...) {
  standardGeneric("plot")
})

##' @import ggplot2
##' @importFrom data.table data.table melt.data.table
##' @rdname plot-methods
setMethod("plot", "prior", function (x, y = NULL, ps = NULL, save = FALSE, ...)
{
  pnames <- slot(x, "pnames")
  if( is.null(pnames) ) stop("Prior object must have name attributes.")

  pD <- NULL
  for(j in pnames) {
    invisible(pD <- rbind(pD, plot_prior(i = j, prior = x, save = TRUE)))
  }

  if (save) {
    return(pD)
  } else if (!is.null(ps) & is.vector(ps)) {
    pveclines <- data.frame(Parameter = pnames, true = ps)
    p0 <- ggplot(pD, aes_string(x = "xpos", y = "ypos")) +
      geom_line() +
      geom_vline(data = pveclines, aes_string(xintercept = "true"),
                 linetype = "dotted", size = 1) +
      xlab("")+ ylab("")+
      theme_bw(base_size = 14) +
      facet_wrap( ~ Parameter, scales = "free")
    print(p0)
    invisible(return(p0))
  } else if ((!is.null(ps) & is.matrix(ps))) {
    wide <- data.table::data.table(ps)
    wide$s <- factor(1:nrow(ps))
    pveclines <- data.table::melt.data.table(
      wide, id.vars = "s", variable.name = "Parameter", value.name = "true")
    p0 <- ggplot(pD, aes_string(x = "xpos", y = "ypos")) +
      geom_line() +
      geom_vline(data = pveclines, aes_string(xintercept = "true"),
                 linetype = "dotted", size = 1) +
      xlab("")+ ylab("")+
      theme_bw(base_size = 14) +
      facet_wrap( ~ Parameter, scales = "free")
    print(p0)
    invisible(return(p0))
  } else {
    p0 <- ggplot(pD, aes_string(x = "xpos", y = "ypos")) +
      geom_line() +
      xlab("")+ ylab("")+ theme_bw(base_size = 14) +
      facet_wrap( ~ Parameter, scales = "free")
    print(p0)
    invisible(return(p0))
  }
})

##' @rdname plot-methods
setMethod("plot", "posterior", function (x, y = NULL, hyper = FALSE, start = 1,
                                         end = NA, pll = TRUE, save = FALSE,
                                         den = FALSE, subchain = FALSE,
                                         nsubchain = 3, chains = NA, ...) {
  if (is.na(end)) end <- x@nmc
  out <- plot_one(x, start, end, pll, save, den, subchain, nsubchain,
                  chains)
  return(out)
})

##' @rdname plot-methods
setMethod("plot", "hyper", function (x, y = NULL, hyper = TRUE, start = 1,
                                     end = NA, pll = TRUE, save = FALSE,
                                     den = FALSE, subchain = FALSE,
                                     nsubchain = 3, chains = NA, ...) {

  if (is.na(end)) end <- x@nmc

  if (hyper) {
    out <- plot_phi(x, start, end, pll, save, den, subchain, nsubchain,
                    chains)
  } else {

    out <- plot_many(x@individuals, start, end, pll, save, den,
                     subchain, nsubchain, chains)
  }
  return(out)
})

##' @rdname plot-methods
setMethod("plot", "list", function (x, y = NULL, start = 1,
                                    end = NA, pll = TRUE, save = FALSE,
                                    den = FALSE, subchain = FALSE,
                                    nsubchain = 3, chains = NA, ...)
{
  out <- plot_many(x, start, end, pll, save, den, subchain, nsubchain,
                   chains)
})

# pairs.model <- function(x, start = 1, end = NA, ...) {
#   ## Make this function invisible;
#   if (x$n.chains == 1) stop ("MCMC needs multiple chains to check convergence")
#   if (is.null(x$theta)) stop("Use hyper mcmc_list")
#   if (is.na(end)) end <- x$nmc
#   if (end <= start) stop("End must be greater than start")
#
#   d <- ConvertChains(x, start, end, FALSE)
#   D_wide <- data.table::dcast.data.table(d, Iteration + Chain ~ Parameter, value.var = "value")
#
#   bracket_names <- names(D_wide)
#   par_cols <- !(bracket_names %in% c("Iteration", "Chain"))
#   p0 <- GGally::ggpairs(D_wide, columnLabels = bracket_names[par_cols],
#                         columns = which(par_cols), ...)
#   print(p0)
#   return(invisible(p0))
# }


