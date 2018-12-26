##' Autocorrelation Plot
##'
##' Plot the autocorrelation of posterior samples,
##'
##' @param x posterior samples
##' @param start start from which iteration.
##' @param end end at which iteration
##' @param nLags the number of lags of the autocorrelation plot.
##' @param nsubchain plot only a chain subset
##' @export
autocor <- function(x, start = 1, end = NA, nLags = 50, nsubchain = NULL) {

  if (x$n.chains == 1) stop ("MCMC needs multiple chains to check convergence")
  if (is.null(x$theta)) stop("Use hyper mcmc_list")
  if ( is.na(end) ) end <- x$nmc
  if ( end <= start ) stop("End must be greater than start")
  if (!is.null(nsubchain)) idx <- sample(1:x$n.chains, nsubchain) else idx <- 1:x$n.chains

  d <- ConvertChains(x, start, end, FALSE)
  DT <- d[, .SD[, .(Lag = 1:nLags,
                    Autocorrelation = ggdmc:::ac_(value, nLags))],
              .(Parameter, Chain)]

  p0 <- ggplot(DT[Chain %in% idx],
               aes(x = Lag, y = Autocorrelation, colour = Chain, fill = Chain)) +
    geom_bar(stat = "identity", position = "identity") +
    ylim(-1, 1) +
    scale_fill_discrete(name = "Chain") +
    scale_colour_discrete(name = "Chain") +
    facet_grid(Parameter ~ Chain) +
    theme(legend.position = "none") +
    theme(axis.text.x  = element_blank())
  print(p0)
  invisible(p0)
}

##' Correlation Matrix
##'
##' Plot the correlation matrix of posterior samples,
##'
##' @param x posterior samples
##' @param start start from which iteration.
##' @param end end at which iteration.
##' @param ... other arguments
##' @importFrom ggmcmc ggs_pairs
##' @export
pairs.model <- function(x, start = 1, end = NA, ...) {

  if (x$n.chains == 1) stop ("MCMC needs multiple chains to check convergence")
  if (is.null(x$theta)) stop("Use hyper mcmc_list")
  if ( is.na(end) ) end <- x$nmc
  if ( end <= start ) stop("End must be greater than start")

  d <- ConvertChains(x, start, end, FALSE)
  ggs_pairs(d, lower = list(continuous = "density"))

}

##' @export
##' @importFrom ggmcmc ggs
plot.model <- function(x, y = NULL, hyper = FALSE, start = 1,
  end = NA, pll = TRUE, save = FALSE, den = FALSE, subchain = FALSE,
  nsubchain = 3, chains = NA, ...) {


  if (hyper) {
    out <- plot_phi(x, start, end, pll, save, den, subchain, nsubchain,
      chains)
  } else if (!is.null(x$theta)) {
    ## single subject
    out <- plot_one(x, start, end, pll, save, den, subchain, nsubchain,
      chains)
  } else {
    out <- plot_many(x, start, end, pll, save, den, subchain, nsubchain,
      chains)
  }
  return(out)
}


##' @import ggplot2
##' @importFrom coda mcmc mcmc.list
preplot_one <- function(x, start, end, pll) {

  if ( is.na(end) ) end <- x$nmc
  if ( end <= start ) stop("End must be greater than start")
  nchain <- x$n.chains
  # thin   <- x$thin

  if (pll) {
    lp <- x$summed_log_prior[start:end,] + x$log_likelihoods[start:end,]
    colnames(lp) <- 1:nchain
    step1 <- lapply(data.frame(lp), function(xx){
      coda::mcmc(as.matrix(xx), start, end, thin = 1)
    })
    d <- coda::mcmc.list(step1) ## log-posterior likelihood

    attr(d, "nchain") <- nchain
    attr(d, "npar")   <- x$n.pars
    attr(d, "thin")   <- x$thin
    attr(d, "iter")   <- start:end
    attr(d, "pnames") <- x$p.names
    attr(d, "nmc")    <- x$nmc
    attr(d, "start")  <- start
    attr(d, "end")    <- end

  } else {
    d <- theta2mcmclist(x, start, end, thin = 1)
  }

  DT <- ConvertChains2(d, pll)
  # DT$Chain <- factor(DT$Chain)
  # if (pll) DT$Parameter <- "lp"
  return(DT)
}

##' @import ggplot2
plot_one <- function(x, start, end, pll, save, den, subchain, nsubchain,
  chains, ...) {

  if ( is.na(end) ) end <- x$nmc
  if ( end <= start ) stop("End must be greater than start")
  DT <- preplot_one(x, start, end, pll)

  if (subchain) {
    if (missing(nsubchain)) stop("Please supply nsubchain")
    if (any(is.na(chains))) chains <- base::sample( unique(DT$Chain), nsubchain)
    cat("Plot chains: ", chains, "\n")
    DT <- DT[ DT$Chain %in% chains, ]
  }

  if (pll) {
    DT$Parameter <- "lp"

    f1 <- ggplot(DT, aes_string(x = "Iteration", y = "value", color = "Chain")) +
      geom_line() +
      ylab("Posterior log-likelihood") +
      theme(legend.position = "none")
    print(f1)

  } else if (den) {

    f1 <- ggplot(DT, aes_string(x = "value", colour = "Chain", fill = "Chain")) +
      geom_density(alpha = 0.3) +
      scale_fill_discrete(name = "Chain") +
      scale_colour_discrete(name = "Chain") +
      facet_wrap(~Parameter, scales = "free") +
      xlab("") + ylab("Density") + theme(legend.position="none") +
      geom_rug(alpha = 0.1)
    print(f1)

  } else {
    f1 <- ggplot(DT, aes_string(x = "Iteration", y = "value", colour = "Chain")) +
      geom_line(alpha = 0.7) +
      scale_colour_discrete(name = "Chain") +
      facet_wrap(~Parameter, scales = "free") +
      ylab("") + theme(legend.position = "none")
    print(f1)
  }
  if (save) { return(DT) } else { return(f1) }
}

##' @import ggplot2
preplot_many <- function(x, start = 1, end = NA, pll = TRUE) {

  xx <- lapply(x, preplot_one, start, end, pll)
  return(xx)
}

# x <- fit
# start <- 1
# end <- 500
# pll <- T
# save <- F
# den <- F
# subchain <- F
# nsubchain <- 3
##' @import ggplot2
plot_many <- function(x, start, end, pll, save, den, subchain, nsubchain,
  chains, ...) {

  nsub <- length(x)
  DT <- ggdmc:::preplot_many(x, start, end, pll)

  if (subchain) {
    if (missing(nsubchain)) stop("Please supply nsubchain")
    if (any(is.na(chains))) chains <- base::sample(x[[1]]$n.chains, nsubchain)
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

##' @importFrom ggmcmc ggs
preplot_phi <- function(x, start = 1, end = NA, pll = TRUE, ...) {
  # x <- fit1
  # start <- 1
  # pll <- FALSE
  xx <- attr(x, "hyper")
  if ( is.na(end) ) end <- xx$nmc
  if ( end <= start ) stop("End must be greater than start")
  nchain <- xx$n.chains
  # thin   <- xx$thin

  if (pll) {
    lp <- xx$h_summed_log_prior[start:end,] + xx$h_log_likelihoods[start:end,]
    colnames(lp) <- 1:nchain
    step1 <- lapply(data.frame(lp), function(xxx){
      coda::mcmc(as.matrix(xxx), start, end, thin = 1) ## thin must be 1 for coda
    })
    d <- coda::mcmc.list(step1) ## log-posterior likelihood
    attr(d, "nchain") <- nchain
    attr(d, "npar")   <- xx$n.pars
    attr(d, "thin")   <- xx$thin
    attr(d, "iter")   <- start:end
    attr(d, "pnames") <- xx$p.names
    attr(d, "nmc")    <- xx$nmc
    attr(d, "start")  <- start
    attr(d, "end")    <- end

  } else {

    d <- phi2mcmclist(xx, start, end)
  }

  DT <- ConvertChains2(d, pll)
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
    DT$Parameter <- "lp"
    ## Output 1
    f1 <- ggplot(DT) +
      geom_line(aes_string(x = "Iteration", y = "value", color = "Chain")) +
      ylab("Log-posterior likelihood") +
      theme(legend.position = "none")
    print(f1)
  } else if (den) {

    f1 <- ggplot(DT, aes_string(x = "value", colour = "Chain", fill = "Chain")) +
      geom_density(alpha = 0.3) +
      scale_fill_discrete(name = "Chain") +
      scale_colour_discrete(name = "Chain") +
      facet_wrap(~Parameter, scales = "free") +
      xlab("") + ylab("Density") + theme(legend.position="none") +
      geom_rug(alpha = 0.1)
    print(f1)

  } else {
    f1 <- ggplot(DT, aes_string(x = "Iteration", y = "value", colour = "Chain")) +
      geom_line(alpha = 0.7) +
      scale_colour_discrete(name = "Chain") +
      facet_wrap(~Parameter, scales = "free") +
      ylab("") + theme(legend.position = "none")
    print(f1)

  }
  if (save) { return(DT) } else { return(f1) }
}


#' @import ggplot2
#' @importFrom graphics plot
plot_subchain <- function(x, nchain, hyper = FALSE, start = 1,
  end = NA, idx = NA) {

  if ( is.na(end) ) end <- x$nmc
  if ( end <= start ) stop("End must be greater than start")
  if (x$n.chain < nchain) stop("nchain is too large")

  if (is.na(idx)) idx <- sample(1:x$n.chains, nchain)
  lp <- x$summed_log_prior[start:end,idx] + x$log_likelihoods[start:end,idx]
  d <- coda::mcmc.list(lapply(data.frame(lp),
    function(xx){coda::mcmc(as.matrix(xx))}))

  DT <- ggmcmc::ggs(d)

  DT$Chain <- factor(DT$Chain)
  levels(DT$Chain)

  DT$Parameter <- "lp"
  f2 <- ggplot(DT, aes_string(x = "Iteration", y = "value", color = "Chain")) +
    geom_line() +
    ylab("Log-posterior likelihood") +
    theme_minimal() +
    theme(legend.position = "none")
  print(f2)
  return(invisible(DT))
}

##' Profile a model object
##'
##' This function produces data for profiling model likelihoods based on a data
##' model instance (ie \code{fitted}).
##'
##' The argument, \code{pname} indicates which model parameter to profile.
##' For example, if we want to profile the boundary separation \emph{a} in a
##' DDM, we extract it from a \code{p.vector} and calculates its marginal
##' likelihoods based on the data model instance on a grid of \code{n.point}.
##'
##' Briefly, the function sets a range for the points on the x axis, which
##' is the values for the profiled parameter (e.g., \emph{a}). Then it
##' initiates a log-likelihodd vector with length of \code{n.point} and 0
##' everywhere.  Next, it keeps the other parameters in the model fixed and
##' changes only the target parameter values to ps[i]. After that, it
##' calculates their sum log-likelihoods. Finally, it stores the
##' log-likelihoods in the \code{i} position of the ll vector.
##'
##' @param fitted a data model instance
##' @param pname indicate which parameter in theta to plot. For example, in a
##' LBA model with a \code{pVec <- c(A="1",B="1",mean_v="M",sd_v="1",t0="1",
##' st0="1")}, one can assign \code{pname <- "A"} to ask \code{profile} to
##' profile \emph{A} parameter
##' @param minp lower bound for pname parameter
##' @param maxp upper bound for pname parameter
##' @param p.vector a parameter vector. Use Lognromal LBA model as an example,
##' \code{pVec <- c(A = .25, B = .35, meanlog_v.true = 1, meanlog_v.false = .25,
##' sdlog_v.true = .5, t0 = .2)}
##' @param npoint grid for p.name parameter
##' @param digits print out how many digits
##' @param ylim y range
##' @param nthread number of thread in a GPU block. This argument works for
##' GPU-based PDA probability density functions)
##' @param ... additional optional arguments.
##' @importFrom graphics plot
##' @examples
##' model <- BuildModel(
##'         p.map     = list(a = "1",v = "1",z = "1", d = "1", sz = "1",
##'                     sv = "1", t0 = "1", st0 = "1"),
##'         constants = c(st0 = 0, d = 0),
##'         match.map = list(M = list(s1 = "r1", s2 = "r2")),
##'         factors   = list(S = c("s1", "s2")),
##'         responses = c("r1", "r2"),
##'         type      = "rd")
##'
##' p.prior <- BuildPrior(
##'   dists = rep("tnorm", 6),
##'   p1=c(a=2,   v=2.5, z=0.5, sz=0.3, sv=1,  t0=0.3),
##'   p2=c(a=0.5, v=.5,  z=0.1, sz=0.1, sv=.3, t0=0.05),
##'   lower=c(0,-5, 0, 0, 0, 0),
##'   upper=c(5, 7, 2, 2, 2, 2))
##'
##' p.vector <- c(a=1,v=1, z=0.5, sz=0.25, sv=0.2,t0=.15)
##' dat <- simulate(model, 1e2, ps = p.vector)
##' dmi <- BuildDMI(dat, model)
##'
##' ## ------------------------------40
##' par(mfrow=c(2,3));
##' profile(dmi, "a",  .1,   2, p.vector)
##' profile(dmi, "v",  .1,   2, p.vector)
##' profile(dmi, "z",  .2,  .8, p.vector)
##' profile(dmi, "sz", .1,  .9, p.vector)
##' profile(dmi, "sv", .1,   2, p.vector)
##' profile(dmi, "t0", .01, .5, p.vector)
##' par(mfrow=c(1,1));
##' @export
profile.model <- function(fitted, pname, minp, maxp, p.vector,
  npoint = 100, digits = 2, ylim = NA, nthread = 32, ...) {

  if (!(pname %in% names(p.vector))) stop("parameter not in p.vector")

  RT        <- fitted$RT
  model     <- attr(fitted, "model")
  ise       <- attr(fitted, "cell.empty")
  allpar    <- attr(model, "all.par")
  parnames  <- attr(model, "par.names")
  type      <- attr(model, "type")
  n1idx     <- attr(model, "n1.order")
  mc        <- attr(model, "match.cell")
  isr1      <- check_rd(type, model)
  cellidx   <- cellIdx2Mat(fitted)
  pnames    <- names(p.vector)
  ps        <- seq(minp, maxp, length.out = npoint)
  nsim      <- attr(fitted, "n.pda")
  bw        <- attr(fitted, "bw")
  gpuid     <- attr(fitted, "gpuid")
  debug     <- attr(fitted, "debug")

  if (type == "norm") {
    posdrift <- attr(model, "posdrift")

    ll <- profile_norm(p.vector, pnames, allpar, parnames, model, type,
      dimnames(model)[[1]], dimnames(model)[[2]], dimnames(model)[[3]],
      n1idx, ise, cellidx, RT, mc, isr1, pname, ps, posdrift)

  } else if (type == "norm_pda") {

    ll <- profile_norm_pda(p.vector, pnames, allpar, parnames, model, type,
      dimnames(model)[[1]], dimnames(model)[[2]], dimnames(model)[[3]], n1idx,
      ise, cellidx, RT, mc, isr1, pname, ps, nsim, bw)

  } else if (type == "norm_pda_gpu") {

    ll <- profile_norm_gpu(p.vector, pnames, allpar, parnames, model, type,
      dimnames(model)[[1]], dimnames(model)[[2]], dimnames(model)[[3]], n1idx,
      ise, cellidx, RT, mc, isr1, pname, ps, nsim, bw)

  } else if (type == "plba0_gpu") {
    ll <- profile_plba0_gpu(p.vector, pnames, allpar, parnames, model, type,
      dimnames(model)[[1]], dimnames(model)[[2]], dimnames(model)[[3]], n1idx,
      ise, cellidx, RT, mc, isr1, pname, ps, nsim, bw, gpuid, nthread,
      debug)

  } else if (type == "plba1") {

    ll <- profile_plba1(p.vector, pnames, allpar, parnames, model, type,
      dimnames(model)[[1]], dimnames(model)[[2]], dimnames(model)[[3]], n1idx,
      ise, cellidx, RT, mc, isr1, pname, ps, nsim, bw)

  } else if (type == "plba1_gpu") {

    ll <- profile_plba1_gpu(p.vector, pnames, allpar, parnames, model, type,
      dimnames(model)[[1]], dimnames(model)[[2]], dimnames(model)[[3]], n1idx,
      ise, cellidx, RT, mc, isr1, pname, ps, nsim, bw, gpuid, nthread,
      debug)

  } else if (type == "rd") {

    ll <- profile_rd(p.vector, pnames, allpar, parnames, model, type,
      dimnames(model)[[1]], dimnames(model)[[2]], dimnames(model)[[3]],
      n1idx, ise, cellidx, RT, mc, isr1, pname, ps)

  } else if (type == "cnorm") {

    ll <- profile_cnorm_pda(p.vector, pnames, allpar, parnames, model, type,
      dimnames(model)[[1]], dimnames(model)[[2]], dimnames(model)[[3]], n1idx,
      ise, cellidx, RT, mc, isr1, pname, ps, nsim, bw)

  } else {
    stop("Type not defined")
  }

  names(ll) <- round(ps, digits)
  plot(ps, ll, type = "l", xlab = pname, ylab = "log-likelihood")
  ll[ll==max(ll)]
}


##' Plot Prior Distributions
##'
##' \code{plot_prior} plots the ith member of the list created by
##' \code{BuildPrior}.  If \code{trans = TRUE}, the function will plot on natural
##' (logarithmic) scale using transform specified in
##' \code{attr(p.prior[[i]], "trans")}. \code{plot.prior} plots all parameters
##' at once.
##'
##' NOTE: the tran function is not thoroughtly checked. Use it with your own risk.
##' \code{plot_prior} checks if any of the elements in \code{trans} is \code{NA}.
##' It then checks if \code{natural} is set \code{TRUE}
##' (by default it's \code{TRUE}). If \code{natural} is set also
##' \code{TRUE} (i.e., the user wants to do transformation), it then checks what
##' has been set in the "untrans" attribute for the i'th parameter.  Otherwise,
##' its default is set all parameters as \code{identity}.
##'
##' @param x prior distributions
##' @param i an integer or a character string indicating to plot which parameter
##' @param prior a list of list storing prior setting.
##' @param xlim set the range of on x axis. This is usually the range for each
##' parameter.
##' @param natural default TRUE
##' @param npoint default to plot 100
##' @param trans default NA. trans can be a scalar or vector.
##' @param save whether to save the data out
##' @param ... other plotting arguments passing through dot dot dot.
##' @import ggplot2
##' @export
##' @examples
##' p.prior <- BuildPrior(
##'            dists = rep("tnorm", 7),
##'            p1    = c(a = 2,   v.f1 = 4,  v.f2 = 3,  z = 0.5, sv = 1,  sz = 0.3, t0 = 0.3),
##'            p2    = c(a = 0.5, v.f1 = .5, v.f2 = .5, z = 0.1, sv = .3, sz = 0.1, t0 = 0.05),
##'            lower = c(0,-5, -5, 0, 0, 0, 0),
##'            upper = c(5, 7,  7, 1, 2, 1, 1))
##'
##' plot_prior("a", p.prior)
##' plot_prior(2, p.prior)
##'
##' print(p.prior)
##' plot(p.prior)
##'
##' ## require(ggplot2)
##' ## p2 <- ggplot(d, aes(x = xpos, y = ypos)) + geom_line() +
##' ##       facet_wrap(~gpvar, scales="free") + theme_bw(base_size =14)
plot_prior <- function(i, prior, xlim = NA, natural = TRUE, npoint = 100,
                       trans = NA, save = FALSE, ... ) {

  if (any(is.na(trans))) {   ### Check 1 ###
    trans <- ifelse(natural, attr(prior[[i]], "untrans"), "identity")
  }

  if (is.numeric(i)) i <- names(prior)[i]   ### Check 2 ###
  ## Stop the function running, if the parameters are out-of-range
  if (!(i %in% names(prior))) stop("Parameter not in prior")
  ## Finish checks. Calculate the data frame
  p <- prior[[i]]    ## Save all parameter setting to another object called p
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
                      Parameter = rep(names(prior[i]), length(x)))
    return(dat)
  } else {
    xpos <- do.call(trans, list(x = x))
    ypos <- do.call(paste0("d", disttype), p)
    dat <- data.frame(xpos  = xpos, ypos  = ypos,
                      Parameter = rep(names(prior[i]), length(x)))

    p0 <- ggplot(dat, aes_string(x="xpos", y="ypos")) + geom_line() +
      xlab("")+ ylab("Density") + facet_grid(. ~Parameter) +
      theme_bw(base_size = 16)
    print(p0)
    invisible(p0)
  }

}

##' @rdname plot_prior
##' @export
plot.prior <- function(x, save = FALSE, ps = NULL, ...) {
  if( is.null(names(x))) stop("Prior object must has name attributes.")

  pD <- NULL
  for(j in names(x)) {
    invisible(pD <- rbind(pD, plot_prior(i = j, prior = x, save = TRUE)))
  }

  if (save) {
    return(pD)
  } else if (!is.null(ps) & is.vector(ps)) {
    pveclines <- data.frame(Parameter = names(x), true = ps)
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
    wide <- data.frame(ps)
    wide$s <- factor(1:nrow(ps))
    pveclines <- melt(wide, id.vars = "s", variable.name = "Parameter",
                 value.name = "true")
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
}



