### Generic  ---------------------------------------------

##' @export
plot.model <- function(x, y = NULL, hyper = FALSE, start = 1,
  end = NA, pll = TRUE, save = FALSE, den = FALSE, subchain = FALSE,
  nsubchain = 3, chains = NA, ...)
{
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


### One Subject  ---------------------------------------------
##' @import ggplot2
##' @importFrom coda mcmc mcmc.list
preplot_one <- function(x, start, end, pll) {

  if ( is.na(end) ) end <- x$nmc
  if ( end <= start ) stop("End must be greater than start")
  nchain <- x$n.chains

  if (pll) {
    ## nchain x nmc
    lp <- x$summed_log_prior[,start:end] + x$log_likelihoods[,start:end]

    rownames(lp) <- 1:nchain
    step1 <- lapply(data.frame(t(lp)), function(xx){
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
    DT$Parameter <- "lp"

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
  DT <- preplot_many(x, start, end, pll)

  if (subchain) {
    if (missing(nsubchain)) stop("Please supply nsubchain")
    if (any(is.na(chains))) chains <- base::sample(x[[1]]$n.chains, nsubchain)
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
  # x <- fit1
  # start <- 1
  # pll <- FALSE
  xx <- attr(x, "hyper")
  if ( is.na(end) ) end <- xx$nmc
  if ( end <= start ) stop("End must be greater than start")
  nchain <- xx$n.chains
  # thin   <- xx$thin

  if (pll) {
    lp <- xx$h_summed_log_prior[,start:end] + xx$h_log_likelihoods[,start:end]
    rownames(lp) <- 1:nchain
    step1 <- lapply(data.frame(t(lp)), function(xxx){
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


### Prior  ---------------------------------------------
##' Plot prior distributions
##'
##' \code{plot_prior} plots one member in a prior object. \code{plot.prior}
##' plots all members in a prior object.
##'
##' @param i an integer or a character string indicating which parameter to plot
##' @param x a prior object
##' @param prior a prior object
##' @param xlim set the range of on x axis. This is usually the range for each
##' parameter.
##' @param natural default TRUE.
##' @param npoint default to plot 100
##' @param trans default NA. trans can be a scalar or vector.
##' @param save whether to save the data out
##' @param ps true parameter vectors or matrix in the case of many observation
##' units
##' @param ... other plotting arguments passing through dot dot dot.
##' @import ggplot2
##' @export
##' @examples
##' p.prior <- BuildPrior(
##'            dists = rep("tnorm", 7),
##'            p1    = c(a = 2,   v.f1 = 4,  v.f2 = 3,  z = 0.5, sv = 1,
##'                      sz = 0.3, t0 = 0.3),
##'            p2    = c(a = 0.5, v.f1 = .5, v.f2 = .5, z = 0.1, sv = .3,
##'                      sz = 0.1, t0 = 0.05),
##'            lower = c(0,-5, -5, 0, 0, 0, 0),
##'            upper = c(5, 7,  7, 1, 2, 1, 1))
##'
##' plot_prior("a", p.prior)
##' plot_prior(2, p.prior)
##' plot(p.prior)
plot_prior <- function(i, prior, xlim = NA, natural = TRUE, npoint = 100,
                       trans = NA, save = FALSE, ... ) {

  # plot_prior("a", p.prior)
  # plot_prior(2, p.prior)
  # i <- 2
  # natural <- F
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

##' @import ggplot2
##' @importFrom data.table data.table melt.data.table
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
}



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


