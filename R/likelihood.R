##' Calculate log likelihoods
##'
##' These function calculate log likelihoods. \code{likelihood_rd} implements
##' the equations in Voss, Rothermund, and Voss (2004). These equations
##' calculate diffusion decision model (Ratcliff & Mckoon, 2008). Specifically,
##' this function implements Voss, Rothermund, and Voss's (2004) equations A1
##' to A4 (page 1217) in C++.
##'
##' @param x a parameter vector
##' @param data data model instance
##' @param min_lik minimal likelihood.
##' @return a vector
##' @references Voss, A., Rothermund, K., & Voss, J. (2004).  Interpreting the
##' parameters of the diffusion model: An empirical validation.
##' \emph{Memory & Cognition}, \bold{32(7)}, 1206-1220. \cr\cr
##' Ratcliff, R. (1978). A theory of memory retrival. \emph{Psychological
##' Review}, \bold{85}, 238-255.
##'
##' @examples
##' model <- BuildModel(
##' p.map     = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1",
##'             st0 = "1"),
##' match.map = list(M = list(s1 = 1, s2 = 2)),
##' factors   = list(S = c("s1", "s2")),
##' constants = c(st0 = 0, sd_v = 1),
##' responses = c("r1", "r2"),
##' type      = "norm")
##'
##' p.vector <- c(A = .25, B = .35,  t0 = .2, mean_v.true = 1, mean_v.false = .25)
##' dat <- simulate(model, 1e3,  ps = p.vector)
##' dmi <- BuildDMI(dat, model)
##' den <- likelihood_norm(p.vector, dmi)
##' ## den <- likelihood_norm_pda(p.vector, dmi) ## This takes more than 1 s, so commented out
##'
##' model <- BuildModel(
##' p.map     = list(a = "1", v = "1", z = "1", d = "1", t0 = "1", sv = "1",
##'             sz = "1", st0 = "1"),
##' constants = c(st0 = 0, d = 0),
##' match.map = list(M = list(s1 = "r1", s2 = "r2")),
##' factors   = list(S = c("s1", "s2")),
##' responses = c("r1", "r2"),
##' type      = "rd")
##'
##' p.vector <- c(a = 1, v = 1, z = 0.5, sz = 0.25, sv = 0.2, t0 = .15)
##' dat <- simulate(model, 1e2, ps = p.vector)
##' dmi <- BuildDMI(dat, model)
##' den <- likelihood_rd(p.vector, dmi)
##'
##' @export
likelihood_norm <- function(x, data, min_lik = 1e-10) {
  model    <- attr(data,  "model")
  ise      <- attr(data,  "cell.empty")
  allpar   <- attr(model, "all.par")
  parnames <- attr(model, "par.names")
  type     <- attr(model, "type")
  n1idx    <- attr(model, "n1.order")
  mc       <- attr(model, "match.cell")
  isr1     <- check_rd(type, model)
  cellidx  <- cellIdx2Mat(data)
  pnames   <- names(x)
  posd     <- attr(model, "posdrift")

  out <- density_norm(x, pnames, allpar, parnames, model, type,
    dimnames(model)[[1]], dimnames(model)[[2]], dimnames(model)[[3]],
    n1idx, ise, cellidx, data$RT, mc, isr1, posd)
  pmax(out[,1], min_lik)
}

##' @rdname likelihood_norm
##' @export
likelihood_norm_pda <- function(x, data, min_lik = 1e-10) {
  model    <- attr(data,  "model")
  ise      <- attr(data,  "cell.empty")
  allpar   <- attr(model, "all.par")
  parnames <- attr(model, "par.names")
  type     <- attr(model, "type")
  n1idx    <- attr(model, "n1.order")
  mc       <- attr(model, "match.cell")
  isr1     <- check_rd(type, model)
  cellidx  <- cellIdx2Mat(data)
  pnames   <- names(x)
  nsim     <- attr(data, "n.pda")
  bw       <- attr(data, "bw")

  out <- density_norm_pda(x, pnames, allpar, parnames, model, type,
              dimnames(model)[[1]], dimnames(model)[[2]], dimnames(model)[[3]],
              n1idx, ise, cellidx, data$RT, mc, isr1, nsim, bw)
  pmax(out[,1], min_lik)
}

##' @rdname likelihood_norm
##' @export
likelihood_rd <- function(x, data, min_lik = 1e-10) {
  model    <- attr(data,  "model")
  ise      <- attr(data,  "cell.empty")
  allpar   <- attr(model, "all.par")
  parnames <- attr(model, "par.names")
  type     <- attr(model, "type")
  n1idx    <- attr(model, "n1.order")
  mc       <- attr(model, "match.cell")
  isr1     <- check_rd(type, model)
  cellidx  <- cellIdx2Mat(data)
  pnames   <- names(x)
  out <- density_rd(x, pnames, allpar, parnames, model, type,
                    dimnames(model)[[1]], dimnames(model)[[2]],
                    dimnames(model)[[3]],
                    n1idx, ise, cellidx, data$RT, mc, isr1)
  pmax(out[,1], min_lik)
}

##' @rdname likelihood_norm
##' @export
likelihood_cnorm <- function(x, data, min_lik = 1e-10) {
  model    <- attr(data,  "model")
  ise      <- attr(data,  "cell.empty")
  allpar   <- attr(model, "all.par")
  parnames <- attr(model, "par.names")
  type     <- attr(model, "type")
  n1idx    <- attr(model, "n1.order")
  mc       <- attr(model, "match.cell")
  isr1     <- check_rd(type, model)
  cellidx  <- cellIdx2Mat(data)
  pnames   <- names(x)
  dim1 <- dimnames(model)[[1]]
  dim2 <- dimnames(model)[[2]]
  dim3 <- dimnames(model)[[3]]
  nsim     <- attr(data, "n.pda")
  bw       <- attr(data, "bw")

  out <- density_cnorm_pda(x, pnames, allpar, parnames, model, type,
    dim1, dim2, dim3, n1idx, ise, cellidx, data$RT, mc, isr1, nsim, bw)
  pmax(out[,1], min_lik)
}

##' Substract nondecision times from response times
##'
##' This function minus t0's from RTs.
##'
##' @param x a response time vector
##' @param t0 a numeric scalar
##' @return a vector
##'
##' @examples
##' rt <- rlnorm(10) + .2
##' dt <- removet0(rt, .2)
##' all.equal(dt, rt - .2)
##' @export
removet0 <- function(x, t0) {
  if (!is.vector(x)) stop("RT must be a vector")
  if (!is.numeric(t0)) stop("t0 must be a number")
  out <- remove_t0(x, t0)
  return(out[,1])
}

n1PDFfixedt0_lnr <- function (dt, meanlog, sdlog) {
  ## incompleted
  if (is.matrix(meanlog) & is.matrix(sdlog)) {
    out <- n1PDFfixedt0_lnr2(dt, meanlog, sdlog)
  } else if (is.vector(meanlog) & is.vector(sdlog)) {
    out <- n1PDFfixedt0_lnr1(dt, meanlog, sdlog)
  } else {
    stop("meanlog/sdlog type not supported.")
  }
  return(out[,1]);
}


