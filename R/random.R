##' Piecewise LBA model
##'
##' Density and random generation of the PLBA Model Type 0, 1, and 2.
##'
##' @param n number of observations.
##' @param x vector of quantiles.
##' @param A upper bound of start point. It can be an integer or a 2-element
##' vector.
##' @param b response threshold. It can be an integer or a 2-element vector.
##' @param mean_v stage 1 mean drift rate. It must be a 2-element vector
##' @param mean_w stage 2 mean drift rate. It must be a 2-element vector
##' @param sd_v common standard deviation of the piece 1 drift rates. If
##' sd_w is not present, this will also be used as the piece 2 drift rate
##' standard deviation, which cannot be negative.
##' @param sd_w standard deviation of the piece 2 drift rates
##' @param rD drift rate delay (in second)
##' @param swt switch time (in second)
##' @param t0 nondecision time (in second)
##' @param h bandwidth for the kernel function
##' @param ncore number of CPU cores for running Open MP.
##' @param debug internal debug switch
##' @param B first stage traveling distance
##' @param C second stage traveling distance
##' @param tD threshold delay time
##' @param pVec PLBA parameter vector
##' @return a [RT R] matrix (C++) or a data frame (R)
##' @references Holmes, R. W., Trueblood, J. & Heathcote, A. (2016). A new
##' framework for modeling decisions about changing information: The Piecewise
##' Linear Ballistic Accumulator model. \emph{Cognitive Psychology}, 85, 1--29,
##' doi: http://dx.doi.org/10.1016/j.cogpsych.2015.11.002.Approximate
##' @examples
##' #########################################################################80
##' ## rplba1
##' #########################################################################80
##' \dontrun{
##' n <- 2^20; n
##' A <- 1.5
##' b <- 2.7
##' mean_v <- c(3.3, 2.2)
##' mean_w <- c(1.5, 3.7)
##' sd_v <- c(1, 1)
##' rD    <- .3
##' swt   <- .5
##' t0    <- .08
##' ncore <- 12
##' dat1 <- rplba1R(n, A, b, t0, mean_v, mean_w, sd_v, rD, swt)
##' dat2 <- rplba1(n, A, b, t0, mean_v, mean_w, sd_v, rD, swt, ncore)
##' dat3 <- ppda::rplba1(n, A, b, t0, mean_v, mean_w, sd_v, rD, swt)
##'
##' dat1r1 <- dat1[dat1[, 2] == 1, 1]
##' dat1r2 <- dat1[dat1[, 2] == 2, 1]
##' dat2r1 <- dat2[dat2[, 2] == 1, 1]
##' dat2r2 <- dat2[dat2[, 2] == 2, 1]
##' dat3r1 <- dat3[dat3[, 2] == 1, 1]
##' dat3r2 <- dat3[dat3[, 2] == 2, 1]
##'
##' xlim <- c(0, 3)
##' ## Check if two methods produce SPDF overlaping with each other
##' par(mfrow = c(4, 2), mar = c(4, 5.3, 0.82, 1))
##' hist(dat1r1, breaks = "fd", freq = FALSE, main = "Choice1 R", xlim = xlim)
##' hist(dat1r2, breaks = "fd", freq = FALSE, main = "Choice2 R", xlim = xlim)
##' hist(dat2r1, breaks = "fd", freq = FALSE, main = "Choice1 C++", xlim = xlim)
##' hist(dat2r2, breaks = "fd", freq = FALSE, main = "Choice2 C++", xlim = xlim)
##' hist(dat3r1, breaks = "fd", freq = FALSE, main = "Choice1 GPU", xlim = xlim)
##' hist(dat3r2, breaks = "fd", freq = FALSE, main = "Choice2 GPU", xlim = xlim)
##'
##' par(mfrow = c(1, 2))
##' hist(dat1r1, breaks = "fd", freq = FALSE, main = "Choice1 R, C++, & GPU",
##'   xlim = xlim, ylim = c(0, 3))
##' hist(dat2r1, breaks = "fd", freq = FALSE, add = TRUE, col = "lightblue")
##' hist(dat3r1, breaks = "fd", freq = FALSE, add = TRUE, col = "lightgreen")
##'
##' hist(dat1r2, breaks = "fd", freq = FALSE, main = "Choice2 R, C++ & GPU",
##'   xlim = xlim, ylim = c(0, 3))
##' hist(dat2r2, breaks = "fd", freq = FALSE, add = TRUE, col = "lightblue")
##' hist(dat3r2, breaks = "fd", freq = FALSE, add = TRUE, col = "lightgreen")
##' }
##'
##' #############20
##' ## rplba2    ##
##' #############20
##' \dontrun{
##' n <- 2^15
##' ncore <- 4
##' A <- c(1.5, 1.5)
##' b <- c(2.7, 2.7)
##' mean_v <- c(3.3, 2.2)
##' mean_w <- c(1.5, 3.7)
##' sd_v <- c(1, 1)
##' sd_w <- c(1, 1)
##' rD <- .3
##' swt <- .5
##' t0 <- .08
##' dat1 <- rplba2R(n, A, b, t0, mean_v, mean_w, sd_v, sd_w, rD, swt)
##' dat2 <- rplba2(n, A, b, t0, mean_v, mean_w, sd_v, sd_w, rD, swt, ncore)
##' dat3 <- rplba2(n, A, b, t0, mean_v, mean_w, sd_v, sd_w, rD, swt)
##' dat4 <- rplba2_test(n, A, b, t0, mean_v, mean_w, sd_v, sd_w, rD, swt)
##'
##' dat1r1 <- dat1[dat1[, 2] == 1, 1]
##' dat1r2 <- dat1[dat1[, 2] == 2, 1]
##' dat2r1 <- dat2[dat2[, 2] == 1, 1]
##' dat2r2 <- dat2[dat2[, 2] == 2, 1]
##' dat3r1 <- dat3[dat3[, 2] == 1, 1]
##' dat3r2 <- dat3[dat3[, 2] == 2, 1]
##' dat4r1 <- dat4[dat4[, 2] == 1, 1]
##' dat4r2 <- dat4[dat4[, 2] == 2, 1]
##'
##' wesanderson::wes_palette("Royal1")
##' palettes  <- wesanderson::wes_palettes$GrandBudapest
##' palettes2 <- wesanderson::wes_palettes$GrandBudapest2
##' xlim <- c(0, 3)
##' ## Check if two methods produce SPDF overlaping with each other
##' par(mfrow = c(4, 2), mar = c(4, 5.3, 0.82, 1))
##' hist(dat1r1, breaks = "fd", freq = FALSE, main = "Choice1 R", xlim = xlim)
##' hist(dat1r2, breaks = "fd", freq = FALSE, main = "Choice2 R", xlim = xlim)
##' hist(dat2r1, breaks = "fd", freq = FALSE, main = "Choice1 C++", xlim = xlim)
##' hist(dat2r2, breaks = "fd", freq = FALSE, main = "Choice2 C++", xlim = xlim)
##' hist(dat3r1, breaks = "fd", freq = FALSE, main = "Choice1 GPU", xlim = xlim)
##' hist(dat3r2, breaks = "fd", freq = FALSE, main = "Choice2 GPU", xlim = xlim)
##' hist(dat4r1, breaks = "fd", freq = FALSE, main = "Choice1 test", xlim = xlim)
##' hist(dat4r2, breaks = "fd", freq = FALSE, main = "Choice2 test", xlim = xlim)
##'
##' par(mfrow = c(1, 2))
##' hist(dat1r1, breaks = "fd", freq = FALSE, main = "Choice1 R, C++, & GPU",
##'   xlim = xlim, ylim = c(0, 3))
##' hist(dat2r1, breaks = "fd", freq = FALSE, add = TRUE, col = palettes[1])
##' hist(dat3r1, breaks = "fd", freq = FALSE, add = TRUE, col = palettes[2])
##' hist(dat4r1, breaks = "fd", freq = FALSE, add = TRUE, col = palettes[4])
##'
##' hist(dat1r2, breaks = "fd", freq = FALSE, main = "Choice2 R, C++ & GPU",
##'   xlim = xlim, ylim = c(0, 3))
##' hist(dat2r2, breaks = "fd", freq = FALSE, add = TRUE, col = palettes2[1])
##' hist(dat3r2, breaks = "fd", freq = FALSE, add = TRUE, col = palettes2[2])
##' hist(dat4r2, breaks = "fd", freq = FALSE, add = TRUE, col = palettes2[3])
##'
##'
##' }
##' @importFrom stats runif
##' @export
rplba1R <- function(n, A, b, t0, mean_v, mean_w, sd_v, rD, swt) {
  if (length(mean_v) != 2) stop("Current version fit only two accumulators.")
  if (length(sd_v) != 2) stop("The length of sd_v and mean_v must match.")
  eswt <- swt + rD
  v1 <- rtnorm(n, mean_v[1], sd_v[1], 0, Inf)[,1] ## Stage 1 LBA
  v2 <- rtnorm(n, mean_v[2], sd_v[2], 0, Inf)[,1]
  sp <- matrix(runif(2*n, 0, A), 2)
  dt1 <- rbind((b - sp[1,])/v1, (b - sp[2,])/v2) ## Race

  ## dt[dt<0] <- Inf
  choice    <- apply(dt1, 2, which.min)
  chosen_dt <- dt1[cbind(choice, 1:n)]  ## convert to vector choose (row, col)

  done <- (chosen_dt <= eswt)   ## Which are finished?
  n2 <- sum(!done)

  ## Distance left to travel for those not finished
  B1 <- b - (sp[1, !done] + eswt*v1[!done])
  B2 <- b - (sp[2, !done] + eswt*v2[!done])

  w1 <- rtnorm(n2, mean_w[1], sd_v[1], 0, Inf)[,1]   ## Stage 2 LBA
  w2 <- rtnorm(n2, mean_w[2], sd_v[2], 0, Inf)[,1]
  dt2 <- rbind(B1/w1, B2/w2)   ## Race

  choice[!done] <- apply(dt2, 2, which.min)
  chosen_dt[!done] <- eswt + dt2[cbind(choice[!done], 1:n2)]

  ## The last object automatically return
  data.frame(cbind(RT=t0 + chosen_dt, R = choice))
}

##' @importFrom stats runif
##' @rdname rplba1R
##' @export
rplba2R <- function(n, A, b, t0, mean_v, mean_w, sd_v, sd_w, rD, swt)
{
  # Calculate effective switch time
  eswt <- swt + rD

  # Stage 1 LBA Race
  v <- t(cbind(rtnorm(n, mean_v[1], sd_v[1], 0, Inf), rtnorm(n, mean_v[2], sd_v[2], 0, Inf)))
  sp <- matrix(runif(2*n, 0, A), 2)
  dt <- (b- sp) / v
  # dt[dt<0] <- Inf
  choice <- apply(dt, 2, which.min)
  chosen_dt <- dt[cbind(choice, 1:n)]

  # Which are finished?
  done <- (chosen_dt <= eswt)
  n2   <- sum(!done)

  # Distance left to travel for those not finished
  B <- b - (sp[, !done] + eswt * v[, !done])

  # Stage 2 LBA Race
  w <- t(cbind(rtnorm(n2, mean_w[1], sd_w[1], 0, Inf), rtnorm(n2, mean_w[2], sd_w[2], 0, Inf)))
  dt <- B / w
  choice[!done] <- apply(dt, 2, which.min)
  chosen_dt[!done] <- eswt + dt[ cbind(choice[!done], 1:n2) ]

  # save results
  ## The last object automatically return
  data.frame(cbind(RT = t0 + chosen_dt, R = choice))
}

##' @importFrom stats runif
##' @rdname rplba1R
##' @export
rplba3R <- function(n=10, pVec=c(A1=1.5, A2=1.5, B1=1.2, B2=1.2, C1=.3, C2=.3,
  v1=3.32, v2=2.24, w1=1.51, w2=3.69, sv1=1, sv2=1,
  sw1=1, sw2=1, rD=0.3, tD=.3, swt=0.5, t0=0.08)) {

  # n <- 10
  # pVec=c(A1=1.5, A2=1.5, B1=1.2, B2=1.2, C1=.3, C2=.3,
  #   v1=3.32, v2=2.24, w1=1.51, w2=3.69, sv1=1, sv2=1,
  #   sw1=1, sw2=1, rD=0.3, tD=.3, swt=0.5, t0=0.08)

  # Stage 1 LBA
  v1 <- rtnorm(n, pVec["v1"], pVec["sv1"], 0, Inf)[,1]
  v2 <- rtnorm(n, pVec["v2"], pVec["sv2"], 0, Inf)[,1]
  sp <- matrix(runif(2*n,0,pVec[c("A1","A2")]),nrow=2)

  # Calcualte thresholds
  b1 <- sum(pVec[c("A1","B1")])
  b2 <- sum(pVec[c("A2","B2")])
  c1 <- b1 + pVec[c("C1")]
  c2 <- b2 + pVec[c("C2")]

  # Race
  dt <- rbind((c(b1,b2)-sp[1,])/v1,(c(b1,b2)-sp[2,])/v2)
  # dt[dt<0] <- Inf
  choice <- apply(dt,2,which.min)
  rt <- dt[cbind(choice,1:n)]

  # Calculate effective switch times
  swt_b <- pVec["swt"] + pVec["tD"]
  swt_r <- pVec["swt"] + pVec["rD"]

  # Which switch is first
  swt <- pmin(swt_b,swt_r)
  if (swt_b==swt_r) {
    change <- "both"
  } else if (swt_r < swt_b) {
    change <- "rate"
  } else {
    change <- "threshold"
  }

  # Which are finished?
  done <- rt <= swt
  n2 <- sum(!done)

  # Stage 2 LBA

  # Distance left to travel for those not finished
  # threshold - distance already travelled
  if ( change=="rate" ) {
    B1 <- b1 - (sp[1,!done] + swt*v1[!done])
    B2 <- b2 - (sp[2,!done] + swt*v2[!done])
  } else {
    B1 <- c1 - (sp[1,!done] + swt*v1[!done])
    B2 <- c2 - (sp[2,!done] + swt*v2[!done])
  }


  # Change rates?
  if ( change=="threshold" ) {
    w1 <- v1[!done]; w2 <- v2[!done]
  } else {
    w1 <- rtnorm(n2, pVec["w1"], pVec["sw1"],0, Inf)[,1]
    w2 <- rtnorm(n2, pVec["w2"], pVec["sw2"],0, Inf)[,1]
  }

  # Race
  dt <- rbind(B1/w1,B2/w2)
  # dt[dt<0] <- Inf
  choice[!done] <- apply(dt,2,which.min)
  rt[!done] <- swt+dt[cbind(choice[!done],1:n2)]

  if ( change != "both" ) { # Stage 3 LBA

    if ( change=="threshold" ) swt1 <- swt_r else swt1 <- swt_b
    t2 <- swt1-swt

    # Which are finished?
    done1 <- rt[!done] < swt1
    n2 <- sum(!done1)

    if ( !all(done1) ) {

      # Distance left to travel for those not finished
      # Distance left at end of stage 1 - further travel
      B1 <- B1[!done1] - t2*w1[!done1]
      B2 <- B2[!done1] - t2*w2[!done1]

      if ( change=="threshold" ) {
        w1 <- rtnorm(n2,pVec["w1"],pVec["sw1"],0, Inf)[,1]
        w2 <- rtnorm(n2,pVec["w2"],pVec["sw2"],0, Inf)[,1]
      }  else {
        w1 <- w1[!done1];
        w2 <- w2[!done1]
        B1 <- B1 + pVec["C1"]
        B2 <- B2 + pVec["C2"]
      }

      # Race
      dt <- rbind(B1/w1,B2/w2)
      # dt[dt<0] <- Inf
      choice[!done][!done1] <- apply(dt,2,which.min)
      rt[!done][!done1] <- swt1+dt[cbind(choice[!done][!done1],1:n2)]
    }

  }

  # save results
  data.frame(R=choice, RT=pVec["t0"] + rt)
}

rlnrR <- function (n, meanlog, sdlog, t0, st0 = 0) {
  n_acc <- ifelse(is.null(dim(meanlog)), length(meanlog), dim(meanlog)[1])
  dt    <- matrix(rlnorm(n = n*n_acc, meanlog = meanlog, sdlog = sdlog),
               nrow = n_acc) + t0
  winner <- apply(dt, 2, which.min)
  if (st0[1] == 0) {
    out <- data.frame(RT = dt[cbind(winner, 1:n)], R = winner)
  } else {
    out <- data.frame(RT = dt[cbind(winner, 1:n)] + runif(n, 0, st0[1]),
      R = winner)
  }
  return(out)
}

##' Generate random number from a correlated accumulator model
##'
##' @param n number of observatoins
##' @param A start point variability
##' @param b decision threshold
##' @param t0 non-decision time
##' @param mean_v mean drift rate vector
##' @param sd_v standard deviation of the drift rates
##' @param st0 non-decision time variability
##' @param corr_v correlations among accumulators
##' @param return_ttf a Boolean switch to return time to finish matrix
##' @importFrom tmvtnorm rtmvnorm
##' @export
rca <- function (n, A, b, t0, mean_v, sd_v, st0 = 0, corr_v = 0,
  return_ttf = FALSE) {

  ## Turn off check to prevent Bayesian stop
  ## if (any(b < A)) stop("b cannot be smaller than A!")
  nv    <- ifelse(is.null(dim(mean_v)), length(mean_v), nrow(mean_v))
  Sigma <- make_sigma(nv, sd_v, corr_v)
  if (!is.vector(mean_v)) mean_v <- mean_v[,1]

  drifts <- tmvtnorm::rtmvnorm(n = n, mean = mean_v, sigma = Sigma,
    lower = rep(0, (nv)), algorithm = "gibbs")

  return(make_r(drifts, A, b, t0, st0, return_ttf))
}

# dca <- function(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
#   lower = rep(-Inf, length = length(mean)), upper = rep(Inf,
#   length = length(mean)), D = diag(length(mean)), H = NULL, ...) {
#
#   tmvtnorm::rtmvnorm(n, mean, sigma, lower, upper, D, H, "gibbs")
# }


##' Canonical Linear Ballistic Accumulation/Accumualtor Model
##'
##' \code{makeR} stands for making/generating/simulating responses from
##' a LBA model. \code{make_r} and \code{make.r} use C++ function. These
##' make \code{r}, \code{_r}, \code{.r} functions are essentially \code{rLBA},
##' including \code{rlba_norm}. They uses a LBA model with parameters, b, A,
##' mean_v, sd_v and t0 (no st0) to generate choice RT random deviates.
##'
##' \code{make_v} draws drift rate from normal or truncated normal distribution.
##' Each trial is stored as a row and each column is a drift rate for an
##' accumulator. You need to transpose drift rates generated by make_v for
##' \code{makeR}.
##'
##' \code{make.r} is a wrapper function of \code{make_r}. You may
##' need to use ":::" to call make.r, because of S3 method naming convention. If
##' you call \code{make_r} directly, beware it returns C index and is only a
##' numeric matrix. It does not carry a string vector for the column names, RTs
##' and responses. See timing test to see why it might be a good idea not to
##' return it as a data frame. \code{rlbaCnorm} is R version of correlated LBA
##' model.
##'
##' \code{rlba_norm} adds checks and calls \code{make_v} and \code{make_r}.
##' \code{rlba_norm} is only slightly quicker than \code{make_r}.
##'
##' \code{n1PDFfixedt0} is defective density function for the fisrt node LBA
##' model. Defective means its probability does not normally normalize to 1.
##' Only the probabilities from all nodes/accumulators add to 1.
##' \code{n1PDFfixedt0} is equation (3) on  page 159 in Brown and
##' Heathcote (2008).  This equation assumes st0 is 0.
##'
##' \code{fptcdf} and \code{fptpdf} are distribution and density functions with
##' four parameters A, b, v and sv, assuming t0 is zero. \code{fptcdf} and
##' \code{fptpdf} are respectively equation (1) and equation (2) on page 159 in
##' Brown and Heathcote (2008).
##'
##' @param drifts a n x n_v drift rate matrix. It can be a vector with 2 or more
##' elements. n is the numbers of observation. n_v is the numbers of
##' response/accumulator.
##' @param b decision threshold, a vector or a scalar.
##' @param A start point upper bound, a vector of a scalar.
##' @param n_v numbers of response/accumulator, an integer. Note n_v must match
##' the length/size of \code{drifts} vector.
##' @param t0 nondecision time, a vector or a scalar.
##' @param st0 nondecision time variation, a vector of a scalar. It is the upper
##' bound of a uniform distribution for t0 variability.
##' @param n numbers of observation/model simulations. This must be a scalar.
##' @param seed an integer specifying if and how the random number generator
##' should be initialized.
##' @param return_ttf a boolean switch indicating if return RTs for all
##' accumulators. When \code{return_ttf} is TRUE, a n_v x n ttf matrix is
##' returned.
##'
##' @return \code{make_r} gives either a time-to-finish (ttf) matrix or a n x 2
##' matrix, storing RTs (first column) and responses (second column). \code{n}
##' equals to number of model simulations. ttf is a n_v x n matrix with RTs from
##' all accumulators.
##
##' @export
maker <- function (drifts, n, b, A, n_v, t0, st0 = 0, seed = NULL,
  return_ttf = FALSE) {

  set.seed(seed)
  tmp <- make_r(drifts, A, b, t0, st0, return_ttf)

  if (return_ttf) {
    out <- tmp
  } else {
    out <- data.frame(RT = tmp[,1], R = tmp[,2])
  }
  attr(out, "seed") <- seed
  return(out)
}


##' Generate random variates of various cognitive models
##'
##' A wrapper function to \code{rd}, \code{norm}, \code{norm_pda},
##' \code{norm_pda_gpu}, \code{plba0_gpu}, \code{plba1}, \code{plba_gpu},
##' \code{plba2}, \code{plba3}, \code{lnr}, and \code{cnorm} models.
##'
##' @param type a character string indicating the model type
##' @param pmat a matrix of accumulator x parameter
##' @param n number of simulations
##' @param seed an integer specifying if and how the random number generator
##' should be initialized.
##' @param regressor independent variables in regression models
##' @importFrom rtdists rdiffusion
##' @export
random <- function(type, pmat, n, seed = NULL, regressors = NULL) {

  set.seed(seed)

  if (type == "rd") {
    out <- rtdists::rdiffusion(n, a = pmat$a[1], v = pmat$v[1],
      t0 = pmat$t0[1],
      z  = pmat$z[1]*pmat$a[1], # convert to absolute
      d  = pmat$d[1],
      sz = pmat$sz[1]*pmat$a[1],
      sv = pmat$sv[1], st0 = pmat$st0[1], stop_on_error = TRUE)

  } else if (type %in% c("norm", "norm_pda", "norm_pda_gpu")) {

    ## pmat: A b t0 mean_v sd_v st0
    # print(n)
    # print(pmat)
    out <- rlba_norm(n, pmat[, 1], pmat[, 2], matrix(pmat[, 4]), pmat[, 5],
      pmat[,3], pmat[1,6])

  } else if (type %in% c("plba0_gpu") ) {

    out <- rplba0(n, pmat[,1], pmat[,2], pmat[1,7], pmat[,3], pmat[,5],
      pmat[, 4], pmat[1,6], pmat[1,8])

  } else if (type %in% c("plba1", "plba1_gpu") ) {

    out <- rplba1(n, pmat[,1], pmat[,2], pmat[1,7], pmat[,3], pmat[,5],
      pmat[, 4], pmat[1,6], pmat[1,8])

  } else if (type == "plba2") {
    out <- rplba2(n, pmat[,1], pmat[,2], pmat[,3], pmat[,4],
      pmat[,5], pmat[,6], pmat[1, 7], pmat[1, 9], pmat[1, 8])

  } else if (type == "plba3") {
    out <- rplba3(n, pmat[,1], pmat[,2], pmat[,3],
      pmat[,4], pmat[,5], pmat[,6], pmat[,7], pmat[1, 8],
      pmat[1, 9], pmat[1, 11], pmat[1, 10])
  } else if (type == "lnr") {

    out <- rlnrDF(n, pmat[,1],  pmat[,2], pmat[,3], pmat[1,4])
  } else if (type == "cnorm") {

    out <- rca(n, pmat[,1], pmat[,2], pmat[,3], pmat[,4], pmat[,5], pmat[,7],
      unique(pmat[,6]))

  } else if (type == "glm") {
    # pmat is a data frame
    #       a   b   tau
    # r1 143 -1.2   8.8

    # X_  <- sample(x = regressors, size = n, replace = TRUE)
    # mu <- pmat[1,1] + pmat[1,2] * X_
    # sd <- 1/sqrt(pmat[1,3])
    # out <- cbind(rep(1, n), X_, rnorm(n, mu, sd))

    nbeta <- ncol(pmat) - 1
    beta <- matrix(unlist(pmat[1, 1:nbeta]), ncol = 1)
    tau  <- pmat[1, ncol(pmat)]
    X_ <- sample(x = regressors, size = n, replace = TRUE)
    X  <- as.matrix(cbind(rep(1, n), X_))
    linpred <- X %*% beta
    sd <- 1/sqrt(tau)
    out <- cbind(X, rnorm(n, linpred, sd))

  } else if (type == "logit") {
    nbeta <- ncol(pmat) - 1
    beta <- matrix(unlist(pmat[1, 1:nbeta]), ncol = 1)
    X  <- as.matrix(cbind(rep(1, n), regressors,
                    regressors[,1] * regressors[,2]))

    b <- rnorm(n, 0, pmat[1, nbeta+1])
    linpred <- X %*% beta + b

    prob <- plogis(linpred, 0, 1, TRUE, FALSE)
    ni <- sample(1:100, n, replace = TRUE)
    out <- cbind(rep(1, n), rbinom(n, size = ni, prob = prob), ni)

  } else {
    stop("Model type yet created")
  }

  attr(out, "seed") <- seed
  return(out)
}

##' Return ns-npar matrix
##'
##' Contructs a ns x npar matrix, indicating the true paramters
##' used to simualte data. Each row represents a set of parameters for a
##' participant. One should enter either a valid vector or matrix for
##' true parameters (i.e., ps) or a list of (parameter) prior distributions
##' (p.prior). When \code{p.prior} is supplied, true parameters are drawn
##' from prior distributions.
##'
##' @param x a model object
##' @param ns number of subjects.
##' @param prior a list of parameter prior distributions
##' @param ps a vector or matirx. Each row indicates a set of true parameters
##' for a participant.
##' @param seed an integer specifying if and how the random number generator
##' should be initialized.
##' @examples
##' model <- BuildModel(
##' p.map     = list(a ="1", v = "1",z = "1", d = "1", sz = "1", sv = "1",
##'             t0 = "1", st0 = "1"),
##' match.map = list(M = list(s1 = "r1", s2 ="r2")),
##' factors   = list(S = c("s1", "s2")),
##' constants = c(st0 = 0, d = 0),
##' responses = c("r1", "r2"),
##' type      = "rd")
##'
##' p.prior <- BuildPrior(
##'   dists = c("tnorm", "tnorm", "beta", "beta", "tnorm", "beta"),
##'   p1    = c(a = 1, v = 0, z = 1, sz = 1, sv = 1, t0 = 1),
##'   p2    = c(a = 1, v = 2, z = 1, sz = 1, sv = 1, t0 = 1),
##'   lower = c(0, -5, NA, NA, 0, NA),
##'   upper = c(2,  5, NA, NA, 2, NA))
##'
##' ## Example 1: Randomly generate 2 sets of true parameters from
##' ## parameter priors (p.prior)
##' GetParameterMatrix(model, 2, p.prior)
##' ##            a         v         z        sz       sv        t0
##' ## [1,] 1.963067  1.472940 0.9509158 0.5145047 1.344705 0.0850591
##' ## [2,] 1.512276 -1.995631 0.6981290 0.2626882 1.867853 0.1552828
##'
##' ## Example 2: Use a user-selected true parameters
##' true.vector  <- c(a=1, v=1, z=0.5, sz=0.2, sv=1, t0=.15)
##' GetParameterMatrix(model, 2, NA, true.vector)
##' ##      a v   z  sz sv   t0
##' ## [1,] 1 1 0.5 0.2  1 0.15
##' ## [2,] 1 1 0.5 0.2  1 0.15
##' GetParameterMatrix(model, 2, ps = true.vector)
##'
##' ## Example 3: When a user enter arbritary sequence of parameters.
##' ## Note sv is before sz. It should be sz before sv
##' ## See correct sequence, by entering "attr(model, 'p.vector')"
##' ## GetParameterMatrix will rearrange the sequence.
##' true.vector  <- c(a=1, v=1, z=0.5, sv=1, sz = .2, t0=.15)
##' GetParameterMatrix(model, 2, NA, true.vector)
##' ##      a v   z  sz sv   t0
##' ## [1,] 1 1 0.5 0.2  1 0.15
##' ## [2,] 1 1 0.5 0.2  1 0.15
##'
##' @export
GetParameterMatrix <- function(x, ns, prior = NA, ps = NA, seed = NULL) {

  message1 <- "Parameters are incompatible with model"
  # message2 <- "Must supply either a list of p.prior or a parameter vector."
  # if(anyNA(prior) & anyNA(ps)) stop(message2)
  pnames <- names(attr(x, "p.vector"))

  if (anyNA(prior)) { ## use ps
    if (is.vector(ps)) {
      if (check_pvec(ps, x)) stop(message1)
      ps    <- ps[pnames]
      pss   <- rep(ps, each = ns)
      psmat <- matrix(pss, ns, dimnames = list(NULL, pnames))
    } else if (is.matrix(ps)) {
      psmat <- matrix(ps, ns, dimnames = list(NULL, pnames))
    } else {
      if ((ns != dim(ps)[1])) stop("ps matrix must have ns rows")
      if (check_pvec(ps[1,], x)) stop(message1)
    }

    rownames(psmat) <- 1:ns
  } else {  ## use prior; random-effect model
    if (!all( pnames %in% names(prior))) stop(message1)
    set.seed(seed)
    psmat <- rprior(prior[pnames], ns)

    ## A nasty way to deal with MG's error gate; ie keep redrawing until we
    ## pass his checks.
    # if (attr(x, "type") == "rd") {
    #   facs <- ggdmc::createfacsdf(x)
    #
    #   for (i in 1:ns) {
    #     for (j in 1:nrow(facs)) {
    #       psmat_allpar <- TableParameters(psmat[i,], j, x, FALSE)
    #       psmat_allpar <- checkddm3(psmat_allpar, j, x, prior)
    #       psmat[i,] <- as.numeric(psmat_allpar[1, pnames])
    #     }
    #   }
    # }
  }

  return(psmat)
}

GetParameterMatrix_bk <- function(x, ns, prior = NA, ps = NA, seed = NULL) {

  message1 <- "Parameters are incompatible with model"
  pnames <- names(attr(x, "p.vector"))

  if (anyNA(prior)) { ## use ps
    if (is.vector(ps)) {
      if (check_pvec(ps, x)) stop(message1)
      ps    <- ps[pnames]
      pss   <- rep(ps, each = ns)
      psmat <- matrix(pss, ns, dimnames = list(NULL, pnames))
    } else if (is.matrix(ps)) {
      psmat <- matrix(ps, ns, dimnames = list(NULL, pnames))
    } else {
      if ((ns != dim(ps)[1])) stop("ps matrix must have ns rows")
      if (check_pvec(ps[1,], x)) stop(message1)
    }

    rownames(psmat) <- 1:ns
  } else {  ## use prior; random-effect model
    if (!all( pnames %in% names(prior))) stop(message1)
    set.seed(seed)
    psmat <- rprior(prior[pnames], ns)
  }

  return(psmat)
}

# ps <- c(a = 242.7, b = 6.185, tau = .01)
# ps <- c(b0 = -.55, b1 = .08, b2 = -.81, b3 = 1.35, sd = .267)
# ps <- c(b1 = .08, b2 = -.81, b3 = 1.35, sd = .267)
# n <- 10
# seed <- NULL
simulate_one <- function(model, n, ps, seed) {
  if (check_pvec(ps, model)) stop("p.vector and model incompatible")
  resp <- attr(model, "responses")
  type <- attr(model, "type")
  levs <- attr(model, "factors")
  facs <- ggdmc:::createfacsdf(model)
  nvec <- ggdmc:::check_n(n, facs)
  dat  <- ggdmc:::nadf(nvec, facs, levs, type)
  row1 <- 1

  dfnames <- names(dat)
  Xnames   <- dfnames[!dfnames %in% c("R", "N", "Y")]

  # i <- 2
  for (i in 1:nrow(facs)) {
    pmat <- TableParameters(ps, i, model, FALSE) ## simulate use n1.order == FALSE
    rown <- row1 + nvec[i] - 1

    if (type == "glm") {
      regressors <- attr(model, "regressors")
      tmp <- random(type, pmat, nvec[i], seed, regressors)
      dat[row1:rown, c("R", "X", "Y")] <- tmp
    } else if (type == "logit") {
      regressors <- dat[row1:rown, Xnames]
      tmp <- random(type, pmat, nvec[i], seed, regressors)
      dat[row1:rown, c("R", "Y", "N")] <- tmp
    } else {
      dat[row1:rown, c("RT", "R")] <- random(type, pmat, nvec[i], seed)
    }
    row1 <- rown+1
  }

  dat$R <- factor(dat$R, levels = 1:length(resp), labels = resp)
  if (type == "rd") dat <- FlipResponse_rd(model, dat, facs)
  return(dat)
}

simulate_many <- function(model, n, ns, prior, ps, seed) {

  n  <- GetNsim(model, n, ns)
  ps <- GetParameterMatrix(model, ns, prior, ps, seed)

  ismanysub <- ismanymodels(model, ns)
  if(ismanysub) modeli <- model[[1]] else modeli <- model

  ndatai <- cumsum(c(0, matrixStats::rowSums2(n))); ## index boundaries
  datr <- (ndatai[1] + 1):(ndatai[2]); ## First subject's trial index
  ## Simulate first subject; modeli must be 'model' class
  dat <- cbind(s = rep.int(1, length(datr)),
    simulate_one(modeli, n[1,], ps[1,], seed))

  if (ns > 1) {
    for (i in 2:ns) {
      if (ismanysub) modeli <- model[[i]] else modeli <- model
      datr <- (ndatai[i] + 1):(ndatai[i + 1]) ## continue to index trials
      dat  <- rbind(dat,
        cbind(s = rep.int(i, length(datr)),
          simulate_one(modeli, n[i,], ps[i,], seed)))
    }
  }

  dat$s <- factor(dat$s)
  attr(dat, "parameters") <- ps
  ## if ps is not auto-created by p.prior, save the user's p.prior in 'attribute'
  if(!anyNA(prior)) attr(dat, "p.prior") <- prior
  return(dat)
}

##' Simulate RT Data
##'
##' Simulate stochastic responses either for one subject or multiple subjects.
##' The simulation is based on the \code{model} object. For one subject, the
##' user must supply true parameters, \code{p.vector} at \code{ps} argument.
##' For multiple subjects, the user can supply a matrix (or a row vector),
##' indicating true parameters for each subject, separately on each row
##' (via \code{ps} argument). This is the fixed-effect model. If the user
##' wants to simulate from a random-effect (i.e., hierarchical) model, in which
##' case p.prior must be supplied and ps will be ignored. Note in some cases,
##' a random-effect model may fail to draw data from the model, because
##' true parameters are drawn from \code{p.prior} and a specific model, like
##' DDM, may has certain ranges from different parameters.
##'
##' \code{ps} can be a row vector, in which case each subject has identical
##' parameters. It can also be a matrix with one row per subject, in which
##' case it must have \code{ns} rows. The true values will be saved as
##' "parameters" attribute.
##'
##' @param object a model object.
##' @param nsim number of trials/responses. \code{n} can be a single number for a
##' balanced design or matrix for an unbalanced design, where rows are
##' subjects and columns are design cells. If the matrix has one row then all
##' subjects have the same \code{n} in each cell, if it has one column then all
##' cells have the same \code{n}; Otherwise each entry specifies the \code{n}
##' for a particular design subject x design cell combination.
##' @param nsub number of subjects
##' @param prior parameter priors. A list of distributions based on which
##' the true parameters fro each subject are drawn.  It is usually created by
##' \code{BuildPrior} and will be saved as "p.prior" attribute.
##' @param ps p.vector matrix. Each row represent a subject.
##' @param seed an integer specifying if and how the random number generator
##' should be initialized.
##' @param ... additional optional arguments.
##' @return a data frame
##' @importFrom stats simulate
##' @examples
##' model <- BuildModel(
##'   p.map     = list(a = "1", v = "1", z = "1", d = "1", sz = "1",
##'   sv = "1", t0 = "1", st0 = "1"),
##'   match.map = list(M = list(s1 = "r1", s2 = "r2")),
##'   factors   = list(S = c("s1", "s2")),
##'   constants = c(st0 = 0, d = 0),
##'   responses = c("r1", "r2"),
##'   type      = "rd")
##'
##'
##' @export
simulate.model <- function(object, nsim = NA, seed = NULL, nsub = NA,
  prior = NA, ps = NA, ...) {
  if (is.na(nsub)) {
    if (is.na(nsim)) stop("How many response you want to generate? Must supply n")
    if (anyNA(ps)) {
      ps <- GetParameterMatrix(object, 1, prior, seed)
      ## stop("Some true parameters missing")
    }
    out <- simulate_one(object, nsim, ps, seed)
    attr(out, "parameters") <- ps
  } else {
    message1 <- "Must supply either sampling distribution or a true vector."
    if (anyNA(prior) & anyNA(ps)) stop(message1)
    out <- simulate_many(object, nsim, nsub, prior, ps, seed)
  }
  return(out)
}

##' Post-predictive Simulation
##'
##' Simulate post-predictive data.
##'
##' @param object a model object of a subject.
##' @param npost number of posterior predictive replications.
##' @param nsub whether randomly select a npost estimated parameter values or
##' select the first npost estimates
##' @param factors experimental factors
##' @param xlim trimming outlier, e.g., xlim = c(0, 5).
##' @param seed an integer specifying if and how the random number generator
##' should be initialized.
##' @return a data frame
##'
##' @export
predict_one <- function(object, npost = 100, rand = TRUE, factors = NA,
                        xlim = NA, seed = NULL) {

  model <- attributes(object$data)$model
  facs <- names(attr(model, "factors"))
  class(object$data) <- c("data.frame", "list")

  if (!is.null(factors)) {
    if (any(is.na(factors))) factors <- facs
    if (!all(factors %in% facs))
      stop(paste("Factors argument must contain one or more of:",
                 paste(facs, collapse=",")))
  }

  resp <- names(attr(model, "responses"))
  ns   <- table(object$data[,facs], dnn = facs)
  npar   <- object$n.pars
  nchain <- object$n.chains
  nmc    <- object$nmc
  ntsample <- nchain * nmc
  pnames   <- object$p.names
  # str(samples$theta) ## nchain x npar x nmc
  # str(thetas)        ## (nchainx nmc) x npar
  thetas <- matrix(aperm(object$theta, c(3,1,2)), ncol = npar)
  colnames(thetas) <- pnames

  if (is.na(npost)) {
    use <- 1:ntsample
  } else {
    if (rand) {
      use <- sample(1:ntsample, npost, replace = F)
    } else {
      use <- round(seq(1, ntsample, length.out = npost))
    }
  }

  npost  <- length(use)
  posts   <- thetas[use, ]
  nttrial <- sum(ns) ## number of total trials
  out <- data.frame(matrix(nrow = npost*nttrial, ncol = dim(object$data)[2]))

  # for (i in 1:npost) {
  #   tmp <- ggdmc:::simulate_one(model, n = ns, ps = posts[i,], seed = NULL)
  #   out[(1+(i-1)*nttrial):(i*nttrial), names(tmp)] <- tmp
  #   if ( (i %% report) == 0) cat(".")
  # }

  ## should replace with parallel
  v <- lapply(1:npost, function(i) {
    ggdmc:::simulate_one(model, n = ns, ps = posts[i,], seed = seed)
  })
  out <- data.table::rbindlist(v)
  # names(out) <- names(object$data)
  reps <- rep(1:npost, each = nttrial)
  out <- cbind(reps, out)

  if (!any(is.na(xlim))) {
    out <- out[RT > xlim[1] & RT < xlim[2]]
  }

  attr(out, "data") <- object$data
  return(out)
}

"%w/o%" <- function(x, y) x[!x %in% y] #--  x without y
