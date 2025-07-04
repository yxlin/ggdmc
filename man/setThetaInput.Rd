% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model-class.R
\name{setThetaInput}
\alias{setThetaInput}
\title{Create a theta_input Object for MCMC Configuration}
\usage{
setThetaInput(
  nmc = 500L,
  nchain = NULL,
  thin = 1L,
  pnames = NULL,
  report_length = 100L,
  max_init_attempts = 1000L,
  is_print = TRUE
)
}
\arguments{
\item{nmc}{Integer specifying total number of MCMC iterations (including burn-in).
Default: 500L}

\item{nchain}{Integer specifying number of MCMC chains facilitating the
DE-MC sampling. If NULL (default), it will be set to 3 times the number of
parameters.}

\item{thin}{Integer specifying thinning interval (keep every nth sample).
Default: 1L (no thinning)}

\item{pnames}{Character vector naming the free parameters. Default: NULL}

\item{report_length}{Integer controlling how often progress is reported.
Default: 100L (report every 100 iterations)}

\item{max_init_attempts}{Integer specifying maximum attempts to initialise
valid starting values. Default: 1000L}

\item{is_print}{Logical controlling whether to print sampling progress.
Default: TRUE}
}
\value{
A \code{theta_input} S4 object containing MCMC configuration
  parameters with the following slots:
  \itemize{
    \item \code{nmc}: Total iterations
    \item \code{nchain}: Number of chains
    \item \code{thin}: Thinning interval
    \item \code{nparameter}: Number of parameters (derived from \code{pnames})
    \item \code{pnames}: Parameter names
    \item \code{report_length}: Progress report frequency
    \item \code{max_init_attempts}: Initialisation attempts
    \item \code{is_print}: Print progress flag
  }
}
\description{
Constructs a \code{theta_input} instance storing the parameter that controls
Markov Chain Monte Carlo (MCMC) sampling.
}
\details{
This constructor:
\itemize{
  \item Validates all input parameters
  \item Automatically sets \code{nchain} to 3 × number of parameters if NULL
  \item Derives \code{nparameter} from length of \code{pnames}
  \item Ensures all numeric parameters are converted to integers
}
}
\examples{
# A minimal LBA model with 5 free parameters
nmc <- 1000L
npar <- 5L
nchain <- npar * 3
thin <- 1L
pnames <- c("A", "B", "mean_v.false", "mean_v.true", "t0")

sub_theta_input <- setThetaInput(
    nmc = nmc, nchain = as.integer(nchain), pnames = pnames,
    thin = thin
)

# theta_input has its own print method
print(sub_theta_input)

}
