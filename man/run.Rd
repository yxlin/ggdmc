% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{run_subject}
\alias{run_subject}
\alias{run_hyper}
\alias{run}
\title{Run MCMC Sampling for Cognitive Models}
\usage{
run_subject(config_r, dmi, samples)

run_hyper(config_r, dmi, samples)

run(config_r, dmis, samples)
}
\arguments{
\item{config_r}{An S4 object of class \code{config} specifying model
settings, priors, tuning parameters, etc.}

\item{dmi}{For \code{run_subject()} and \code{run_hyper()}, a single S4
object containing a data-model instance for one subject or one set of
pressumed true theta values (in a R matrix format). For \code{run()}, a
list of \code{dmi} objects (one per subject).}

\item{samples}{For \code{run_subject()} and \code{run_hyper()}, an
S4 \code{samples} object containing starting parameter values.
For \code{run()}, a list of such \code{samples} objects for multiple
subjects and a phi `samples` for the hyper-level.}

\item{dmis}{For \code{run()}, a list of \code{dmi} objects (one per
subject).}
}
\value{
\itemize{
  \item \code{run_subject()} returns an S4 \code{posterior} object
for a single subject.
  \item \code{run_hyper()} returns an S4 \code{posterior} object
         for hyperparameters (group-level parameters).
  \item \code{run()} returns a named \code{list} of \code{posterior}
          objects, one per subject.
}
}
\description{
These functions manage the MCMC sampling process for Linear
Ballistic Accumulator (LBA) models via C++ backends. The sampling is
configured using an S4 `config` object, which carries sampler settings
and model specifications. The `dmi` (data-model instance) contains
trial-level data and likelihood structure, and `samples` represents an
initial or resumed parameter state.
}
\details{
Executes MCMC-based posterior sampling for LBA or decision-diffusion
models. The model can be structured as either one or two level(s),
including subject-only, group-level, and hierarchical
multi-subject models.


These functions are typically called from higher-level model fitting
workflows. Internally, they coordinate log-likelihood evaluation,
parameter proposal, acceptance checks, and diagnostics.

\code{run_subject()} and \code{run_hyper()} differ in that the former
operates on trial-level likelihoods for an individual subject, while the
latter works with summary statistics and priors for group-level parameters.

\code{run()} enables full hierarchical modeling across multiple subjects
by looping through each subject's `dmi` and `samples`,
applying \code{run_subject()}.
}
\examples{
\dontrun{
# Setup (assuming valid config_r, dmi, and samples objects exist)
posterior_subject <- run_subject(config_r, dmi, samples)

# Group-level inference
posterior_hyper <- run_hyper(config_r, dmi, samples)

# Hierarchical model
posteriors <- run(config_r, dmis = list(dmi1, dmi2), samples = list(s1, s2))
}

}
\keyword{internal}
