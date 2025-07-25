% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model-class.R
\name{.compute_psrf_components}
\alias{.compute_psrf_components}
\title{Compute PSRF components}
\usage{
.compute_psrf_components(W, B, niter, nchain, conf)
}
\arguments{
\item{W}{Within-chain covariance matrix}

\item{B}{Between-chain covariance matrix}

\item{niter}{Number of iterations}

\item{nchain}{Number of chains}

\item{conf}{Confidence level for upper CI (default = 0.95)}
}
\value{
List with R-hat components
}
\description{
Internal function for R-hat calculations
}
\keyword{internal}
