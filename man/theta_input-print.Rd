% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model-class.R
\docType{methods}
\name{theta_input-print}
\alias{theta_input-print}
\alias{print,theta_input-method}
\title{Print Method for theta_input Objects}
\usage{
\S4method{print}{theta_input}(x)
}
\arguments{
\item{x}{A \code{theta_input} object to be printed}
}
\value{
Invisibly returns the input object \code{x}. Called for its side effect
  of printing to the console.
}
\description{
Displays a human-readable summary of a \code{theta_input} object's contents.
}
\details{
This method provides a concise console representation of \code{theta_input} objects,
showing all key configuration parameters in a standardised format:
\itemize{
  \item Basic MCMC dimensions (iterations, chains, thinning)
  \item Parameter information (count and names)
  \item Runtime control settings (reporting, initialization attempts)
  \item Operational flags (printing, accumulation)
}

The output is designed to give users immediate visibility into the MCMC
configuration while avoiding overwhelming technical details.
}
\examples{
# Create and print a theta_input object
theta <- setThetaInput(
    nmc = 1000L,
    nchain = 3L,
    pnames = c("alpha", "beta", "sigma")
)
print(theta)

# typing the object name to use the default print method
theta
}
