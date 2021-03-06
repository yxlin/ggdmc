% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model-class.R
\docType{class}
\name{dmi-class}
\alias{dmi-class}
\title{An S4 class of the Data-model Instance}
\description{
The class is to represent a data-model instance, which joins a model object
with a data frame. The process of BuildDMI also generates cell.index and
cell.empty.
}
\section{Slots}{

\describe{
\item{\code{data}}{A data frame storing the would-be fit data set}

\item{\code{model}}{A 3-D model array.  Dimension one stores the combinations
of the factor levels and response types, dimension two stores parameters,
and dimension three stores response types.}

\item{\code{cell.index}}{A ncell-element list. Each element represents one cell.
Each element stores \code{nobs} Boolean indicators, showing whether a
particular observation belongs to this cell.}

\item{\code{cell.empty}}{A ncell-element logical vector, indicating whether this
cell has no observation.}
}}

