% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{dataframe2list}
\alias{dataframe2list}
\title{converts multivariate dataframe to list}
\usage{
dataframe2list(x, R.res, Trait, Pop)
}
\arguments{
\item{x}{Data frame or list containing summary statistics for
multiple parameters measured in both sexes in two or more populations.}

\item{R.res}{Pooled within correlational matrix, Default: NULL}

\item{Trait}{Number of the column containing names of measured parameters,
Default: 1}

\item{Pop}{Number of the column containing populations' names, Default: 2}
}
\description{
helper function for multivariate analysis
}
\keyword{internal}
