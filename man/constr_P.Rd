% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{constr_P}
\alias{constr_P}
\title{Construct the P matrix for the FJLT transform}
\usage{
constr_P(p, k, d, n)
}
\arguments{
\item{p}{the norm}

\item{k}{dimension we reduce to}

\item{d}{dimension of the input matrix}

\item{n}{number of cells (columns) in the input matrix}
}
\value{
matrix P of the FJLT transform
}
\description{
Construct the P matrix for the FJLT transform
}

