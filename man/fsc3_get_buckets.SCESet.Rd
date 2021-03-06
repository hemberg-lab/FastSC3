% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoreMethods.R
\docType{methods}
\name{fsc3_get_buckets.SCESet}
\alias{fsc3_get_buckets}
\alias{fsc3_get_buckets,SCESet-method}
\alias{fsc3_get_buckets.SCESet}
\title{Define buckets of cells based on their signatures}
\usage{
fsc3_get_buckets.SCESet(object, common_bits = NULL, runs = 50)

\S4method{fsc3_get_buckets}{SCESet}(object, common_bits = NULL, runs = 50)
}
\arguments{
\item{object}{an object of 'SCESet' class}

\item{common_bits}{number of common bits which is enough to put two signatures
in the same bucket}
}
\value{
an object of 'SCESet' class
}
\description{
Two cells are put in the same buckets if their binary signatures have a 
specific number of common bits.
}

