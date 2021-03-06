% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoreMethods.R
\docType{methods}
\name{fsc3_select_repr_cells.SCESet}
\alias{fsc3_select_repr_cells}
\alias{fsc3_select_repr_cells,SCESet-method}
\alias{fsc3_select_repr_cells.SCESet}
\title{Reduce the original dataset by using cell buckets}
\usage{
fsc3_select_repr_cells.SCESet(object, d.region.min = 0.04,
  d.region.max = 0.07)

\S4method{fsc3_select_repr_cells}{SCESet}(object, d.region.min = 0.04,
  d.region.max = 0.07)
}
\arguments{
\item{object}{an object of 'SCESet' class}

\item{d.region.min}{defines the minimum number of eigenvectors used for 
kmeans clustering as a fraction of the total number of cells. Default is 0.04.}

\item{d.region.max}{defines the maximum number of eigenvectors used for 
kmeans clustering as a fraction of the total number of cells. Default is 0.07.}
}
\value{
an object of 'SCESet' class
}
\description{
Reduce the original dataset by using cell buckets
}

