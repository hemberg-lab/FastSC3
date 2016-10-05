#' Fast-Johnson-Lindenstrauss-Transform (FJLT)
#'
#' This function calculates the FJLT and project onto d dimensions. The FJLT 
#' is faster than standard random pro-jections and just as easy to implement. 
#' It is based upon the preconditioning of a sparse projection matrix with a 
#' randomized Fourier transform.
#' 
#' Ailon, N. & Chazelle, B. Approximate nearest neighbors and the fast 
#' Johnson-Lindenstrauss transform. in Proceedings of the thirty-eighth 
#' annual ACM symposium on Theory of computing 557–563 (ACM, 2006).
#' 
#' Functions adapted from: http://www.cs.ubc.ca/~jaquesn/MachineLearningTheory.pdf
#'
#' @param x input expression matrix
#' @param k number of dimension to reduce to
#' 
#' @useDynLib FastSC3
#' 
#' @return transformed reduced expression matrix
#' @export
fjlt <- function(x, k = 100) {
    # pad data to a power of 2
    i <- 1
    while(2^i < nrow(x)){
        i <- i + 1
    }
    # return FJLT of x
    # p - the norm we are using (Euclidean, p = 2)
    calc_fjlt(x = x, p = 2, k = k, d = 2^i, n = ncol(x))
}


#' Distance matrix transformation
#'
#' All distance matrices are transformed to either covariance matrix 
#' (in case of PCA) or to graph Laplacian (in case of Laplacian).
#'
#' @param dists distance matrix
#' @param method transformation method: either 'pca' or
#' 'laplacian'
#' @return transformed distance matrix
#' 
#' @importFrom SC3 norm_laplacian
ftransformation <- function(dists, method) {
    if (method == "pca") {
        covar <- cov(scale(dists, center = TRUE, scale = TRUE))
        return(covar)
    } else if (method == "laplacian") {
        return(norm_laplacian(dists))
    }
}
