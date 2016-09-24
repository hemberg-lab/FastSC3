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
#' @param epsilon error tolerance parameter
#' 
#' @useDynLib FastSC3
#' 
#' @return transformed reduced expression matrix
#' @export
FJLT <- function(x, k = 100) {
    # p - the norm we are using (Euclidean, p = 2)
    p <- 2
    n <- ncol(x)
    # pad data to a power of 2
    d_orig <- nrow(x)
    d <- findClosestPowerOf2(d_orig)
    # pad matrix with 0 if d > d_orig
    x <- rbind(x, matrix(0, nrow = d - d_orig, ncol = n))
    # compute FJLT
    res <- calculateFJLT(x, p, k, d, n)
    return(res)
}

findClosestPowerOf2 <- function(d) {
    i <- 1
    repeat{
        if(2^i >= d){
            return(2^i)
        }
        i <- i + 1
    }
}
