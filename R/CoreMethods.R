#' Fast-Johnson-Lindenstrauss-Transform (FJLT)
#'
#' This function calculates the FJLT and project onto d dimensions. The FJLT 
#' is faster than standard random projections and just as easy to implement. 
#' It is based upon the preconditioning of a sparse projection matrix with a 
#' randomized Fourier transform.
#' 
#' Ailon, N. & Chazelle, B. Approximate nearest neighbors and the fast 
#' Johnson-Lindenstrauss transform. in Proceedings of the thirty-eighth 
#' annual ACM symposium on Theory of computing 557–563 (ACM, 2006).
#' 
#' Calculate Fast-Johnson-Lindenstrauss-Transform (FJLT)
#' 
#' This function calculates FJLT projection of
#' the processed_dataset item of the object@sc3 slot. It then
#' creates and populates the following items of the object@sc3 slot:
#' \itemize{
#'   \item processed_dataset_fjlt - contains the FJLT transform of the processed_dataset
#' }
#'
#' @param object an object of 'SCESet' class
#' @param k number of dimensions to reduce to
#' 
#' @export
fsc3_fjlt.SCESet <- function(object, k = 100) {
    dataset <- object@sc3$processed_dataset
    if (is.null(dataset)) {
        warning(paste0("Please run sc3_prepare() first!"))
        return(object)
    }
    object@sc3$processed_dataset <- fjlt(dataset, k)
    return(object)
}

#' @rdname fsc3_fjlt.SCESet
#' @aliases fsc3_fjlt
#' @export
setMethod("fsc3_fjlt", signature(object = "SCESet"), function(object, 
                                                k = 100) {
    fsc3_fjlt.SCESet(object, k)
})
