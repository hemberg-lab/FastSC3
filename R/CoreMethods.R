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

#' Calculate transformations of the distance matrices.
#' 
#' This function calculates transforamtions of the distance matrices contained in 
#' the sc3_distances item of the object@sc3 slot. It then
#' creates and populates the following items of the object@sc3 slot:
#' \itemize{
#'   \item transformations - contains a list of transformations of the 
#'   distance matrices corresponding to covariance and graph laplacian.
#' }
#' 
#' @param object an object of 'SCESet' class
#' 
#' @return an object of 'SCESet' class
#' 
#' @importFrom doRNG %dorng%
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' 
#' @export
fsc3_calc_transfs.SCESet <- function(object) {
    dists <- object@sc3$distances
    if (is.null(dists)) {
        warning(paste0("Please run sc3_calc_dists() first!"))
        return(object)
    }
    
    # NULLing the variables to avoid notes in R CMD CHECK
    i <- NULL
    
    distances <- names(dists)
    transformations <- c("pca", "laplacian")
    
    hash.table <- expand.grid(dists = distances, transfs = transformations, stringsAsFactors = FALSE)
    
    message("Performing transformations...")
    
    if (object@sc3$n_cores > nrow(hash.table)) {
        n.cores <- nrow(hash.table)
    } else {
        n.cores <- object@sc3$n_cores
    }
    
    cl <- parallel::makeCluster(n.cores, outfile = "")
    doParallel::registerDoParallel(cl, cores = n.cores)
    
    # calculate the 6 distinct transformations in parallel
    transfs <- foreach::foreach(i = 1:nrow(hash.table)) %dopar% {
        try({
            ftransformation(get(hash.table[i, 1], dists), hash.table[i, 2])
        })
    }
    
    # stop local cluster
    parallel::stopCluster(cl)
    
    names(transfs) <- paste(hash.table[, 1], hash.table[, 2], sep = "_")
    
    object@sc3$transformations <- transfs
    return(object)
}

#' @rdname fsc3_calc_transfs.SCESet
#' @aliases fsc3_calc_transfs
#' @importClassesFrom scater SCESet
#' @export
setMethod("fsc3_calc_transfs", signature(object = "SCESet"), function(object) {
    fsc3_calc_transfs.SCESet(object)
})

#' Calculate transformations of the distance matrices.
#' 
#' This function calculates transforamtions of the distance matrices contained in 
#' the sc3_distances item of the object@sc3 slot. It then
#' creates and populates the following items of the object@sc3 slot:
#' \itemize{
#'   \item transformations - contains a list of transformations of the 
#'   distance matrices corresponding to covariance and graph laplacian.
#' }
#' 
#' @param object an object of 'SCESet' class
#' 
#' @return an object of 'SCESet' class
#' 
#' @importFrom doRNG %dorng%
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' 
#' @export
fsc3_calc_eigenv.SCESet <- function(object) {
    transfs <- object@sc3$transformations
    if (is.null(transfs)) {
        warning(paste0("Please run fsc3_calc_transfs() first!"))
        return(object)
    }
    
    # NULLing the variables to avoid notes in R CMD CHECK
    i <- NULL
    
    transformations <- names(transfs)
    n.dim <- max(object@sc3$n_dim)

    message("Calculating eigenvectors...")
    
    if (object@sc3$n_cores > length(transformations)) {
        n.cores <- length(transformations)
    } else {
        n.cores <- object@sc3$n_cores
    }
    
    cl <- parallel::makeCluster(n.cores, outfile = "")
    doParallel::registerDoParallel(cl, cores = n.cores)
    
    # calculate the 6 distinct transformations in parallel
    eigenv <- foreach::foreach(i = 1:length(transformations)) %dopar% {
        try({
            ssNystrom(transfs[[i]], r = n.dim)
        })
    }
    
    # stop local cluster
    parallel::stopCluster(cl)
    
    names(eigenv) <- transformations
    
    object@sc3$transformations <- eigenv
    return(object)
}

#' @rdname fsc3_calc_eigenv.SCESet
#' @aliases fsc3_calc_eigenv
#' @importClassesFrom scater SCESet
#' @export
setMethod("fsc3_calc_eigenv", signature(object = "SCESet"), function(object) {
    fsc3_calc_eigenv.SCESet(object)
})
