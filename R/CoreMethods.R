#' Fast-Johnson-Lindenstrauss-Transform (FJLT)
#'
#' This function calculates the FJLT and project onto d dimensions. The FJLT 
#' is faster than standard random projections and just as easy to implement. 
#' It is based upon the preconditioning of a sparse projection matrix with a 
#' randomized Fourier transform.
#' 
#' Ailon, N. and Chazelle, B. Approximate nearest neighbors and the fast 
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
    dataset <- SC3::get_processed_dataset(object)
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

#' Calculate eigenvectors of the transformed distance matrix.
#' 
#' This function calculates eigenvectors of the transformed distance matrices contained in 
#' the transformation item of the object@sc3 slot. It then
#' populates back the transformations item the calculated eigenvectors.
#' 
#' Eigenvectors are calculated by using the SS-Nystrom approximation.
#' Wang, S., Zhang, C., Qian, H. and Zhang, Z. Improving the Modified NystroM 
#' Method Using Spectral Shifting. in Proceedings of the 20th ACM SIGKDD 
#' International Conference on Knowledge Discovery and Data Mining 611–620 
#' (ACM, 2014).
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
            ssNystrom(transfs[[i]], c = n.dim)
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

#' Calculate doubly-stochastic approximation of the PSD matrix
#' 
#' Zass, R. and Shashua, A. in Advances in Neural Information Processing 
#' Systems 19 (eds. Scholkopf, B., Platt, J. C. and Hoffman, T.) 1569-1576 
#' (MIT Press, 2007).
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
fsc3_norm_kernel.SCESet <- function(object) {
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
    normals <- foreach::foreach(i = 1:length(transformations)) %dopar% {
        try({
            normalise_kernel(transfs[[i]])
        })
    }
    
    # stop local cluster
    parallel::stopCluster(cl)
    
    names(normals) <- transformations
    
    object@sc3$transformations <- normals
    return(object)
}

#' @rdname fsc3_norm_kernel.SCESet
#' @aliases fsc3_norm_kernel
#' @importClassesFrom scater SCESet
#' @export
setMethod("fsc3_norm_kernel", signature(object = "SCESet"), function(object) {
    fsc3_norm_kernel.SCESet(object)
})

#' Define a set of genes used for creating cell binary signatures
#' 
#' The important genes (hyperplanes) are defined as differentially expressed 
#' genes identified using M3Drop (Michaelis-Menten Modelling of Dropouts for scRNASeq,
#' http://bioconductor.org/packages/M3Drop). A 'data.frame' with the genes
#' and corresponding p and q values is written to the 'hyperplanes' item of the
#' 'sc3' slot of the input object.
#' 
#' @param object an object of 'SCESet' class
#' @param n_genes number of the genes to be returned
#' 
#' @return an object of 'SCESet' class
#' 
#' @importFrom M3Drop M3DropDifferentialExpression
#' @importFrom scater get_exprs
#' @importFrom SC3 get_processed_dataset
#' 
#' @export
fsc3_get_hyperplanes.SCESet <- function(object, n_genes = 1000) {
    data <- SC3::get_processed_dataset(object)
    if(object@sc3$take_log) {
        data <- 2^(data) - 1
    }
    if (is.null(data)) {
        warning(paste0("Please run sc3_prepare() first!"))
        return(object)
    }
    tmp <- M3Drop::M3DropDifferentialExpression(data, "fdr", 2, suppress.plot = T)
    tmp <- tmp[!is.na(tmp$p.value) & !is.na(tmp$q.value), ]
    tmp <- tmp[order(tmp$q.value), ]
    object@sc3$hyperplanes <- as.character(tmp[1:n_genes, ]$Gene)
    return(object)
}

#' @rdname fsc3_get_hyperplanes.SCESet
#' @aliases fsc3_get_hyperplanes
#' @importClassesFrom scater SCESet
#' @export
setMethod("fsc3_get_hyperplanes", signature(object = "SCESet"), function(object, n_genes = 1000) {
    fsc3_get_hyperplanes.SCESet(object, n_genes)
})

#' Create a binary signature for each cell
#' 
#' 
#' @param object an object of 'SCESet' class
#' 
#' @importFrom SC3 get_processed_dataset
#' 
#' @return an object of 'SCESet' class
#' 
#' @export
fsc3_get_signatures.SCESet <- function(object) {
    hyperps <- object@sc3$hyperplanes
    if (is.null(hyperps)) {
        warning(paste0("Please run fsc3_get_hyperplanes() first!"))
        return(object)
    }
    data <- SC3::get_processed_dataset(object)
    data <- data[rownames(data) %in% hyperps, ]
    object@sc3$signatures <- signature_mapper(data)
    names(object@sc3$signatures) <- colnames(data)
    return(object)
}

#' @rdname fsc3_get_signatures.SCESet
#' @aliases fsc3_get_signatures
#' @importClassesFrom scater SCESet
#' @export
setMethod("fsc3_get_signatures", signature(object = "SCESet"), function(object) {
    fsc3_get_signatures.SCESet(object)
})

#' Define buckets of cells based on their signatures
#' 
#' Two cells are put in the same buckets if their binary signatures have a 
#' specific number of common bits.
#' 
#' @param object an object of 'SCESet' class
#' @param common_bits number of common bits which is enough to put two signatures
#' in the same bucket
#' 
#' @importFrom scater pData<-
#' @importFrom methods new
#' 
#' @return an object of 'SCESet' class
#' 
#' @export
fsc3_get_buckets.SCESet <- function(object, common_bits = NULL) {
    sigs <- object@sc3$signatures
    if (is.null(sigs)) {
        warning(paste0("Please run fsc3_get_signatures() first!"))
        return(object)
    }
    if (is.null(common_bits)) {
        warning(paste0("Please provide the 'common_bits' argument!"))
        return(object)
    }
    if (nchar(sigs[1]) < common_bits) {
        warning(paste0("Number of 'common_bits' must be less than the signature length!"))
        return(object)
    }
    buckets <- get_buckets(sigs, common_bits)
    message(paste0("Number of buckets is ", length(unique(buckets)), ". On average each bucket contains ", round(length(sigs) / length(unique(buckets))), " cells."))
    p_data <- object@phenoData@data
    p_data$fsc3_buckets <- buckets
    pData(object) <- new("AnnotatedDataFrame", data = p_data)
    object@sc3$buckets <- buckets
    return(object)
}

#' @rdname fsc3_get_buckets.SCESet
#' @aliases fsc3_get_buckets
#' @importClassesFrom scater SCESet
#' @export
setMethod("fsc3_get_buckets", signature(object = "SCESet"), function(object, common_bits = NULL) {
    fsc3_get_buckets.SCESet(object, common_bits)
})

#' Reduce the original dataset by using cell buckets
#' 
#' @param object an object of 'SCESet' class
#' @param d.region.min defines the minimum number of eigenvectors used for 
#' kmeans clustering as a fraction of the total number of cells. Default is 0.04.
#' @param d.region.max defines the maximum number of eigenvectors used for 
#' kmeans clustering as a fraction of the total number of cells. Default is 0.07.
#' 
#' @return an object of 'SCESet' class
#' 
#' @export
fsc3_select_repr_cells.SCESet <- function(object, d.region.min = 0.04, d.region.max = 0.07) {
    buckets <- object@sc3$buckets
    if (is.null(buckets)) {
        warning(paste0("Please run fsc3_get_buckets() first!"))
        return(object)
    }
    
    inds <- NULL
    for(i in 1:length(unique(buckets))) {
        if(length(which(buckets == i)) != 1) {
            inds <- c(inds, sample(which(buckets == i), 1))
        } else {
            inds <- c(inds, which(buckets == i))
        }
    }
    # the concept is similar to SVM, when clustering is performed on a subset of
    # cells, therefore the 'svm_train_inds' slot is used to store the indices of the
    # cells selected for clustering
    object@sc3$svm_train_inds <- inds
    # we also need to update number of dimensions used for clustering
    object@sc3$n_dim <- floor(d.region.min * length(inds)):ceiling(d.region.max * length(inds))
    # also update kmeans_nstart
    object@sc3$kmeans_nstart <- 50
    return(object)
}

#' @rdname fsc3_select_repr_cells.SCESet
#' @aliases fsc3_select_repr_cells
#' @importClassesFrom scater SCESet
#' @export
setMethod("fsc3_select_repr_cells", signature(object = "SCESet"), function(object, d.region.min = 0.04, d.region.max = 0.07) {
    fsc3_select_repr_cells.SCESet(object, d.region.min, d.region.max)
})

#' Reindex cell buckets based on the merged buckets results
#' 
#' @param object an object of 'SCESet' class
#' @param k number of clusters
#' 
#' @return an object of 'SCESet' class
#' 
#' @export
fsc3_get_clusters.SCESet <- function(object, k) {
    if (is.null(object@sc3$svm_train_inds)) {
        warning(paste0("Please run fsc3_select_repr_cells() first!"))
        return(object)
    }
    
    hc <- object@sc3$consensus[[as.character(k)]]$hc
    bucket_clusts <- SC3:::get_clusts(hc, k)
    
    names(bucket_clusts) <- object@sc3$svm_train_inds
    buckets <- object@sc3$buckets
    
    res <- rep(0, length(buckets))
    for(i in unique(sort(bucket_clusts))) {
        to_merge <- buckets[as.numeric(names(bucket_clusts[bucket_clusts == i]))]
        res[!is.na(match(buckets, to_merge))] <- i
    }
    
    object@sc3$svm_result <- res
    return(object)
}

#' @rdname fsc3_get_clusters.SCESet
#' @aliases fsc3_get_clusters
#' @importClassesFrom scater SCESet
#' @export
setMethod("fsc3_get_clusters", signature(object = "SCESet"), function(object, k) {
    fsc3_get_clusters.SCESet(object, k)
})

