#' Fast-Johnson-Lindenstrauss-Transform (FJLT)
#'
#' This function calculates the FJLT and project onto d dimensions. The FJLT 
#' is faster than standard random pro-jections and just as easy to implement. 
#' It is based upon the preconditioning of a sparse projection matrix with a 
#' randomized Fourier transform.
#' 
#' Ailon, N. and Chazelle, B. Approximate nearest neighbors and the fast 
#' Johnson-Lindenstrauss transform. in Proceedings of the thirty-eighth 
#' annual ACM symposium on Theory of computing 557–563 (ACM, 2006).
#' 
#' Functions adapted from: \url{http://www.cs.ubc.ca/~jaquesn/MachineLearningTheory.pdf}
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
#' @importFrom stats cov
ftransformation <- function(dists, method) {
    if (method == "pca") {
        covar <- cov(scale(dists, center = TRUE, scale = TRUE))
        return(covar)
    } else if (method == "laplacian") {
        return(norm_laplacian(dists))
    }
}

#' @importFrom dplyr %>% summarise group_by
#' @importFrom reshape2 melt
create_data_for_sankey <- function(labs.ref, consens) {
    res.all <- NULL
    for(j in unique(labs.ref)){
        res <- NULL
        for(i in names(table(consens))) {
            tmp <- length(intersect(which(consens == i), which(labs.ref == j)))#/length(which(labs.ref == j))
            res <- c(res, tmp)
        }
        res.all <- rbind(res.all, res)
    }
    # res.all <- res.all*100
    colnames(res.all) <- names(table(consens))
    rownames(res.all) <- unique(labs.ref)
    
    res.all <- res.all[order(as.numeric(rownames(res.all))), order(as.numeric(table(consens)), decreasing = T)]
    # res.all <- log10(res.all + 1)
    res <- reshape2::melt(res.all)
    res <- res[res$value != 0, ]
    
    maxs <- res %>%
        dplyr::group_by(Var1) %>%
        dplyr::summarise(max = max(value))
    
    res <- merge(res, maxs)
    maxs <- res[res$value == res$max, ]
    res <- res[res$value != res$max, ]
    
    res <- rbind(maxs, res)
    res <- res[,1:3]
    
    # res <- res[,c(2,1,3)]
    mat<-function(x){
        return(paste0('["', x[1], ' "," ', x[2], '",', x[3], '],' ))
    }
    res.js <- apply(res,1,mat)
    fileConn<-file("~/Desktop/sankey.txt")
    writeLines(res.js, fileConn)
    close(fileConn)
}
