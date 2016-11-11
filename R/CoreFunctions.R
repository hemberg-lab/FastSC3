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

#' Plot Sankey diagram comparing two clusterings
#' 
#' Sometimes it is useful to see how the clusters in two different clustering
#' solutions correspond to each other. Sankey diagram is a good way to visualize
#' them. This function takes as input two clustering solutions and visualizes them
#' using a Sankey diagram. The order of the reference clusters is defined by their
#' labels in increasing order.
#' 
#' @param reference reference clustering labels
#' @param clusters clustering labels under investigations
#' @param plot_width width of the output plot in pixels
#' @param colors colors of the links between two clusterings. If defined please
#' note that each cluster in the reference clustering has to have its own color.
#' This should be a normal text vector, e.g. c("#FF0000", "#FFA500", "#008000")
#' 
#' @importFrom dplyr %>% summarise group_by
#' @importFrom reshape2 melt
#' @importFrom googleVis gvisSankey
plot_sankey <- function(reference, clusters, plot_width = 400, colors = "#808080") {
    res.all <- NULL
    for(j in unique(reference)){
        res <- NULL
        for(i in names(table(clusters))) {
            tmp <- length(intersect(which(clusters == i), which(reference == j)))
            res <- c(res, tmp)
        }
        res.all <- rbind(res.all, res)
    }
    colnames(res.all) <- names(table(clusters))
    rownames(res.all) <- unique(reference)
    
    res.all <- res.all[order(as.numeric(rownames(res.all))), 
                       order(as.numeric(table(clusters)), decreasing = T)]
    
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
    
    # remove cycles from the data
    res[, 1] <- paste0(res[, 1], " ")
    res[, 2] <- paste0(" ", res[, 2])
    
    colnames(res) <- c("From", "To", "Weight")
    
    colors <- paste(colors, collapse = "', '")
    colors <- paste0("['", colors, "']")

    Sankey <- gvisSankey(
        res,
        from="From",
        to="To",
        weight="Weight",
        options = list(
            width = plot_width, 
            sankey = paste0("{
                node:{
                    label:{
                        fontName:'Arial',
                        fontSize:11,color:
                        '#000000',
                        bold:true,
                        italic:false
                    },
                    colors:'#FFFFFF',
                    nodePadding:12
                },
                link:{
                    colorMode: 'source',
                    colors: ", colors, "
                },
                iterations:0
            }"
        ))
    )
    
    plot(Sankey)
}
