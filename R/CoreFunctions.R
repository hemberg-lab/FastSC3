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
#' 
#' @export
plot_sankey <- function(reference, clusters, plot_width = 400, plot_height = 600, colors = NULL) {
    res.all <- NULL
    for(j in names(table(reference))){
        res <- NULL
        for(i in names(table(clusters))) {
            tmp <- length(intersect(which(clusters == i), which(reference == j)))
            res <- c(res, tmp)
        }
        res.all <- rbind(res.all, res)
    }
    colnames(res.all) <- names(table(clusters))
    rownames(res.all) <- names(table(reference))
    
    res.all <- res.all[order(as.numeric(table(reference)), decreasing = T), 
                       order(as.numeric(table(clusters)), decreasing = T)]
    
    res <- reshape2::melt(res.all)
    res <- res[res$value != 0, ]
    
    maxs <- res %>%
        dplyr::group_by(Var1) %>%
        dplyr::summarise(max = max(value))
    
    res <- merge(res, maxs)
    maxs <- res[res$value == res$max, ]
    maxs <- maxs[order(maxs$value, decreasing = T), ]
    res <- res[res$value != res$max, ]
    res <- rbind(maxs, res)
    res <- res[,1:3]
    
    # remove cycles from the data
    res[, 1] <- paste0(res[, 1], " ")
    res[, 2] <- paste0(" ", res[, 2])
    
    colnames(res) <- c("From", "To", "Weight")
    
    if(!is.null(colors)) {
        colors <- paste(colors, collapse = "', '")
        colors <- paste0("['", colors, "']")
    }

    Sankey <- gvisSankey(
        res,
        from="From",
        to="To",
        weight="Weight",
        options = list(
            width = plot_width,
            height = plot_height,
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
                },", if(!is.null(colors)) {
                paste0("link:{
                    colorMode: 'source',
                    colors: ", colors, "
                },")},
                "iterations:0
            }"
        ))
    )
    
    plot(Sankey)
}



get_common_fsc3_features <- function(objects) {
    common_features <- NULL
    for(object in objects) {
        f_data <- object@featureData@data
        if(is.null(f_data$feature_symbol)) {
            message("Every input object has to contain 'feature_symbol' column in the featureData slot! Please check your input objects!")
            return(NULL)
        }
        object_features <- as.character(f_data$feature_symbol)[f_data$fsc3_features]
        common_features <- c(common_features, object_features)
    }
    common_features <- unique(common_features)
    common_inds <- rep(TRUE, length(common_features))
    for(object in objects) {
        f_data <- object@featureData@data
        inds <- match(common_features, f_data$feature_symbol)
        common_inds[is.na(inds)] <- FALSE
    }
    return(common_features[common_inds])
}

fsc3_assign_signatures <- function(object_to_assign, object_ref, threshold = 0.7) {
    # reference object
    buckets_ref <- unique(pData(object_ref)[, c("fsc3_buckets", "fsc3_buckets_signatures")])
    # convert factors to strings
    if(is.factor(buckets_ref$fsc3_buckets)) {
        buckets_ref$fsc3_buckets <- levels(buckets_ref$fsc3_buckets)[buckets_ref$fsc3_buckets]
    }
    # object to assign
    sigs_to_assign <- pData(object_to_assign)$fsc3_signatures
    
    buckets_assigned <- NULL
    for(sig in sigs_to_assign) {
        sig <- unlist(strsplit(sig, ""))
        common_bits <- c()
        for(sig_ref in buckets_ref$fsc3_buckets_signatures) {
            tmp <- unlist(strsplit(sig_ref, ""))
            common_bits <- c(common_bits, length(which(tmp == sig))/length(which(tmp != "_")))
        }
        if(max(common_bits) >= threshold) {
            tmp <- nnet::which.is.max(common_bits)
            buckets_assigned <- c(buckets_assigned, buckets_ref$fsc3_buckets[tmp])
        } else {
            buckets_assigned <- c(buckets_assigned, "unassigned")
        }
    }
    
    p_data <- object_to_assign@phenoData@data
    p_data$fsc3_buckets_assigned <- buckets_assigned
    pData(object_to_assign) <- new("AnnotatedDataFrame", data = p_data)

    return(object_to_assign)
}

fsc3_plot_sig_expression <- function(object, fColumn = NULL, n_cells = 10, hc = NULL) {
    if(is.null(fColumn)) {
        warning("Please provide the fColumn argument!")
        return(object)
    }
    p_data <- object@phenoData@data
    sigs <- p_data$fsc3_signatures[order(p_data[[fColumn]])]
    mat <- do.call(cbind, lapply(strsplit(sigs, ""), as.numeric))
    colnames(mat) <- p_data[[fColumn]][order(p_data[[fColumn]])]
    to_plot <- NULL
    gaps <- NULL
    gaps_it <- 0
    for(i in unique(colnames(mat))) {
        tmp <- mat[,colnames(mat) == i]
        if(ncol(tmp) > n_cells) {
            tmp <- tmp[, sample(1:ncol(tmp), n_cells)]
        }
        colnames(tmp) <- rep(i, ncol(tmp))
        to_plot <- cbind(to_plot, tmp)
        gaps <- c(gaps, gaps_it + ncol(tmp))
        gaps_it <- gaps_it + ncol(tmp)
    }
    if(is.null(hc)) {
        phm <- pheatmap::pheatmap(to_plot, cluster_cols = F, gaps_col = gaps)
    } else {
        phm <- pheatmap::pheatmap(to_plot, cluster_cols = F, cluster_rows = hc, gaps_col = gaps)
    }
}

