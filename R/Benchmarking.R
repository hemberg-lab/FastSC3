#' Run support vector machines (\code{SVM}) prediction
#'
#' Train an \code{SVM} classifier on a training dataset (\code{train}) and then
#' classify a study dataset (\code{study}) using the classifier.
#'
#' @param train training data.frame with colnames, corresponding to training labels
#' @param study study data.frame
#' @param kern kernel to be used with SVM
#' @return classification of the study dataset
#'
#' @importFrom e1071 svm
#' @importFrom stats predict
#' @export
support_vector_machines <- function(train, study, kern = "linear") {
    train <- t(train)
    labs <- factor(rownames(train))
    rownames(train) <- NULL
    model <- tryCatch(e1071::svm(train, labs, kernel = kern), error = function(cond) return(NA))
    pred <- stats::predict(model, t(study))
    return(pred = pred)
}

#' Run random forest prediction
#'
#' Create a random forest calssifier based on a training dataset (\code{train}) and then
#' classify a study dataset (\code{study}) using the classifier.
#'
#' @param train training data.frame with colnames, corresponding to training labels
#' @param study study data.frame
#' @param ntree number of trees to be used in random forest
#' @return classification of the study dataset
#' @importFrom randomForest randomForest
#' @importFrom stats predict
#' @export
random_forest <- function(train, study, ntree = 50) {
    train <- t(train)
    train <- as.data.frame(train)
    y <- as.factor(rownames(train))
    study <- t(study)
    study <- as.data.frame(study)
    rownames(train) <- NULL
    rownames(study) <- NULL
    train_rf <- randomForest::randomForest(x = train, y = y, ntree = ntree)
    Prediction <- stats::predict(train_rf, study)
    return(Prediction)
}

#' Calculate median feature expression in every bucket
#' 
#' When the cell buckets are defined, it is possible to calculate a median
#' expression based on the expressions of the cells contained in a given
#' bucket.
#' 
#' @param object an object of 'SCESet' class
#' 
#' @importFrom scater pData<-
#' @importFrom dplyr group_by summarise %>%
#' @importFrom reshape2 melt dcast
#' 
#' @return an object of 'SCESet' class
#' 
#' @export
fsc3_get_buckets_med_exprs <- function(object) {
    buckets <- object@phenoData@data$fsc3_buckets
    if (is.null(buckets)) {
        warning(paste0("Please run fsc3_get_buckets() first!"))
        return(object)
    }
    dat <- exprs(object)[fData(object)$fsc3_features, ]
    colnames(dat) <- buckets
    dat <- reshape2::melt(dat)
    colnames(dat) <- c("gene", "bucket", "exprs")
    dat <- dat %>%
        group_by(gene, bucket) %>%
        summarise(
            med_exprs = median(exprs)
        )
    dat <- reshape2::dcast(dat, gene ~ bucket, value.var = "med_exprs")
    rownames(dat) <- dat$gene
    dat <- dat[ , 2:ncol(dat)]
    return(dat)
}

#' @importFrom scater pData<-
#' @importFrom proxy dist
#' @export
fsc3_map_by_similarity <- function(object_to_map, buckets_ref, similarity, scale, threshold) {
    dat <- exprs(object_to_map)[fData(object_to_map)$fsc3_features, ]
    if(scale) {
        dat <- scale(dat, center = TRUE, scale = TRUE)
        buckets_ref <- scale(buckets_ref, center = TRUE, scale = TRUE)
    }
    dat <- dat[order(rownames(dat)), ]
    buckets_ref <- buckets_ref[order(rownames(buckets_ref)), ]
    dat <- t(dat)
    buckets_ref <- t(buckets_ref)
    
    if(similarity == "euclidean") {
        res <- proxy::dist(buckets_ref, dat, method = "euclidean")
    }
    if(similarity == "cosine") {
        res <- proxy::dist(buckets_ref, dat, method = "cosine")
    }
    res <- matrix(res, ncol = nrow(buckets_ref), byrow = T)
    if(similarity == "euclidean") {
        res <- res / max(res, na.rm = T)
    }
    min_inds <- unlist(apply(-res, 1, nnet::which.is.max))
    mins <- unlist(apply(res, 1, min))
    buckets_assigned <- rownames(buckets_ref)[min_inds]
    buckets_assigned[mins > threshold] <- "unassigned"
    p_data <- object_to_map@phenoData@data
    p_data$fsc3_buckets_assigned <- buckets_assigned
    pData(object_to_map) <- new("AnnotatedDataFrame", data = p_data)
    return(object_to_map)
}

