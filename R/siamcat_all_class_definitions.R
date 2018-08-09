#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

# ##############################################################################
# label object

# check label object for validity
#'@keywords internal
check.label <- function(object){
    errors <- character()
    # check that all entries are there
    if (!all(names(object) == c('label', 'info', 'type'))){
        msg <- 'Label object does not contain all needed entries!'
        errors <- c(errors, msg)
    }
    # check that label is binary or test
    if (object$type != 'BINARY' & object$type != 'TEST'){
        msg <- 'Label object is neither binary nor a test label!'
        errors <- c(errors, msg)
    }
    # check that info and label match up
    if (!all(sort(unique(object$label)) == object$info)){
        msg <- 'label info does not match to label entries!'
        errors <- c(errors, msg)
    }
    # check that info has names
    if (is.null(names(object$info))){
        msg <- 'Label info does not contain group names!'
        errors <- c(errors, msg)
    }
    if (length(errors) == 0) TRUE else errors
}

#' The S4 class for storing label info.
#' @name label-class
#' @rdname label-class
#' @slot .Data inherited from \code{\link{list}} class, contains
#' a list with:
#' \itemize{
#'     \item \code{label} numeric vector, specifying to which category samples
#'     belong, usualy 1 and -1
#'     \item \code{type} contains information about the label type
#'     \item \code{info} list with additional informations about the dataset
#' }
#' @exportClass label
setClass("label",
         contains = "list",
         validity=check.label)

# ##############################################################################
# Filtered data/parameters

# check filtered data/parameters for validity
#'@keywords internal
check.filt.feat <- function(object){
    errors <- character()
    # TODO
    if (length(errors) == 0) TRUE else errors
}

#' The S4 class for storing the filter features/paramters
#' @name filt_feat-class
#' @rdname filt_feat-class
#' @slot filt.feat An object of class \link[phyloseq]{otu_table-class} storing
#' the filtered features
#' @slot filt.param A list storing the parameters of the feature filtering
#' @exportClass filt_feat
setClass("filt_feat",
         representation(filt.param= "list",
                        filt.feat="otu_table"),
         validity=check.filt.feat,
         prototype(filt.param = list(),
                   filt.feat = otu_table(NA, taxa_are_rows=TRUE,
                                         errorIfNULL=FALSE)))

# ##############################################################################
# Associations

# check associations for validity
#'@keywords internal
check.assoc <- function(object){
    errors <- character()
    # TODO
    if (length(errors) == 0) TRUE else errors
}

#' The S4 class for storing the results of the association testing
#' @name associations-class
#' @rdname associations-class
#' @slot assoc.results a data.frame containing the results of the association
#'  testing
#' @slot assoc.param a list containing the parameters for the association
#' testing
#' @exportClass associations
setClass("associations",
         representation(assoc.results='data.frame',
                        assoc.param = "list"),
        validity=check.assoc)

# ##############################################################################
# Normalization data/parameters

# check normalization data/parameters for validity
#'@keywords internal
check.norm.feat <- function(object){
    errors <- character()
    # TODO
    if (length(errors) == 0) TRUE else errors
}

#' The S4 class for storing the normalization data/parameters
#' @name norm_feat-class
#' @rdname norm_feat-class
#' @slot norm.feat An object of class \link[phyloseq]{otu_table-class} storing
#' the normalized features
#' @slot norm.param A list with:
#'     \itemize{
#'     \item \code{norm.method} the normalization method used
#'     \item \code{retained.feat} the names of features retained after filtering
#'     \item \code{log.n0} pseudocount
#'     \item \code{n.p} vector norm
#'     \item \code{norm.margin} margin for the normalization
#' } and additional entries depending on the normalization method used.
#' @exportClass norm_feat
setClass("norm_feat",
         representation(norm.param= "list",
                        norm.feat='otu_table'),
         validity=check.norm.feat,
         prototype(norm.param = list(),
                   norm.feat = otu_table(NA, taxa_are_rows=TRUE,
                                         errorIfNULL=FALSE)))

# ##############################################################################
# Data split

# check data split for validity
#'@keywords internal
check.data.split <- function(object){
    errors <- character()
    # check names
    if (!all(names(object) == c('training.folds', 'test.folds',
                                'num.resample', 'num.folds'))){
        msg <- 'Data split does not contain all needed entries!'
        errors <- c(errors, msg)
    }
    # check that num.resample and num.folds are numbers
    if (length(object$num.resample) != 1 |
        class(object$num.resample) != 'numeric'){
        msg <- 'num.resample should be numeric and of length 1!'
        errors <- c(errors, msg)
    }
    if (length(object$num.folds) != 1 |
        class(object$num.folds) != 'numeric'){
        msg <- 'num.folds should be numeric and of length 1!'
        errors <- c(errors, msg)
    }
    # check that training.folds is a list (of the right length)
    if (class(object$training.folds) != 'list'){
        msg <- 'training.folds should be a list!'
        errors <- c(errors, msg)
    }
    if (length(object$training.folds) != object$num.resample){
        msg <- 'training.folds should be of length num.resample!'
        errors <- c(errors, msg)
    }
    if (!all(vapply(object$training.folds, length,
             FUN.VALUE = numeric(1)) == object$num.folds)){
        msg <- 'All training.folds should be of length num.folds!'
        errors <- c(errors, msg)
    }
    # same for test.folds
    if (class(object$test.folds) != 'list'){
        msg <- 'test.folds should be a list!'
        errors <- c(errors, msg)
    }
    if (length(object$test.folds) != object$num.resample){
        msg <- 'test.folds should be of length num.resample!'
        errors <- c(errors, msg)
    }
    if (!all(vapply(object$test.folds, length,
             FUN.VALUE = numeric(1)) == object$num.folds)){
        msg <- 'All test.folds should be of length num.folds!'
        errors <- c(errors, msg)
    }
    # check that no samples are in training and test fold at the same time
    test.consistency <- any(c(vapply(seq_len(object$num.resample),
        FUN=function(x){vapply(seq_len(object$num.folds),
            FUN=function(y){
                any(object$test.folds[[x]][[y]] %in%
                    object$training.folds[[x]][[y]])
                }, FUN.VALUE = logical(1))
            }, FUN.VALUE=logical(object$num.folds))))
    if (test.consistency){
        msg <- 'Some samples are in both training and test folds!'
        errors <- c(errors, msg)
    }


    if (length(errors) == 0) TRUE else errors
}

#' The S4 class for storing data splits
#' @name data_split-class
#' @rdname data_split-class
#' @slot .Data inherited from \code{\link{list}} class, contains
#' a list with:
#'     \itemize{
#'     \item \code{training.folds} a list - for each cv fold contains ids of
#' samples used for training
#'     \item \code{test.folds} a list - for each cv fold contains ids of
#' samples used for testing
#'     \item \code{num.resample} number of repetition rounds for cv
#'     \item \code{num.folds} number of folds for cv}
#' @exportClass data_split
setClass("data_split",
         contains = "list",
         validity=check.data.split)

# ##############################################################################
# Model list

# check model list for validity
#'@keywords internal
check.model.list <- function(object){
    errors <- character()
    # check model.type
    length.model.type <- length(object@model.type)
    if (length.model.type != 1){
        msg <- paste0('Model type is of length ', length.model.type, '.',
            ' Should be length 1!')
        errors <- c(errors, msg)
    }
    # check that all models in the list are mlr models
    if (any(vapply(object@models, class,
            FUN.VALUE = character(1)) != 'WrappedModel')){
        msg <- 'Models are supposed to be mlr-WrappedModels!'
        errors <- c(errors, msg)
    }
    if (length(errors) == 0) TRUE else errors
}

#' The S4 class for storing models.
#' @name model_list-class
#' @rdname model_list-class
#' @slot models a list with models obtained from \link{train.model}
#' @slot model.type name of the method used by \link{train.model}
#' @exportClass model_list
setClass("model_list",
         representation(models = "list",
                        model.type = "character"),
         validity=check.model.list)


# ##############################################################################
# Prediction matrix

#' The S4 class for storing predictions.
#' @name pred_matrix-class
#' @rdname pred_matrix-class
#' @slot .Data inherited from \code{\link{matrix}} class, contains
#' a matrix with predictions made by \link{make.predictions} function
#' @exportClass pred_matrix
setClass("pred_matrix", contains = "matrix")

# ##############################################################################
# Evaluation data

# check evaluation data for validity
#'@keywords internal
check.eval.data <- function(object){
    errors <- character()
    # TODO
    if (length(errors) == 0) TRUE else errors
}

#' The S4 class for storing evaluation data.
#' @name eval_data-class
#' @rdname eval_data-class
#' @slot .Data inherited from \code{\link{list}} class, contains
#' a list with:
#'     \itemize{
#'     \item \code{$roc} average ROC-curve across repeats or a single
#'     ROC-curve on the complete dataset (object of class \link[pROC]{roc});
#'     \item \code{$auroc} AUC value for the average ROC-curve;
#'     \item \code{$prc} average Precision Recall curve across repeats or
#'     a single PR-curve on the complete dataset;
#'     \item \code{$auprc} AUC value for the average PR-curve;
#'     \item \code{$ev} list containing for different decision thresholds the
#'     number of false positives, false negatives, true negatives, and
#'     true positives;
#' }. If \code{prediction} had more than one column, i.e. if the models has
#'     been trained with several repeats, the function will additonally
#'     return \itemize{
#'     \item \code{$roc.all} list of roc objects (see \link[pROC]{roc}) for
#'     every repeat;
#'     \item \code{$auroc.all} vector of AUC values for the ROC curves for
#'     every repeat;
#'     \item \code{$prc.all} list of PR curves for every repeat;
#'     \item \code{$auprc.all} vector of AUC values for the PR curves for every
#'     repeat;
#'     \item \code{$ev.all} list of false positive, false negatives, true
#'     negatives, true positives, and thresholds for the different repeats.
#'}
#' @exportClass eval_data
setClass("eval_data", contains = "list", validity=check.eval.data)

# ##############################################################################
# SIAMCAT object

#' The S4 class for storing taxa-abundance information and models.
#' @name siamcat-class
#' @rdname siamcat-class
#' @slot phyloseq object of class \link[phyloseq]{phyloseq-class}
#' @slot label an object of class \link{label-class}
#' @slot data_split an object of class \link{data_split-class}
#' @slot norm_param a list of normalzation parameters, see
#'     \link{normalize.features} for more details
#' @slot model_list an object of class \link{model_list-class}
#' @slot eval_data an object of class \link{eval_data-class}
#' @slot pred_matrix an object of class \link{pred_matrix-class}
#' @slot assoc an object of class \link{assoc-class}
#' @exportClass siamcat
setClass(
    "siamcat",
    representation(
        phyloseq = "phyloseq",
        label = "label",
        filt_feat = "filt_feat",
        associations = "associations",
        norm_feat = "norm_feat",
        data_split = "data_split",
        model_list = "model_list",
        pred_matrix = "matrix",
        eval_data = "list"
    )
)
