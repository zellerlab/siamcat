#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' The S4 class for storing models.
#' @name model_list-class
#' @rdname model_list-class
#' @slot models a list with models obtained from \link{train.model}
#' @slot model.type name of the method used by \link{train.model}
#' @exportClass model_list
setClass("model_list",
    representation(models = "list",
        model.type = "character"))

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
setClass("data_split", contains = "list")

#' The S4 class for storing label info.
#' @name label-class
#' @rdname label-class
#' @slot .Data inherited from \code{\link{list}} class, contains
#' a list with:
#' \itemize{
#'     \item \code{label} numeric vector, specifying to which category samples
#'     belong, usualy made of 1s and -1s
#'     \item \code{header} contains information from the header of the label
#'     file
#'     \item \code{info} list with additional informations about the dataset
#'     \item \code{positive.lab} specifies which of two numbers in label is a
#'     positive label
#'     \item \code{negative.lab} specifies which of two numbers in label is a
#'     negative label
#'     \item \code{n.idx} numeric vector - on which positions in the label there
#'     are samples with negative label
#'     \item \code{p.idx} numeric vector - on which positions in the label there
#'     are samples with positive label
#'     \item \code{n.lab} character string with a name for the negative label
#'     (e.g. 'healthy')
#'     \item \code{p.lab} character string with a name for the positive label
#'     (e.g. 'cancer')}
#' @exportClass label
setClass("label", contains = "list")

#' The S4 class for storing label info.
#' @name pred_matrix-class
#' @rdname pred_matrix-class
#' @slot .Data inherited from \code{\link{matrix}} class, contains
#' a matrix with
#'     predictions made by \link{make.predictions} function
#' @exportClass pred_matrix
setClass("pred_matrix", contains = "matrix")

#' The S4 class for storing original features info.
#' @name orig_feat-class
#' @rdname orig_feat-class
#' @slot taxa_are_rows A single logical specifying the orientation
#' of the abundance table
#' @slot .Data inherited from \code{\link{matrix}} class, contains
#' a matrix with predictions made by \link{make.predictions} function
#' @exportClass orig_feat
setClass("orig_feat", contains = "otu_table")

#' The S4 class for storing evaluation data.
#' @name eval_data-class
#' @rdname eval_data-class
#' @slot .Data inherited from \code{\link{list}} class, contains
#' a list with:
#'     \itemize{
#'     \item \code{$roc.average} average ROC-curve across repeats or a single
#'     ROC-curve on complete dataset;
#'     \item \code{$auc.average} AUC value for the average ROC-curve;
#'     \item \code{$ev.list} list of \code{length(num.folds)}, containing for
#'     different decision thresholds the number of false positives, false
#'     negatives, true negatives, and true positives;
#'     \item \code{$pr.list} list of \code{length(num.folds)}, containing the
#'     positive predictive value (precision) and true positive rate (recall)
#'     values used to plot the PR curves;
#' }. If \code{prediction} had more than one column, i.e. if the models has
#'     been trained with several repeats, the function will additonally
#'     return \itemize{
#'     \item \code{$roc.all} list of roc objects (see \link[pROC]{roc}) for
#'     every repeat;
#'     \item \code{$aucspr} vector of AUC values for the PR curves for every
#'     repeat;
#'     \item \code{$auc.all} vector of AUC values for the ROC curves for every
#'     repeat
#'}
#' @exportClass eval_data
setClass("eval_data", contains = "list")

#' The S4 class for storing the normalization paramters
#' @name norm_param-class
#' @rdname norm_param-class
#' @slot .Data inherited from \code{\link{list}} class, contains
#' a list with:
#'     \itemize{
#'     \item \code{norm.method} the normalization method used
#'     \item \code{retained.feat} the names of features retained after filtering
#'     \item \code{log.n0} pseudocount
#'     \item \code{n.p} vector norm
#'     \item \code{norm.margin} margin for the normalization
#' } and additional entries depending on the normalization method used.
#' @exportClass norm_param
setClass("norm_param", contains = "list")

#' The S4 class for storing taxa-abundance information and models.
#' @name siamcat-class
#' @rdname siamcat-class
#' @slot phyloseq object of class \link[phyloseq]{phyloseq-class}
#' @slot label an object of class \link{label-class}
#' @slot orig_feat an object of class \link[phyloseq]{otu_table-class}
#' @slot data_split an object of class \link{data_split-class}
#' @slot norm_param a list of normalzation parameters, see
#'     \link{normalize.features} for more details
#' @slot model_list an object of class \link{model_list-class}
#' @slot eval_data an object of class \link{eval_data-class}
#' @slot pred_matrix an object of class \link{pred_matrix-class}
#' @exportClass siamcat
setClass(
    "siamcat",
    representation(
        model_list = "model_list",
        phyloseq = "phyloseq",
        orig_feat = "orig_feat",
        eval_data = "list",
        label = "label",
        norm_param = "list",
        data_split = "data_split",
        pred_matrix = "matrix"
    )
)
