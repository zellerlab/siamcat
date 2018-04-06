#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' The S4 class for storing models.
#' @name model_list-class
#' @rdname model_list-class
#' @slot models a list with models obtained from \link{train.model}
#' @slot model.type name of the method used by \link{train.model}
#' @exportClass model_list
setClass("model_list", representation(models = "list", model.type = "character"))

#' The S4 class for storing data splits
#' @name data_split-class
#' @rdname data_split-class
#' @slot training.folds a list - for each cv fold contains ids of
#' samples used for training
#' @slot test.folds a list - for each cv fold contains ids of
#' samples used for testing
#' @slot num.resample number of repetition rounds for cv
#' @slot num.folds number of folds for cv
#' @exportClass data_split
setClass("data_split", representation(training.folds = "list", test.folds = "list", num.resample = "numeric", num.folds = "numeric"))

#' The S4 class for storing label info.
#' @name label-class
#' @rdname label-class
#' @slot label numeric vector, specifying to which category samples belong,
#' usualy made of 1s and -1s
#' @slot header contains information from the header of the label file
#' @slot info list with additional informations about the dataset
#' @slot positive.lab specifies which of two numbers in label is a positive label
#' @slot negative.lab specifies which of two numbers in label is a negative label
#' @slot n.idx numeric vector - on which positions in the label there are samples
#' withc negative label
#' @slot p.idx numeric vector - on which positions in the label there are samples
#' withc positive label
#' @slot n.lab character string with a name for the negative label (e.g. 'healthy')
#' @slot p.lab character string with a name for the positive label (e.g. 'cancer')
#' @exportClass label
setClass("label", representation(label = "vector", header = "character", info = "list", positive.lab = "numeric", negative.lab = "numeric", 
    n.idx = "vector", p.idx = "vector", n.lab = "character", p.lab = "character"))

#' The S4 class for storing taxa-abundance information and models.
#' @name siamcat-class
#' @rdname siamcat-class
#' @slot phyloseq object of class \link[phyloseq]{phyloseq-class}
#' @slot label an object of class \link{label-class}
#' @slot orig_feat an object of class \link[phyloseq]{otu_table-class}
#' @slot data_split an object of class \link{data_split-class}
#' @slot norm_param a list of normalzation parameters, see \link{normalize.features} for more details
#' @slot model_list an object of class \link{model_list-class}
#' @slot eval_data list containing \itemize{
#'  \item \code{$roc.average} average ROC-curve across repeats or a single ROC-curve on complete dataset;
#'  \item \code{$auc.average} AUC value for the average ROC-curve;
#'  \item \code{$ev.list} list of \code{length(num.folds)}, containing for different decision thresholds the number of false positives, false negatives, true negatives, and true positives;
#'  \item \code{$pr.list} list of \code{length(num.folds)}, containing the positive predictive value (precision) and true positive rate (recall) values used to plot the PR curves;
#' }. If \code{prediction} had more than one column, i.e. if the models has been trained with several repeats, the function will additonally return \itemize{
#'  \item \code{$roc.all} list of roc objects (see \link[pROC]{roc}) for every repeat;
#'  \item \code{$aucspr} vector of AUC values for the PR curves for every repeat;
#'  \item \code{$auc.all} vector of AUC values for the ROC curves for every repeat
#'}
#' @slot pred_matrix a matrix with predictions made by \link{make.predictions} function
#' @exportClass siamcat
setClass("siamcat", representation(model_list = "model_list", phyloseq = "phyloseq", orig_feat = "otu_table", eval_data = "list", 
    label = "label", norm_param = "list", data_split = "data_split", pred_matrix = "matrix"))
