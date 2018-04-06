#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between Microbial
### Communities And host phenoTypes
### EMBL Heidelberg 2012-2018 GNU GPL 3.0

################################################################################
#' Assign a new phyloseq object to \code{x}
#'
#' @usage phyloseq(x) <- value
#'
#' @param x an object of class \link{siamcat-class}
#' @param value an object of class \link[phyloseq]{phyloseq-class}
#' @export
#' @docType methods
#' @rdname assign-phyloseq
#' @aliases assign-phyloseq
#'
#' @examples
#' # data(siamcat_example)
setGeneric("phyloseq<-", function(x, value) standardGeneric("phyloseq<-"))
#' @rdname assign-phyloseq
#' @aliases phyloseq<-
setMethod("phyloseq<-", c("siamcat","phyloseq"), function(x, value){
  siamcat(value, x@model_list, x@eval_data, x@label, x@norm_param, x@data_split, x@pred_matrix)
})


################################################################################
#' Assign a new label object to \code{x}
#'
#' @usage label(x) <- value
#'
#' @param x an object of class \link{siamcat-class}
#' @param value an object of class \link{label-class}
#' @export
#' @docType methods
#' @rdname assign-label
#' @aliases assign-label
#'
#' @examples
#' # data(siamcat_example)
setGeneric("label<-", function(x, value) standardGeneric("label<-"))
#' @rdname assign-label
#' @aliases label<-
setMethod("label<-", c("siamcat","label"), function(x, value){
  siamcat(value, x@model_list, x@eval_data, x@phyloseq, x@orig_feat, 
          x@norm_param, x@data_split, x@pred_matrix)
})

################################################################################
#' Assign a new model_list object to \code{x}
#'
#' @usage model_list(x) <- value
#'
#' @param x an object of class \link{siamcat-class}
#' @param value an object of class \link{model_list-class}
#' @export
#' @docType methods
#' @rdname assign-model_list
#' @aliases assign-model_list
#'
#' @examples
#' # data(siamcat_example)
setGeneric("model_list<-", function(x, value) standardGeneric("model_list<-"))
#' @rdname assign-model_list
#' @aliases model_list<-
setMethod("model_list<-", c("siamcat","model_list"), function(x, value){
  siamcat(value, x@label, x@eval_data, x@phyloseq, x@orig_feat, 
          x@norm_param, x@data_split, x@pred_matrix)
})

################################################################################
#' Assign a new eval_data object to \code{x}
#'
#' @usage eval_data(x) <- value
#'
#' @param x an object of class \link{siamcat-class}
#' @param value an eval_data list
#' @export
#' @docType methods
#' @rdname assign-eval_data
#' @aliases assign-eval_data
#'
#' @examples
#' # data(siamcat_example)
setGeneric("eval_data<-", function(x, value) standardGeneric("eval_data<-"))
#' @rdname assign-eval_data
#' @aliases eval_data<-
setMethod("eval_data<-", c("siamcat","list"), function(x, value){
  siamcat(value, x@label, x@model_list, x@phyloseq, x@orig_feat, 
          x@norm_param, x@data_split, x@pred_matrix)
})

################################################################################
#' Assign a new norm_param object to \code{x}
#'
#' @usage norm_param(x) <- value
#'
#' @param x an object of class \link{siamcat-class}
#' @param value an norm_param list
#' @export
#' @docType methods
#' @rdname assign-norm_param
#' @aliases assign-norm_param
#'
#' @examples
#' # data(siamcat_example)
setGeneric("norm_param<-", function(x, value) standardGeneric("norm_param<-"))
#' @rdname assign-norm_param
#' @aliases norm_param<-
setMethod("norm_param<-", c("siamcat","list"), function(x, value){
  siamcat(value, x@label, x@model_list, x@phyloseq, x@orig_feat, 
          x@eval_data, x@data_split, x@pred_matrix)
})

################################################################################
#' Assign a new pred_matrix object to \code{x}
#'
#' @usage pred_matrix(x) <- value
#'
#' @param x an object of class \link{siamcat-class}
#' @param value an pred_matrix matrix
#' @export
#' @docType methods
#' @rdname assign-pred_matrix
#' @aliases assign-pred_matrix
#'
#' @examples
#' # data(siamcat_example)
setGeneric("pred_matrix<-", function(x, value) standardGeneric("pred_matrix<-"))
#' @rdname assign-pred_matrix
#' @aliases pred_matrix<-
setMethod("pred_matrix<-", c("siamcat","matrix"), function(x, value){
  siamcat(value, x@label, x@model_list, x@phyloseq, x@orig_feat, 
          x@eval_data, x@data_split, x@norm_param)
})