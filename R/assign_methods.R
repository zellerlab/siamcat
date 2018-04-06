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