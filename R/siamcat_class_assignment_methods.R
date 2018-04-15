#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

###############################################################################
#' Assign a new phyloseq object to \code{x}
#'
#' @usage physeq(x) <- value
#'
#' @param x an object of class \link{siamcat-class}
#' @param value an object of class \link[phyloseq]{phyloseq-class}
#' @export
#' @docType methods
#' @rdname assign-physeq
#' @aliases assign-physeq
#'
#' @examples
#' # data(siamcat_example)
setGeneric("physeq<-", function(x, value) standardGeneric("physeq<-"))
#' @rdname assign-physeq
#' @aliases physeq<-
setMethod("physeq<-", c("siamcat","phyloseq"), function(x, value){
  siamcat(value, x@model_list, x@eval_data, x@label, x@norm_param, 
    x@data_split, x@pred_matrix)
})
#' @rdname assign-physeq
#' @aliases physeq<-
setMethod("physeq<-", c("siamcat","otu_table"), function(x, value){
  phyloseq <- physeq(x)
  otu_table(phyloseq) <- value
  siamcat(phyloseq, x@model_list, x@eval_data, x@label, x@norm_param, 
    x@data_split, x@pred_matrix)
})


###############################################################################
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


###############################################################################
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

###############################################################################
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

###############################################################################
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

###############################################################################
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

###############################################################################
#' Assign a new data_split object to \code{x}
#'
#' @usage data_split(x) <- value
#'
#' @param x an object of class \link{siamcat-class}
#' @param value an object of class \link{data_split-class}
#' @export
#' @docType methods
#' @rdname assign-data_split
#' @aliases assign-data_split
#'
#' @examples
#' # data(siamcat_example)
setGeneric("data_split<-", function(x, value) standardGeneric("data_split<-"))
#' @rdname assign-data_split
#' @aliases data_split<-
setMethod("data_split<-", c("siamcat","data_split"), function(x, value){
  siamcat(value, x@label, x@model_list, x@phyloseq, x@orig_feat,
          x@eval_data, x@pred_matrix, x@norm_param)
})

###############################################################################
#' Assign a new otu_table object to \code{x} orig_feat slot
#'
#' @usage orig_feat(x) <- value
#'
#' @param x an object of class \link{siamcat-class}
#' @param value an object of class \link[phyloseq]{otu_table-class}
#' @export
#' @docType methods
#' @rdname assign-orig_feat
#' @aliases assign-orig_feat
#'
#' @examples
#' # data(siamcat_example)
setGeneric("orig_feat<-", function(x, value) standardGeneric("orig_feat<-"))
#' @rdname assign-orig_feat
#' @aliases orig_feat<-
setMethod("orig_feat<-", c("siamcat","otu_table"), function(x, value){
  siamcat(value, x@label, x@model_list, x@phyloseq, x@data_split,
          x@eval_data, x@pred_matrix, x@norm_param)
})

###############################################################################
#' Assign a new otu_table object to \code{x} features slot
#'
#' @usage features(x) <- value
#'
#' @param x an object of class \link{siamcat-class}
#' @param value an object of class \link[phyloseq]{otu_table-class}
#' @export
#' @docType methods
#' @rdname assign-features
#' @aliases assign-features
#'
#' @examples
#' # data(siamcat_example)
setGeneric("features<-", function(x, value) standardGeneric("features<-"))
#' @rdname assign-features
#' @aliases features<-
setMethod("features<-", c("siamcat","otu_table"), function(x, value){
  otu_table(physeq(x)) <-  value
  return(x)
})

###############################################################################
#' Assign a new sam_data object to \code{x}
#'
#' @usage meta(x) <- value
#'
#' @param x an object of class \link{siamcat-class}
#' @param value an object of class \link[phyloseq]{sample_data-class}
#' @export
#' @docType methods
#' @rdname assign-meta
#' @aliases assign-meta
#'
#' @examples
#' # data(siamcat_example)
setGeneric("meta<-", function(x, value) standardGeneric("meta<-"))
#' @rdname assign-meta
#' @aliases meta<-
setMethod("meta<-", c("siamcat","sample_data"), function(x, value){
  sample_data(physeq(x)) <-  value
  return(x)
})
