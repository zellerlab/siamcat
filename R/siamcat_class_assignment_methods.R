#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

################################################################################
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
#' @return none
#'
#' @examples
#' data(siamcat_example)
#' physeq(siamcat_example) <- physeq(siamcat_example)
setGeneric("physeq<-", function(x, value)
    standardGeneric("physeq<-"))
#' @rdname assign-physeq
#' @aliases physeq<-
setReplaceMethod("physeq", c("siamcat", "phyloseq"), function(x, value) {
    x@phyloseq <- value
    validObject(x)
    return(x)
})

################################################################################
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
#' @return none
#'
#' @examples
#' data(siamcat_example)
#' orig_feat(siamcat_example) <- orig_feat(siamcat_example)
setGeneric("orig_feat<-", function(x, value)
    standardGeneric("orig_feat<-"))
#' @rdname assign-orig_feat
#' @aliases orig_feat<-
setMethod("orig_feat<-", c("siamcat", "otu_table"), function(x, value) {
    temp <- physeq(x)
    args.list <- list('otu_table'=value, 'sam_data'=meta(x),
        'tax_table'=tax_table(temp, errorIfNULL=FALSE),
        'phy_tree'=phy_tree(temp, errorIfNULL=FALSE))
    physeq.new <- do.call("new", c(list(Class = "phyloseq"), args.list))
    physeq(x) <- physeq.new
    return(x)
})

################################################################################
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
#' @return none
#'
#' @examples
#' data(siamcat_example)
#' meta(siamcat_example) <- meta(siamcat_example)
setGeneric("meta<-", function(x, value)
    standardGeneric("meta<-"))
#' @rdname assign-meta
#' @aliases meta<-
setMethod("meta<-", c("siamcat", "sample_data"), function(x, value) {
    sample_data(physeq(x)) <-  value
    return(x)
})

################################################################################
#' @title Assign a new label object to a SIAMCAT object
#'
#' @usage label(x) <- value
#'
#' @param x an object of class \link{siamcat-class}
#' @param value an list (in label format)
#' @export
#' @docType methods
#' @rdname assign-label
#' @aliases assign-label
#' @return none
#'
#' @examples
#' data(siamcat_example)
#' label(siamcat_example) <- label(siamcat_example)
setGeneric("label<-", function(x, value)
    standardGeneric("label<-"))
#' @rdname assign-label
#' @aliases label<-
setReplaceMethod("label", c("siamcat", "list"), function(x, value) {
    x@label <- value
    validObject(x)
    return(x)
})

################################################################################
#' Assign a new filt_feat object to \code{x}
#'
#' @usage filt_feat(x) <- value
#'
#' @param x an object of class \link{siamcat-class}
#' @param value an filt_feat object
#' @export
#' @docType methods
#' @rdname assign-filt_feat
#' @aliases assign-filt_feat
#' @return none
#'
#' @examples
#' data(siamcat_example)
#' filt_feat(siamcat_example) <- list(
#'     filt.feat=filt_feat(siamcat_example),
#'     filt.param=filt_params(siamcat_example))
setGeneric("filt_feat<-", function(x, value)
    standardGeneric("filt_feat<-"))
#' @rdname assign-filt_feat
#' @aliases filt_feat<-
setReplaceMethod("filt_feat", c("siamcat", "list"), function(x, value) {
    x@filt_feat <- value
    validObject(x)
    return(x)
})


################################################################################
#' Assign a new assocications object to \code{x}
#'
#' @usage associations(x) <- value
#'
#' @param x an object of class \link{siamcat-class}
#' @param value an associations object
#' @export
#' @docType methods
#' @rdname assign-associations
#' @aliases assign-associations
#' @return none
#'
#' @examples
#' data(siamcat_example)
#' associations(siamcat_example) <- list(
#'     'assoc.results'=associations(siamcat_example),
#'     'assoc.param'=assoc_param(siamcat_example))
setGeneric("associations<-", function(x, value)
    standardGeneric("associations<-"))
#' @rdname assign-associations
#' @aliases associations<-
setReplaceMethod("associations", c("siamcat", "list"), function(x, value) {
    x@associations <- value
    validObject(x)
    return(x)
})


################################################################################
#' Assign a new list containing normalziation parameters and normalized
#' features to a SIAMCAT object
#'
#' @usage norm_feat(x) <- value
#'
#' @param x an object of class \link{siamcat-class}
#' @param value a list containing normaliation parameters and features
#' @export
#' @docType methods
#' @rdname assign-norm_feat
#' @aliases assign-norm_feat
#' @return none
#'
#' @examples
#' data(siamcat_example)
#' norm_feat(siamcat_example) <- norm_feat(siamcat_example)
setGeneric("norm_feat<-", function(x, value)
    standardGeneric("norm_feat<-"))
#' @rdname assign-norm_feat
#' @aliases norm_feat<-
setReplaceMethod("norm_feat", c("siamcat", "list"), function(x, value) {
    x@norm_feat <- value
    validObject(x)
    return(x)
})

################################################################################
#' @title Assign a new list containing a cross-validation split to a
#' SIAMCAT object
#'
#' @usage data_split(x) <- value
#'
#' @param x an object of class \link{siamcat-class}
#' @param value list containing a cross-validation split
#' @export
#' @docType methods
#' @rdname assign-data_split
#' @aliases assign-data_split
#' @return none
#'
#' @examples
#' data(siamcat_example)
#' data_split(siamcat_example) <- data_split(siamcat_example)
setGeneric("data_split<-", function(x, value)
    standardGeneric("data_split<-"))
#' @rdname assign-data_split
#' @aliases data_split<-
setReplaceMethod("data_split", c("siamcat", "list"), function(x, value) {
    x@data_split <- value
    validObject(x)
    return(x)
})

################################################################################
#' Assign a new list containing trained models to a SIAMCAT object
#'
#' @usage model_list(x) <- value
#'
#' @param x an object of class \link{siamcat-class}
#' @param value list containing trained models, type of models and of features
#' @export
#' @docType methods
#' @rdname assign-model_list
#' @aliases assign-model_list
#' @return none
#'
#' @examples
#' data(siamcat_example)
#' siamcat_example <- train.model(siamcat_example, method='lasso')
#' model_list(siamcat_example) <- model_list(siamcat_example)
setGeneric("model_list<-", function(x, value)
    standardGeneric("model_list<-"))
#' @rdname assign-model_list
#' @aliases model_list<-
setReplaceMethod("model_list", c("siamcat", "list"), function(x, value) {
    x@model_list <- value
    validObject(x)
    return(x)
})

################################################################################
#' Assign a new matrix with predictions to a SIAMCAT object
#'
#' @usage pred_matrix(x) <- value
#'
#' @param x an object of class \link{siamcat-class}
#' @param value a matrix containing predictions
#' @export
#' @docType methods
#' @rdname assign-pred_matrix
#' @aliases assign-pred_matrix
#' @return none
#'
#' @examples
#' data(siamcat_example)
#' pred_matrix(siamcat_example) <- pred_matrix(siamcat_example)
setGeneric("pred_matrix<-", function(x, value)
    standardGeneric("pred_matrix<-"))
#' @rdname assign-pred_matrix
#' @aliases pred_matrix<-
setReplaceMethod("pred_matrix", c("siamcat", "matrix"), function(x, value) {
    x@pred_matrix <- value
    validObject(x)
    return(x)
})


################################################################################
#' Assign a new list with evaluation data to a SIAMCAT object
#'
#' @usage eval_data(x) <- value
#'
#' @param x an object of class \link{siamcat-class}
#' @param value a list of evaluation data
#' @export
#' @docType methods
#' @rdname assign-eval_data
#' @aliases assign-eval_data
#' @return none
#'
#' @examples
#' data(siamcat_example)
#' eval_data(siamcat_example) <- eval_data(siamcat_example)
setGeneric("eval_data<-", function(x, value)
    standardGeneric("eval_data<-"))
#' @rdname assign-eval_data
#' @aliases eval_data<-
setReplaceMethod("eval_data", c("siamcat", "list"), function(x, value) {
    x@eval_data <- value
    validObject(x)
    return(x)
})
