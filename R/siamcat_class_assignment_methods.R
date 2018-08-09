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
#' @return none
#'
#' @examples
#' data(siamcat_example)
#' physeq(siamcat_example) <- physeq(siamcat_example)
setGeneric("physeq<-", function(x, value)
    standardGeneric("physeq<-"))
#' @rdname assign-physeq
#' @aliases physeq<-
setMethod("physeq<-", c("siamcat", "phyloseq"), function(x, value) {
    siamcat(
        'phyloseq'=value,
        'label'=x@label,
        x@filt_feat,
        x@associations,
        x@norm_feat,
        x@data_split,
        x@model_list,
        x@pred_matrix,
        x@eval_data,
        validate=FALSE,
        verbose=0
    )
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
#' @return none
#'
#' @examples
#' data(siamcat_example)
#' features(siamcat_example) <- features(siamcat_example)
setGeneric("features<-", function(x, value)
    standardGeneric("features<-"))
#' @rdname assign-features
#' @aliases features<-
setMethod("features<-", c("siamcat", "otu_table"), function(x, value) {

    if (!is.null(norm_feat(x, verbose=0))){
        norm_feat(x) <- new('norm_feat', norm.feat=value,
            norm.param=norm_params(x))
    } else if (!is.null(filt_feat(x, verbose=0))){
        filt_feat(x) <- new('filt_feat', filt.feat=value,
            filt.param=filt_params(x))
    } else {
        orig_feat(x) <- value
    }
    return(x)
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
#' @return none
#'
#' @examples
#' data(siamcat_example)
#' label(siamcat_example) <- label(siamcat_example)
setGeneric("label<-", function(x, value)
    standardGeneric("label<-"))
#' @rdname assign-label
#' @aliases label<-
setMethod("label<-", c("siamcat", "label"), function(x, value) {
    siamcat(
        'phyloseq'=x@phyloseq,
        'label'=value,
        x@filt_feat,
        x@associations,
        x@norm_feat,
        x@data_split,
        x@model_list,
        x@pred_matrix,
        x@eval_data,
        validate=FALSE,
        verbose=0
    )
})

###############################################################################
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
#' filt_feat(siamcat_example) <- filt_feat(siamcat_example)
setGeneric("filt_feat<-", function(x, value)
    standardGeneric("filt_feat<-"))
#' @rdname assign-filt_feat
#' @aliases filt_feat<-
setMethod("filt_feat<-", c("siamcat", "filt_feat"), function(x, value) {
    siamcat(
        'label'=x@label,
        'phyloseq'=x@phyloseq,
        value,
        x@associations,
        x@norm_feat,
        x@data_split,
        x@model_list,
        x@pred_matrix,
        x@eval_data,
        validate=FALSE,
        verbose=0
    )
})


###############################################################################
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
setGeneric("associations<-", function(x, value)
    standardGeneric("associations<-"))
#' @rdname assign-associations
#' @aliases associations<-
setMethod("associations<-", c("siamcat", "associations"), function(x, value) {
    siamcat(
        'phyloseq'=x@phyloseq,
        'label'=x@label,
        x@filt_feat,
        'assoc'=associations(value),
        x@norm_feat,
        x@data_split,
        x@model_list,
        x@pred_matrix,
        x@eval_data,
        validate=FALSE,
        verbose=0
    )
})


###############################################################################
#' Assign a new norm_feat object to \code{x}
#'
#' @usage norm_feat(x) <- value
#'
#' @param x an object of class \link{siamcat-class}
#' @param value an norm_feat object
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
setMethod("norm_feat<-", c("siamcat", "norm_feat"), function(x, value) {
    siamcat(
        'label'=x@label,
        'phyloseq'=x@phyloseq,
        x@filt_feat,
        x@associations,
        value,
        x@data_split,
        x@model_list,
        x@pred_matrix,
        x@eval_data,
        validate=FALSE,
        verbose=0
    )
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
#' @return none
#'
#' @examples
#' data(siamcat_example)
#' data_split(siamcat_example) <- data_split(siamcat_example)
setGeneric("data_split<-", function(x, value)
    standardGeneric("data_split<-"))
#' @rdname assign-data_split
#' @aliases data_split<-
setMethod("data_split<-", c("siamcat", "data_split"), function(x, value) {
    siamcat(
        'label'=x@label,
        'phyloseq'=x@phyloseq,
        x@filt_feat,
        x@associations,
        x@norm_feat,
        value,
        x@model_list,
        x@pred_matrix,
        x@eval_data,
        validate=FALSE,
        verbose=0
    )
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
#' @return none
#'
#' @examples
#' data(siamcat_example)
#' model_list(siamcat_example) <- model_list(siamcat_example)
setGeneric("model_list<-", function(x, value)
    standardGeneric("model_list<-"))
#' @rdname assign-model_list
#' @aliases model_list<-
setMethod("model_list<-", c("siamcat", "model_list"), function(x, value) {
    siamcat(
        'label'=x@label,
        'phyloseq'=x@phyloseq,
        x@filt_feat,
        x@associations,
        x@norm_feat,
        x@data_split,
        value,
        x@pred_matrix,
        x@eval_data,
        validate=FALSE,
        verbose=0
    )
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
#' @return none
#'
#' @examples
#' data(siamcat_example)
#' pred_matrix(siamcat_example) <- pred_matrix(siamcat_example)
setGeneric("pred_matrix<-", function(x, value)
    standardGeneric("pred_matrix<-"))
#' @rdname assign-pred_matrix
#' @aliases pred_matrix<-
setMethod("pred_matrix<-", c("siamcat", "matrix"), function(x, value) {
    siamcat(
        'label'=x@label,
        'phyloseq'=x@phyloseq,
        x@filt_feat,
        x@associations,
        x@norm_feat,
        x@data_split,
        x@model_list,
        value,
        x@eval_data,
        validate=FALSE,
        verbose=0
    )
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
#' @return none
#'
#' @examples
#' data(siamcat_example)
#' eval_data(siamcat_example) <- eval_data(siamcat_example)
setGeneric("eval_data<-", function(x, value)
    standardGeneric("eval_data<-"))
#' @rdname assign-eval_data
#' @aliases eval_data<-
setMethod("eval_data<-", c("siamcat", "list"), function(x, value) {
    siamcat(
        'label'=x@label,
        'phyloseq'=x@phyloseq,
        x@filt_feat,
        x@associations,
        x@norm_feat,
        x@data_split,
        x@model_list,
        x@pred_matrix,
        value,
        validate=FALSE,
        verbose=0
    )
})
