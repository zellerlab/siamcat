#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' Build siamcat-class objects from their components.
#' @title siamcat
#' @name siamcat
#' @description Function to construct an object of class \link{siamcat-class}
#' @param ... list of arguments needed in order to construct a SIAMCAT object
#' @return A new \link{siamcat-class} object
#' @export
#' @examples
#' # example with package data
#' fn.in.feat    <- system.file('extdata',
#'     'feat_crc_zeller_msb_mocat_specI.tsv',
#'     package = 'SIAMCAT')
#' fn.in.label <- system.file('extdata',
#'     'label_crc_zeller_msb_mocat_specI.tsv',
#'     package = 'SIAMCAT')
#' fn.in.meta    <- system.file('extdata',
#'     'num_metadata_crc_zeller_msb_mocat_specI.tsv',
#'     package = 'SIAMCAT')
#'
#' feat    <- read.features(fn.in.feat)
#' label <- read.labels(fn.in.label)
#' meta    <- read.meta(fn.in.meta)
#' siamcat <- siamcat(feat, label, meta)
siamcat <- function(..., feat=NULL, label=NULL, meta=NULL, phyloseq=NULL,
        validate=TRUE, verbose=3) {

    if (is.null(phyloseq) && is.null(feat)){
        stop(paste0('SIAMCAT needs either a feature matrix or a phyloseq',
            ' object!!! Exiting...'))
    }

    other.args <- list(...)

    # keep case info if the user wants to create the label from metadata
    if ('case' %in% names(other.args)){
        case <- other.args$case
    } else {
        case <- NULL
    }

    # Remove names from arglist. Will replace them based on their class
    names(other.args) <- NULL

    if (!is.null(other.args)){
        # ignore all but component data classes.
        component_classes <- get.component.classes("both")

        for (argNr in seq_along(other.args)) {
            classOfArg <- class(other.args[[argNr]])[1]
            if (classOfArg %in% names(component_classes)) {
                names(other.args)[argNr] <- component_classes[classOfArg]
            }
        }
    }

    if (!is.null(phyloseq)){
        if (!is.null(feat)){
            stop(paste0('Both features matrix and phyloseq object provided. ',
                'Please provide only one of them!'))
        }
        if (class(phyloseq) != 'phyloseq'){
            stop('Please provide an object of class phyloseq for SIAMCAT!')
        }
        feat <- otu_table(phyloseq)
        if (!is.null(sample_data(phyloseq, errorIfNULL=FALSE))) {
            meta <- sample_data(phyloseq)
        } else {
            meta <- NULL
        }
    } else {
        feat <- validate.features(feat)
        meta <- validate.metadata(meta)
    }

    if (any(vapply(names(other.args), is.component.class, "phyloseq",
        FUN.VALUE = logical(1)))){
            arglistphyloseq <- other.args[vapply(names(other.args),
                is.component.class,
                "phyloseq", FUN.VALUE = logical(1))]
            arglistphyloseq <- list(arglistphyloseq, 'otu_table'=feat,
                'sam_data'=meta)
    } else {
        arglistphyloseq <- list('otu_table'=feat, 'sam_data'=meta)
    }
    other.args$phyloseq <- do.call("new", c(list(Class = "phyloseq"),
        arglistphyloseq))

    label <- validate.label(label, feat, meta, case, verbose)

    if (is.null(other.args$orig_feat)) {
        other.args$orig_feat <- orig_feat(otu_table(other.args$phyloseq))
    }

    other.args$label <- label

    other.args <-
        other.args[vapply(names(other.args),
            is.component.class,
            "siamcat",
            FUN.VALUE = logical(1))]
    sc <- do.call("new", c(list(Class = "siamcat"), other.args))

    if (validate){
        sc <- validate.data(sc, verbose=verbose)
    } else {
        if (verbose > 0){
            message(paste0('### Not validating the SIAMCAT object!!!\n',
            '\tPlease be advised that some functions may not work correctly!'))
        }
    }
    return(sc)
}

# source: https://github.com/joey711/phyloseq/blob/master/R/phyloseq-class.R
#' Show the component objects classes and slot names.
#' @keywords internal
#' @return list of component classes
get.component.classes <- function(class) {
    # define classes vector the names of component.classes needs to be the slot
    # names to match getSlots / splat

    #slot names
    component.classes.siamcat <-
        c(
            "model_list",
            "orig_feat",
            "label",
            "norm_param",
            "data_split",
            "phyloseq",
            "eval_data",
            "pred_matrix"
        )

    #class names
    names(component.classes.siamcat) <-
        c(
            "model_list",
            "orig_feat",
            "label",
            "norm_param",
            "data_split",
            "phyloseq",
            "eval_data",
            "pred_matrix"
        )

    #slot names
    component.classes.phyloseq <-
        c("otu_table", "sam_data", "phy_tree",
            "tax_table", "refseq")

    #class names
    names(component.classes.phyloseq) <-
        c("otu_table",
            "sample_data",
            "phylo",
            "taxonomyTable",
            "XStringSet")

    if (class == "siamcat") {
        return(component.classes.siamcat)
    } else if (class == "phyloseq") {
        return(component.classes.phyloseq)
    } else if (class == "both") {
        return(c(
            component.classes.siamcat,
            component.classes.phyloseq
        ))
    }
}

# Returns TRUE if x is a component class, FALSE otherwise.
#' @keywords internal
is.component.class = function(x, class) {
    x %in% get.component.classes(class)
}

# check and convert feature object to otu_table
#' @keywords internal
validate.features <- function(feat){
    # check if NA and stop
    if (is.null(feat)){
        stop('SIAMCAT needs features!!! Exiting...')
    }
    # check class of feature input
    if (class(feat) == 'otu_table'){
        # can either be an otu_table (then to nothing)
        return(feat)
    } else if (class(feat) == 'matrix'){
        # or a matrix (then check if it is numeric or not)
        # and convert to otu_table
        if (any(!is.numeric(feat))){
            stop(paste0('SIAMCAT expects numerical features!.\n',
            'Please check your feature matrix! Exiting...'))
        }
        feat <- otu_table(feat, taxa_are_rows=TRUE)
        return(feat)
    } else if (class(feat) == 'data.frame'){
        # or a dataframe (then do the same as above)
        if (any(!is.numeric(unlist(feat)))){
            stop(paste0('SIAMCAT expects numerical features!.\n',
            'Please check your feature data.frame! Exiting...'))
        }
        feat <- otu_table(feat, taxa_are_rows=TRUE)
        return(feat)
    }
}

# check label object
#' @keywords internal
validate.label <- function(label, feat, meta, case, verbose){
    # if NA, return simple label object which contains only one class
    if (is.null(label)){
        message(paste0('No label information given! Generating SIAMCAT object ',
        'with placeholder label!\n\tThis SIAMCAT object is not suitable for ',
        'the complete workflow...'))
        label <- label(list(label = sample(c(1, -1),
            size=ncol(feat), replace=TRUE),
            type="BINARY", info=c('TEST-A'=1, 'TEST-B'=-1)))
        names(label$label) <- colnames(feat)
        return(label)
    }
    if (class(label) == "label"){
        return(label)
    } else if (is.character(label)){
        if(is.null(meta)) stop('Metadata needed to generate label! Exiting...')
        if(is.null(case)) stop('Case information needed! Exiting...')
        label <- create.label.from.metadata(meta, label, case=case,
            verbose=verbose)
        label <- label(label)
        return(label)
    }
}

# check meta-data object
#' @keywords internal
validate.metadata <- function(meta){
    if (is.null(meta)){
        return(NULL)
    }
    if (class(meta) == 'sample_data'){
        return(meta)
    }
    if (class(meta) == 'data.frame'){
        meta <- sample_data(meta)
        return(meta)
    }
}
