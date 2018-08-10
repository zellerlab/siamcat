#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' Build siamcat-class objects from their components.
#' @title siamcat
#' @name siamcat
#' @description Function to construct an object of class \link{siamcat-class}
#' @param ... additional arguments
#' @param feat feature information for SIAMCAT (see details)
#' @param label label information for SIAMCAT (see details)
#' @param meta (optional) metadata information for SIAMCAT (see details)
#' @param phyloseq (optional) a phyloseq object for the creation of an SIAMCAT
#'  object (see details)
#' @param validate boolean, should the newly constructed SIAMCAT object be
#'  validated? defaults to TRUE (\strong{we strongly recommend against
#'  setting this parameter to FALSE})
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'     for only information about progress and success, \code{2} for normal
#'     level of information and \code{3} for full debug information,
#'     defaults to \code{1}
#' @return A new \link{siamcat-class} object
#' @export
#' @details This functions creates a SIAMCAT object (see \link{siamcat-class}).
#'  In order to do so, the function needs \itemize{
#'  \item feat the feature information for SIAMCAT, should be either a matrix,
#'      a data.frame, or a \link[phyloseq]{otu_table}. The columns should
#'      correspond to the different samples (e.g. patients) and the rows the
#'      different features (e.g. taxa). Columns and rows should be named.
#'  \item meta metadata information for the different samples in the feature
#'      matrix. Metadata is optional for the SIAMCAT workflow. Should be
#'      either a data.frame (with the rownames corresponding to the sample
#'      names of the feature matrix) or an object of class
#'      \link[phyloseq]{sample_data}
#'  \item phyloseq Alternatively to supplying both feat and meta, SIAMCAT can
#'      also work with a phyloseq object containing an otu_table and other
#'      optional slots (like sample_data for meta-variables).}
#'
#' Notice: do supply \strong{either} the feature information as
#' matrix/data.frame/otu_table (and optionally metadata) \strong{or} a
#' phyloseq object, but not both.
#'
#' The label information for SIAMCAT can take several forms:\itemize{
#'  \item metadata column: if there is metadata (either via meta or as
#'      sample_data in the phyloseq object), the label object can be created
#'      by taking the information in a specific metadata column. In order to
#'      do so, \code{label} should be the name of the column, and \code{case}
#'      should indicate which group(s) should be the positive group(s). A
#'      typical example could look like that:
#'
#'\code{siamcat <- siamcat(feat=feat.matrix, meta=metadata,
#'                         label='DiseaseState', case='CRC')}
#'
#'      for the construction of a label to predict CRC status (which is encoded
#'      in the column \code{"DiseaseState"} of the metadata). For more control
#'      (e.g. specific labels for plotting or specific control state), the
#'      label can also be created outside of the \code{siamcat} function using
#'      the \link{create.label} function (see below).
#'  \item named vector: the label can also be supplied as named vector which
#'      encodes the label either as characters (e.g. "Healthy" and "Diseased"),
#'      as factor, or numerically (e.g. -1 and 1). The vector must be named
#'      with the names of samples (corresponding to the samples in features).
#'      Also here, the information about the positive group(s) is needed via
#'      the \code{case} parameter. Internally, the vector is given to the
#'      \link{create.label} function (see for more details).
#'  \item label object: A label object can be created with the
#'      \link{create.label} function or by reading a dedicated label file
#'      with \link{read.labels}.
#'  }
#' @examples
#' # example with package data
#' data("feat_crc_zeller", package="SIAMCAT")
#' data("meta_crc_zeller", package="SIAMCAT")
#'
#' siamcat <- siamcat(feat=feat.crc.zeller,
#'                    meta=meta.crc.zeller,
#'                    label='Group',
#'                    case='CRC')
siamcat <- function(..., feat=NULL, label=NULL, meta=NULL, phyloseq=NULL,
        validate=TRUE, verbose=3) {

    if (is.null(phyloseq) && is.null(feat)){
        stop(paste0('SIAMCAT needs either a feature matrix or a phyloseq',
            ' object!!! Exiting...'))
    }

    # optional arguments
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

    # if phyloseq object has been given
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
        # get all other slots from the phyloseq object
        for (x in setdiff(get.component.classes('phyloseq'),
            c('otu_table', 'sam_data'))){
                if (!is.null(access(phyloseq, x))){
                    other.args[[x]] <- access(phyloseq, x)
                }
            }
    }

    # validate features and metadata
    feat <- validate.features(feat)
    meta <- validate.metadata(meta)

    # make Phyloseq object properly
    if (any(vapply(names(other.args), is.component.class, "phyloseq",
        FUN.VALUE = logical(1)))){
            arglistphyloseq <- other.args[vapply(names(other.args),
                is.component.class,
                "phyloseq", FUN.VALUE = logical(1))]
            arglistphyloseq$otu_table = feat
            arglistphyloseq$sam_data = meta
    } else {
        arglistphyloseq <- list('otu_table'=feat, 'sam_data'=meta)
    }
    other.args$phyloseq <- do.call("new", c(list(Class = "phyloseq"),
        arglistphyloseq))

    # label object
    label <- validate.label(label, feat, meta, case, verbose)
    other.args$label <- label

    # any other slots
    other.args <-
        other.args[vapply(names(other.args),
            is.component.class,
            "siamcat",
            FUN.VALUE = logical(1))]
    sc <- do.call("new", c(list(Class = "siamcat"), other.args))

    # validate
    if (validate){
        sc <- validate.data(sc, verbose=verbose)
    } else {
        if (verbose > 0){
            warning(paste0('### Not validating the SIAMCAT object!!!\n',
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
            "phyloseq",
            "label",
            "filt_feat",
            "associations",
            "norm_feat",
            "data_split",
            "model_list",
            "pred_matrix",
            "eval_data"
        )
    #class names
    names(component.classes.siamcat) <-
        c(
            "phyloseq",
            "label",
            "filt_feat",
            "associations",
            "norm_feat",
            "data_split",
            "model_list",
            "pred_matrix",
            "eval_data"
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
        # can either be an otu_table (then do nothing)
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
        warning(paste0('No label information given! Generating SIAMCAT object ',
        'with placeholder label!\n\tThis SIAMCAT object is not suitable for ',
        'the complete workflow...'))
        label <- label(list(label = rep(-1, ncol(feat)),
             info=c('TEST'=-1), type="TEST"))
        names(label$label) <- colnames(feat)
        return(label)
    } else if (class(label) == "label"){
        return(label)
    } else if (is.character(label) & length(label) == 1){
        if(is.null(meta)) stop('Metadata needed to generate label! Exiting...')
        if(is.null(case)) stop('Case information needed! Exiting...')
        label <- create.label(meta=meta, label=label, case=case,
            verbose=verbose)
        label <- label(label)
        return(label)
    } else if (is.atomic(label)) {
        if(is.null(case)) stop('Case information needed! Exiting...')
        label <- create.label(label=label, case=case)
    } else {
        stop(paste0('Cannot interpret the label object!\nPlease ',
            'provide either a label object, a column in your metadata, or a
            vector to distinguish cases and controls!.'))
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
