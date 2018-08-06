#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title create a label object from metadata
#'
#' @description This function creates a label object from metadata
#'
#' @usage create.label.from.metadata(meta, column, case,
#' control=NULL, p.lab = NULL, n.lab = NULL, verbose=1)
#'
#' @param meta metadata as read by \link{read.meta}
#' of \link[phyloseq]{sample_data-class}
#'
#' @param column name of column that will be used
#' to create the label
#'
#' @param case name of a label that will be used as a positive label. If the
#' variable is binary, the other label will be used as a negative one. If the
#' variable has multiple values, all the other values will be used a negative
#' label (testing one vs rest).
#'
#' @param control name of a label or vector with names that will be used as a
#' negative label. All values that are nor equal to case and control will be
#' dropped. Default to NULL in which case: If the variable is binary, the value
#' not equal to case will be used as negative. If the variable has multiple
#' values, all the values not equal to cases will be used a negative label
#' (testing one vs rest).
#'
#' @param p.lab name of the positive label (useful mostly for visualizations).
#' Default to NULL in which case the value of the positive label will be used.
#'
#' @param n.lab name of the negative label (useful mostly for visualizations).
#' Default to NULL in which case the value of the negative label will be used
#' for binary variables and "rest" will be used for variables with multiple
#' values.
#'
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'     for only information about progress and success, \code{2} for normal
#'     level of information and \code{3} for full debug information,
#'     defaults to \code{1}
#'
#' @keywords create.label.from.metadata
#'
#' @return an object of class \link{label-class}
#'
#' @examples
#'     data(siamcat_example)
#'     label <- create.label.from.metadata(meta(siamcat_example),"fobt",
#'     case = 1, control = 0)
#'
#' @export
create.label.from.metadata <- function(meta, column, case, control = NULL,
                                       verbose=1) {
    if (verbose > 1)
        message("+ starting create.label.from.metadata")

    s.time <- proc.time()[3]

    if (!column %in% colnames(meta))
        stop("ERROR: Column ", column, " not found in the metadata\n")

    if (class(meta) == 'sample_data'){
        metaColumn <- vapply(meta[, column], as.character,
            FUN.VALUE = character(nrow(meta)))
    } else if (class(meta) == 'data.frame'){
        metaColumn <- vapply(meta[, column], as.character,
            FUN.VALUE = character(1))
    } else {
        stop('Please provide either a data.frame or a sample_data object!')
    }

    names(metaColumn) <- rownames(meta)

    # remove NAs in label column
    if (any(is.na(metaColumn))){
        if (verbose > 1) message(paste0('+ removing ', sum(is.na(metaColumn)),
            ' instances of NA in chosen column'))
        metaColumn <- metaColumn[!is.na(metaColumn)]
    }

    labels <- unique(metaColumn)

    ### checking case
    if(!all(case %in% labels)){
        stop("Column ", column, " does not contain values:",
                paste(case,collapse=","),"\n")
    }


    ### checking control
    if (is.null(control)) {
        if((length(labels)-length(case))>1){
            control <- "rest"
        }else{
            control <- setdiff(labels, case)
        }
    }else{
        if(!control%in%labels){
            stop("Column ", column, " does not contain value:",control,"\n")
        }
        ### dropping unused values
        if(any(!labels%in%c(case, control))){
            metaColumn <- metaColumn[which(metaColumn%in%c(case, control))]
            warning("Dropping values: ",
                labels[which(!labels%in%c(case, control))])
        }
    }


    if (verbose > 0)
            message("Label used as case:\n   ",paste(case,collapse=","),
                "\nLabel used as control:\n   ",paste(control,collapse=","))
    label <- list(label = rep(-1, length(metaColumn)))

    n.lab <- gsub("[_.-]", " ", control)
    if (length(case) > 1) {
        p.lab <- "Case"
    } else {
        p.lab <- gsub("[_.-]", " ", case)
    }

    info <- c(-1, 1)
    names(info) <- c(n.lab, p.lab)

    names(label$label) <- names(metaColumn)

    if (length(case) > 1){
        label$label[which(metaColumn %in% case)] <- 1
    } else {
        label$label[which(metaColumn == case)] <- 1
    }

    label$type <- "BINARY"
    label$info <- info

    labelRes <- label(label)

        e.time <- proc.time()[3]
    if (verbose > 0)
        message(paste(
            "+ finished create.label.from.metadata in",
            formatC(e.time - s.time, digits = 3),
            "s"
        ))
    return(labelRes)
}
