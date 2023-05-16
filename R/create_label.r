#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Create a label list
#'
#' @description This function creates a label object from metadata
#' or an atomic vector
#'
#' @usage create.label(label, case, meta=NULL, control=NULL, 
#' p.lab = NULL, n.lab = NULL, remove.meta.column=FALSE, verbose=1)
#'
#' @param label named vector to create the label or the name of the metadata
#' column that will be used to create the label
#'
#' @param case name of the group that will be used as a positive label. If the
#' variable is binary, the other label will be used as a negative one. If the
#' variable has multiple values, all the other values will be used a negative
#' label (testing one vs rest).
#'
#' @param meta metadata dataframe object or an object of class 
#' \link[phyloseq]{sample_data-class}
#'
#' @param control name of a label or vector with names that will be used as a
#' negative label. All values that are nor equal to case and control will be
#' dropped. Default to NULL in which case: If the variable is binary, the value
#' not equal to case will be used as negative. If the variable has multiple
#' values, all the values not equal to cases will be used a negative label
#' (testing one vs rest).
#'
#' @param p.lab name of the positive group (useful mostly for visualizations).
#' Default to NULL in which case the value of the positive group will be used.
#'
#' @param n.lab name of the negative group (useful mostly for visualizations).
#' Default to NULL in which case the value of the negative group will be used
#' for binary variables and "rest" will be used for variables with multiple
#' values.
#'
#' @param remove.meta.column boolean indicating if the label column in the
#' metadata should be retained. Please note that if this is set to
#' \code{TRUE}, the function will return a list as result. Defaults to
#' \code{FALSE}
#'
#' @param verbose integer, control output: \code{0} for no output at all, 
#' \code{1} for only information about progress and success, \code{2} for 
#' normal level of information and \code{3} for full debug information, 
#' defaults to \code{1}
#'
#' @keywords create.label
#'
#' @details The function creates a list to be used as label in a SIAMCAT
#' object. Mainly for interal use, but it can be used to customize your label
#' (\code{p.lab} and \code{n.lab} will be used as labels during plotting,
#' for example).
#'
#' The input for the function can be either a named vector encoding the label
#' or the name of a column in the metadata (needs to be provided as well) which
#' contains the label information.
#'
#' @export
#' 
#' @encoding UTF-8
#'
#' @return return either \itemize{
#' \item a list to be used in a SIMCAT object \strong{OR}
#' \item a list with entries \code{meta} and \code{label},
#' if \code{remove.meta.column} is set to \code{TRUE}
#' }
#'
#' @examples
#' data('meta_crc_zeller')
#'
#' label <- create.label(label='Group', case='CRC', meta=meta.crc.zeller)
create.label <- function(label, case=NULL, meta=NULL, control=NULL,
    p.lab=NULL, n.lab=NULL, remove.meta.column=FALSE, verbose=1) {
    if (verbose > 1)
        message("+ starting create.label")
    s.time <- proc.time()[3]

    #if metadata has been supplied and the label is of length 1
    if (!is.null(meta) & length(label) == 1){
        if (!label %in% colnames(meta))
            stop("Column ", label, " not found in the metadata\n")
        if (is(meta,'sample_data')){
            label.vec <- as.matrix(as.data.frame(meta[,label]))[,1]
        } else if (is.data.frame(meta)){
            label.vec <- meta[, label]
        } else {
            stop(paste0('Please provide the metadata either as a data.frame',
                ' or a sample_data object!'))
        }
        names(label.vec) <- rownames(meta)
    #if the label is a vector
    } else if (is.atomic(label) & length(label) > 1){
        label.vec <- label
        if (is.null(names(label.vec))){
            stop('The label vectors needs to be named!')
        }
        if (is.factor(label.vec)){
            names.old <- names(label.vec)
            label.vec <- as.character(label.vec)
            names(label.vec) <- names.old
        }
    } else {
        stop(paste0('Could not interpret your input.\n Please provide ',
        'either a column name and metadata or a label vector.\nExiting...!'))
    }

    # remove NAs in the label
    if (any(is.na(label.vec))){
        if (verbose > 1) message(paste0('+ removing ', sum(is.na(label.vec)),
            ' instances of NA in the label'))
        label.vec <- label.vec[!is.na(label.vec)]
    }

    # find out if binary or regression!
    n.groups <- length(unique(label.vec))

    if (is.character(label.vec) |
        is.factor(label.vec) |
        length(n.groups) == 2 | !is.null(case)){

        if (!is.character(label.vec)){
            x <- names(label.vec)
            label.vec <- as.character(label.vec)
            names(label.vec) <- x
        }

        # get different groups
        groups <- unique(label.vec)
        if (is.null(case)) stop("Case information is needed for binary label")
        ### checking case
        if(!all(case %in% groups)){
            stop("The chosen label does not contain values: ",
                paste(case,collapse=","),"\nInstead, contains: ",
                paste(groups, collapse=','))
        }

        ### checking control
        if (is.null(control)) {
            if((length(groups)-length(case))>1){
                control <- "rest"
            }else{
                control <- setdiff(groups, case)
            }
        } else {
            if(!control%in%groups){
                stop("The chose label does not contain value:",control,
                    "\nInstead, contains: ", paste(groups, collapse=','))
            }
            ### dropping unused values
            if(any(!groups%in%c(case, control))){
                label.vec <- label.vec[which(label.vec%in%c(case, control))]
                warning("Dropping values: ",
                        paste(groups[which(!groups%in%c(case, control))],
                            collapse=', '), '\n')
            }
        }

        # message status
        if (verbose > 0)
            message("Label used as case:\n   ",paste(case,collapse=","),
                    "\nLabel used as control:\n   ",paste(control,collapse=","))

        # create new label.object
        label.new <- list(label = rep(-1, length(label.vec)))

        n.lab <- ifelse(is.null(n.lab), gsub("[_.-]", ".", control), n.lab)
        p.lab <- ifelse(is.null(p.lab),
                        ifelse(length(case) > 1, 'Case',
                            gsub("[_.-]", ".", case)), p.lab)

        info <- c(-1, 1)
        names(info) <- c(n.lab, p.lab)

        names(label.new$label) <- names(label.vec)

        if (length(case) > 1){
            label.new$label[which(label.vec %in% case)] <- 1
        } else {
            label.new$label[which(label.vec == case)] <- 1
        }

        label.new$info <- info
        label.new$type <- "BINARY"

        label.new <- label.new
    } else if (is.double(label.vec)){
        if (!is.null(case) | !is.null(control)){
            warning(paste0("Case and control parameters will be ignored for",
                        " continuous labels!"))
        }

        label.new <- list(label=label.vec, info=range(label.vec), 
                            type='CONTINUOUS')

    }

    e.time <- proc.time()[3]
    if (verbose > 0)
        message(paste(
            "+ finished create.label.from.metadata in",
            formatC(e.time - s.time, digits = 3),
            "s"
        ))
    if (!remove.meta.column){
        return(label.new)
    } else {
        meta.new <- meta[,-which(colnames(meta) == label)]
        return(list(label=label.new,
                    meta=meta.new))
    }
}
