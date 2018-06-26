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
#' 'num_metadata_crc_zeller_msb_mocat_specI.tsv',
#'     package = 'SIAMCAT')
#'
#' feat    <- read.features(fn.in.feat)
#' label <- read.labels(fn.in.label)
#' meta    <- read.meta(fn.in.meta)
#' siamcat <- siamcat(feat, label, meta)
siamcat <- function(...) {
    arglist <- list(...)
    
    # Remove names from arglist. Will replace them based on their class
    names(arglist) <- NULL
    
    # ignore all but component data classes.
    component_classes <- get.component.classes("both")
    
    for (argNr in seq_along(arglist)) {
        classOfArg <- class(arglist[[argNr]])[1]
        if (classOfArg %in% names(component_classes)) {
            names(arglist)[argNr] <- component_classes[classOfArg]
        }
    }
    
    if (is.null(arglist$phyloseq)) {
        arglistphyloseq <-
            arglist[vapply(names(arglist),
                is.component.class,
                "phyloseq",
                FUN.VALUE = logical(1))]
        arglist$phyloseq <-
            do.call("new", c(list(Class = "phyloseq"),
                arglistphyloseq))
    }
    if (is.null(arglist$orig_feat)) {
        arglist$orig_feat <- orig_feat(otu_table(arglist$phyloseq))
    }
    arglist <-
        arglist[vapply(names(arglist),
            is.component.class,
            "siamcat",
            FUN.VALUE = logical(1))]
    sc <- do.call("new", c(list(Class = "siamcat"), arglist))
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
