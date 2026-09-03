#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL

#' @title Reset the label of a siamcat object to a new metadata column, to a vector, or to a test label
#'
#' @description This function reset the label of a siamcat object to a new metadata column, to a vector, or to a test label
#'
#' @usage reset.label(siamcat, "disease_state")
#'
#' @param siamcat an object of class \link{siamcat-class}
#'
#' @param label named vector to create the label or the name of the metadata
#' column that will be used to create the label. If \code{NULL}, a placeholder TEST label is used.
#' 
#' @param ... all additional named parameters are passed to the create.label function
#'

reset.label <- function(siamcat, label=NULL, ...) {
    meta <- meta(siamcat)
    otu_table <- orig_feat(siamcat)
    taxa <- tax_table(siamcat@phyloseq)

    if (is.character(label) && length(label) == 1) {
        label <- create.label(label, meta=meta, ...)
    } else if (is.vector(label)) {
        label <- create.label(label, ...)
    } else {
        stop("The parameter label is not a charachter vector of length 1 or a vector.")
    }
    
    sc <- siamcat(feat = otu_table, label = label, meta = meta, taxonomy = taxa)
    return(sc)
}