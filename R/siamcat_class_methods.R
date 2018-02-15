#!/usr/bin/Rscript
###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# R flavor
# EMBL Heidelberg 2012-2018
# GNU GPL 3.0
###

#' @title Filter samples from \code{siamcat@label}
#' @description This functions filters \code{siamcat@label}.
#' @param siamcat an object of class \link{siamcat}
#' @param ids names of samples to be left in the \code{siamcat@label}
#' @keywords filter.label
#' @export
#' @return siamcat an object of class \link{siamcat}
#'
filter.label <- function(siamcat,ids, verbose=1){
  labels_new      <- new("label", label = siamcat@label@label[ids],
                         header = siamcat@label@header,
                         info = siamcat@label@info,
                         positive.lab = siamcat@label@positive.lab,
                         negative.lab = siamcat@label@negative.lab,
                         n.lab = siamcat@label@n.lab,
                         p.lab = siamcat@label@p.lab)
  labels_new@n.idx <- labels_new@label == labels_new@negative.lab
  labels_new@p.idx <- labels_new@label == labels_new@positive.lab

  if (verbose > 0) cat('Keeping labels of',length(labels_new@label),'sample(s).\n')

  siamcat@label <- labels_new
  return(siamcat)
}

#based on https://github.com/joey711/phyloseq/blob/master/R/show-methods.R
#' @rdname show-methods
setMethod("show", "siamcat", function(object){
  cat("siamcat-class object", fill=TRUE)
  if(!is.null(object@label)) cat(paste("label()                label:           ",sum(object@label@n.idx),
            object@label@n.lab,"and",sum(object@label@p.idx),object@label@p.lab,"samples", sep = " "), fill = TRUE)
  if(!is.null(object@modelList)){
    cat(paste("modelList()            modelList:       ",length(object@modelList@models),
              object@modelList@model.type,"models", sep = " "), fill = TRUE)
  }
  
  # print otu_table (always there).
  cat("\ncontains phyloseq-class experiment-level object @phyloseq:", fill=TRUE)
  cat(paste("phyloseq@otu_table()   OTU Table:         [ ", ntaxa(otu_table(object@phyloseq)), " taxa and ", 
            nsamples(otu_table(object@phyloseq)), " samples ]", sep = ""), fill = TRUE)	
  
  # print Sample Data if there
  if(!is.null(sample_data(object@phyloseq, FALSE))){
    cat(paste("phyloseq@sample_data() Sample Data:       [ ", dim(sample_data(object@phyloseq))[1], " samples by ", 
              dim(sample_data(object@phyloseq))[2], 
              " sample variables ]", sep = ""), fill = TRUE)
  }
  
  # print tax Tab if there	
  if(!is.null(tax_table(object@phyloseq, FALSE))){
    cat(paste("phyloseq@tax_table()   Taxonomy Table:    [ ", dim(tax_table(object@phyloseq))[1], " taxa by ", 
              dim(tax_table(object@phyloseq))[2], 
              " taxonomic ranks ]", sep = ""), fill = TRUE)
  }
  
  # print tree if there
  if(!is.null(phy_tree(object@phyloseq, FALSE))){
    cat(paste("phyloseq@phy_tree()    Phylogenetic Tree: [ ", ntaxa(phy_tree(object@phyloseq)), " tips and ", 
              phy_tree(object@phyloseq)$Nnode,
              " internal nodes ]", sep = ""),
        fill = TRUE
    ) 
  }
  
  # print refseq summary if there
  if(!is.null(refseq(object@phyloseq, FALSE))){
    cat(paste("phyloseq@refseq()      ", class(refseq(object@phyloseq))[1], ":      [ ", ntaxa(refseq(object@phyloseq)), " reference sequences ]", sep = ""), fill=TRUE)
  }
})