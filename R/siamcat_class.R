#!/usr/bin/Rscript
###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# R flavor
# EMBL Heidelberg 2012-2018
# GNU GPL 3.0
###

#' The S4 class for storing models.
#' @name modelList-class
#' @rdname modelList-class
#' @exportClass modelList
setClass("modelList", representation(models = "list", model.type = "character"))

#' The S4 class for storing label info.
#' @name label-class
#' @rdname label-class
#' @exportClass label
setClass("label", representation(label = "vector", header = "character",
                                 info="list", positive.lab="numeric",
                                 negative.lab="numeric", 
                                 n.idx="vector", p.idx="vector",
                                 n.lab="character", p.lab="character"))

#' The S4 class for storing taxa-abundance information and models.
#' @name siamcat-class
#' @rdname siamcat-class
#' @exportClass siamcat
setClass("siamcat", representation(modelList = "modelList", phyloseq = "phyloseq", orig_feat="otu_table", 
                                   label="label"))

#' Build siamcat-class objects from their components.
#' @name siamcat
#' @export
siamcat <- function(...){
  arglist   <- list(...)
  
  # Remove names from arglist. Will replace them based on their class
  names(arglist) <- NULL
  
  # ignore all but component data classes.
  arglist <- arglist[sapply(arglist, is.component.class)]
  
  for(argNr in 1:length())
  
  ps <- phyloseq(arglist$otu_table, arglist$sample_data, arglist$phylo, 
                 arglist$taxonomyTable, arglist$XStringSet) 
  
  sc <- new("siamcat", modelList = arglist$modelList, phyloseq = ps,
                       orig_feat = arglist$otu_table, label = arglist$label)
  return(sc)
}

# source: https://github.com/joey711/phyloseq/blob/master/R/phyloseq-class.R
#' Show the component objects classes and slot names.
#' @usage get.component.classes()
#' @keywords internal
get.component.classes <- function(){
  # define classes vector
  component.classes <- c("otu_table", "sample_data", "phylo", "taxonomyTable", "XStringSet", "modelList",
                         "orig_feat", "label")
  # the names of component.classes needs to be the slot names to match getSlots / splat
  names(component.classes) <- c("otu_table", "sam_data", "phy_tree", "tax_table", "refseq", "modelList",
                                "orig_feat", "label")	
  return(component.classes)
}

# source: https://github.com/joey711/phyloseq/blob/master/R/phyloseq-class.R
# Returns TRUE if x is a component class, FALSE otherwise.
# This shows up over and over again in data infrastructure
#' @keywords internal
is.component.class = function(x){
  inherits(x, get.component.classes())
}
##