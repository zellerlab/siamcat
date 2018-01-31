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
                                   label="label", norm.param="list"))

#' Build siamcat-class objects from their components.
#' @name siamcat
#' @export
siamcat <- function(...){
  arglist   <- list(...)
  
  # Remove names from arglist. Will replace them based on their class
  names(arglist) <- NULL
  
  # ignore all but component data classes.
  component_classes <- get.component.classes()
  
  for(argNr in 1:length(arglist)){
    classOfArg <- class(arglist[[argNr]])[1]
    if(classOfArg%in%names(component_classes)){
      names(arglist)[argNr] <- classOfArg
    }
  }
  
  if(is.null(arglist$modelList))  arglist$modelList <- new("modelList",models=list(NULL), model.type="empty")
  if(is.null(arglist$norm.param)) arglist$norm.param <- list(NULL)
  
  if(is.null(arglist$phyloseq)) ps <- phyloseq(arglist$otu_table, arglist$sample_data, arglist$phylo, 
                 arglist$taxonomyTable, arglist$XStringSet) 
  
  sc <- new("siamcat", modelList = arglist$modelList, phyloseq = ps,
                       orig_feat = arglist$otu_table, label = arglist$label, norm.param = arglist$norm.param)
  return(sc)
}

# source: https://github.com/joey711/phyloseq/blob/master/R/phyloseq-class.R
#' Show the component objects classes and slot names.
#' @usage get.component.classes()
#' @keywords internal
get.component.classes <- function(){
  # define classes vector
  # the names of component.classes needs to be the slot names to match getSlots / splat
  component.classes <- c("otu_table", "sam_data", "phy_tree", "tax_table", "refseq", "modelList",
                                "orig_feat", "label", "list", "phyloseq")	#slot names
  names(component.classes) <- c("otu_table", "sample_data", "phylo", "taxonomyTable", "XStringSet", "modelList",
                                 "orig_feat", "label","norm.param", "phyloseq") #class names
  return(component.classes)
}
