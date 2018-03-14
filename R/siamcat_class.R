#!/usr/bin/Rscript
###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# R flavor
# EMBL Heidelberg 2012-2018
# GNU GPL 3.0
###

#' The S4 class for storing models.
#' @name model_list-class
#' @rdname model_list-class
#' @exportClass model_list
setClass("model_list", representation(models = "list", model.type = "character"))

#' The S4 class for storing data splits
#' @name data_split-class
#' @rdname data_split-class
#' @exportClass data_split
setClass("data_split", representation(training.folds = "list", 
                                     test.folds = "list", 
                                     num.resample = "numeric", 
                                     num.folds = "numeric"))

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
setClass("siamcat", representation(model_list = "model_list", phyloseq = "phyloseq", orig_feat="otu_table", eval_data = "list", 
                                   label="label", norm_param="list", data_split="data_split", predMatrix="matrix"))

#' Build siamcat-class objects from their components.
#' @name siamcat
#' @export
siamcat <- function(...){
  arglist   <- list(...)
  
  # Remove names from arglist. Will replace them based on their class
  names(arglist) <- NULL
  
  # ignore all but component data classes.
  component_classes <- get.component.classes("both")
  
  for(argNr in 1:length(arglist)){
    classOfArg <- class(arglist[[argNr]])[1]
    if(classOfArg%in%names(component_classes)){
      names(arglist)[argNr] <- component_classes[classOfArg]
    }
  }
  
  if(is.null(arglist$phyloseq)){
    arglistphyloseq <- arglist[sapply(names(arglist), is.component.class, "phyloseq")]
    arglist$phyloseq <- do.call("new", c(list(Class="phyloseq"), arglistphyloseq))
  }
  arglist     <- arglist[sapply(names(arglist), is.component.class, "siamcat")]
  sc          <- do.call("new", c(list(Class="siamcat"), arglist))
  return(sc)
}

# source: https://github.com/joey711/phyloseq/blob/master/R/phyloseq-class.R
#' Show the component objects classes and slot names.
#' @usage get.component.classes()
#' @keywords internal
get.component.classes <- function(class){
  # define classes vector
  # the names of component.classes needs to be the slot names to match getSlots / splat
  component.classes.siamcat <- c("model_list", "orig_feat", "label", "norm_param", "data_split","phyloseq")	#slot names
  names(component.classes.siamcat) <- c("model_list", "orig_feat", "label","norm_param", "data_split", "phyloseq") #class names
  
  component.classes.phyloseq <- c("otu_table", "sam_data", "phy_tree", "tax_table", "refseq")	#slot names
  names(component.classes.phyloseq) <- c("otu_table", "sample_data", "phylo", "taxonomyTable", "XStringSet") #class names
  
  if(class=="siamcat"){
    return(component.classes.siamcat)
  }else if(class=="phyloseq"){
    return(component.classes.phyloseq)
  }else if(class=="both"){
    return(c(component.classes.siamcat,component.classes.phyloseq))
  }
}

# Returns TRUE if x is a component class, FALSE otherwise.
# This shows up over and over again in data infrastructure
#' @keywords internal
is.component.class = function(x,class){
  x%in%get.component.classes(class)
}
