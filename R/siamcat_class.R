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

#' The S4 class for storing data splits
#' @name dataSplit-class
#' @rdname dataSplit-class
#' @exportClass dataSplit
setClass("dataSplit", representation(training.folds = "list", 
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
setClass("siamcat", representation(modelList = "modelList", phyloseq = "phyloseq", orig_feat="otu_table", evalData = "list", 
                                   label="label", norm.param="list", dataSplit="dataSplit", predMatrix="matrix"))

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
    print(names(arglistphyloseq))
    arglist$phyloseq <- do.call("new", c(list(Class="phyloseq"), arglistphyloseq))
  }
  arglist     <- arglist[sapply(names(arglist), is.component.class, "siamcat")]
  sc          <- do.call("new", c(list(Class="siamcat"), arglist))
    #new("siamcat", modelList = arglist$modelList, phyloseq = ps, predMatrix = arglist$predMatrix,
     #                  orig_feat = arglist$otu_table, label = arglist$label, norm.param = arglist$norm.param)
  return(sc)
}

# source: https://github.com/joey711/phyloseq/blob/master/R/phyloseq-class.R
#' Show the component objects classes and slot names.
#' @usage get.component.classes()
#' @keywords internal
get.component.classes <- function(class){
  # define classes vector
  # the names of component.classes needs to be the slot names to match getSlots / splat
  component.classes.siamcat <- c("modelList", "orig_feat", "label", "norm.param", "dataSplit","phyloseq")	#slot names
  names(component.classes.siamcat) <- c("modelList", "orig_feat", "label","norm.param", "dataSplit", "phyloseq") #class names
  
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
