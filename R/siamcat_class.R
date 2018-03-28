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
#' @slot phyloseq object of class \link[phyloseq]{phyloseq-class}
#' @slot label an object of class \link{label-class}
#' @slot orig_feat an object of class \link[phyloseq]{otu_table}
#' @slot data_split a list
#' @slot norm_param a list
#' @slot model_list an object of class \link{model_list-class}
#' @slot eval_data list containing \itemize{
#'  \item \code{$roc.average} average ROC-curve across repeats or a single ROC-curve on complete dataset;
#'  \item \code{$auc.average} AUC value for the average ROC-curve;
#'  \item \code{$ev.list} list of \code{length(num.folds)}, containing for different decision thresholds the number of false positives, false negatives, true negatives, and true positives;
#'  \item \code{$pr.list} list of \code{length(num.folds)}, containing the positive predictive value (precision) and true positive rate (recall) values used to plot the PR curves;
#' }. If \code{prediction} had more than one column, i.e. if the models has been trained with several repeats, the function will additonally return \itemize{
#'  \item \code{$roc.all} list of roc objects (see \link[pROC]{roc}) for every repeat;
#'  \item \code{$aucspr} vector of AUC values for the PR curves for every repeat;
#'  \item \code{$auc.all} vector of AUC values for the ROC curves for every repeat
#'}
#' @slot pred_matrix a matrix with predictions made by \link{make.predictions} function
#' @exportClass siamcat
setClass("siamcat", representation(model_list = "model_list", phyloseq = "phyloseq", orig_feat="otu_table", eval_data = "list",
                                   label="label", norm_param="list", data_split="data_split", pred_matrix="matrix"))

#' Build siamcat-class objects from their components.
#' @title construct.siamcat
#' @description Function to construct an object of class \link{siamcat-class}
#' @param ... list of arguments needed in order to construct a SIAMCAT object
#' @export
construct.siamcat <- function(...){
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
