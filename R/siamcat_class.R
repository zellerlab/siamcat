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


#' The S4 class for storing taxa-abundance information and models.
#' @name siamcat-class
#' @rdname siamcat-class
#' @exportClass siamcat
setClass("siamcat", representation(modelList = "modelList", phyloseq = "phyloseq"))

#' Build siamcat-class objects from their components.
#' @name siamcat
#' @export

siamcat <- function(...){
  arglist   <- list(...)
  ps <- phyloseq(arglist$otu_table,arglist$sample_data,arglist$tax_table,arglist$otu_table) 
  sc <- new("siamcat", modelList = arglist$modelList, phyloseq = ps)
  return(sc)
}