#!/usr/bin/Rscript
###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# R flavor
# EMBL Heidelberg 2012-2018
# GNU GPL 3.0
###

#' The S4 class for storing taxa-abundance information and models.
#' @name siamcat-class
#' @rdname siamcat-class
#' @exportClass siamcat
setClass("siamcat", representation(modelList = "list", phyloseq = "phyloseq"))

#' Build siamcat-class objects from their components.
#' @name siamcat
#' @export

siamcat <- function(...){
  arglist   <- list(...)
  modelList <- list(NULL)
  
  sc        <- new("siamcat", name = "Hadley", phyloseq = ps)
  return(sc)
}