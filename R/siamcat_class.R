#!/usr/bin/Rscript
###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# R flavor
# EMBL Heidelberg 2012-2018
# GNU GPL 3.0
###

#' @name siamcat-class
#' @rdname siamcat-class
#' @exportClass siamcat
setClass("siamcat", representation(modelList = "list", phyloseq = "phyloseq"))