#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between Microbial
### Communities And host phenoTypes
### EMBL Heidelberg 2012-2018 GNU GPL 3.0

#' Reset features in siamcat@phylose@otu_table to those in siamcat@orig_feat
#' @title reset.features
#' @name reset.features
#' @description Function reset features in siamcat@phylose@otu_table to those in
#' siamcat@orig_feat in an object of class \link{siamcat-class}
#' @param siamcat an object of class \link{siamcat-class}t
#' @return A new \link{siamcat-class} object
#' @export
#' @examples
#'  data(siamcat_example)
#'  siamcat_example <- reset.features(siamcat_example)
reset.features <- function(siamcat) {
  siamcat@phyloseq@otu_table <- siamcat@orig_feat
  return(siamcat)
}