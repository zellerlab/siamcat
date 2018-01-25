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
filter.label <- function(siamcat,ids){
  labels_new      <- new("label", label = siamcat@label@label[ids],
                         header = siamcat@label@header,
                         info = siamcat@label@info,
                         positive.lab = siamcat@label@positive.lab,
                         negative.lab = siamcat@label@negative.lab,
                         n.lab = siamcat@label@n.lab,
                         p.lab = siamcat@label@p.lab)
  labels_new@n.idx <- labels_new@label == labels_new@negative.lab
  labels_new@p.idx <- labels_new@label == labels_new@positive.lab
  
  cat('Keeping labels of',length(labels_new@label),'sample(s).\n')
  
  siamcat@label <- labels_new
}