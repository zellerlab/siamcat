#!/usr/bin/Rscript
###
# SIAMCAT -  Statistical Inference of Associations between
#   Microbial Communities And host phenoTypes
# R flavor
# EMBL Heidelberg 2012-2018
# GNU GPL 3.0
###

#' @title Add metadata as predictors
#' @description This function adds metadata to the feature matrix to be
#'        later used as predictors
#' @param siamcat object of class \link{siamcat-class}
#' @param pred.names vector of names of the variables within the metadata to be
#'        added to the feature matrix as predictors
#' @param std.meta boolean, should added metadata features be standardized?,
#'        defaults to \code{TRUE}
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'        for only information about progress and success, \code{2} for normal
#'        level of information and \code{3} for full debug information,
#'        defaults to \code{1}
#' @keywords SIAMCAT add.meta.pred
#' @export
#' @return an object of class \link{siamcat-class} with metadata added to the
#'        features
#' @examples
#'  data(siamcat_example)
#'  # Add the Age of the patients as potential predictor
#'  siamcat_age_added <- add.meta.pred(siamcat_example, pred.names=c('age'))
#'
#'  # Add Age, BMI, and Gender as potential predictors
#'  # Additionally, prevent standardization of the added features
#'  siamcat_meta_added <- add.meta.pred(siamcat_example, pred.names=c('age',
#'      'bmi', 'gender'), std.meta=FALSE)
add.meta.pred <- function(siamcat, pred.names=NULL, std.meta=TRUE, verbose=1){
  if(verbose>1) cat("+ starting add.meta.pred\n")
  s.time <- proc.time()[3]
  ### add metadata as predictors to the feature matrix
  cnt <- 0

  if (pred.names != '' && !is.null(pred.names)) {
    if(verbose>2) cat("+ starting to add metadata predictors\n")
    for (p in pred.names) {
      if(verbose>2) cat("+++ adding metadata predictor:",p,"\n")
      if(!p%in%colnames(siamcat@phyloseq@sam_data)) {
        stop("There is no metadata variable called ",p,"\n")
      }
      idx <- which(colnames(siamcat@phyloseq@sam_data) == p)
      if(length(idx) != 1) stop(p, "matches multiple columns in the metada\n")

      m   <-  unlist(siamcat@phyloseq@sam_data[,idx])

      if (!all(is.finite(m))) {
        na.cnt <- sum(!is.finite(m))
        if (verbose > 1) {
          cat('++++ filling in', na.cnt, 'missing values by mean imputation\n')
        }
        mn     <- mean(m, na.rm=TRUE)
        m[!is.finite(m)] <- mn
      }

      if (std.meta) {
        if (verbose > 1) cat('++++ standardizing metadata feature', p, '\n')
        m.mean <- mean(m, na.rm = TRUE)
        m.sd   <- stats::sd(m, na.rm = TRUE)
        stopifnot(!m.sd == 0)
        m      <- (m - m.mean)/m.sd
      }

      siamcat@phyloseq@otu_table <- phyloseq::otu_table(
        rbind(siamcat@phyloseq@otu_table, m),taxa_are_rows=TRUE)
      rownames(siamcat@phyloseq@otu_table)[nrow(
        siamcat@phyloseq@otu_table)] <- paste('META_', toupper(p), sep='')
      cnt <- cnt + 1
    }
      if (verbose > 1) cat('+++ added', cnt, 'meta-variables as predictor',
        'to the feature matrix\n')
  } else {
      if (verbose > 0) cat('+++ Not adding any of the meta-variables as',
      'predictor to the feature matrix\n')
  }
  stopifnot(all(!is.na(siamcat@phyloseq@otu_table)))
  e.time <- proc.time()[3]
  if(verbose>1) cat("+ finished add.meta.pred in",e.time-s.time,"s\n")
  if(verbose==1)cat("Adding metadata as predictor finished\n")
  return(siamcat)
}
