#!/usr/bin/Rscript
###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# R flavor
# EMBL Heidelberg 2012-2018
# GNU GPL 3.0
###

#' @title Select samples based on metadata
#' @description This functions filters features, labels, and metadata based on
#'        a specific column in the metadata. Provided with a column in the
#'        metadata and a range or a set of allowed values, the function will
#"        filter the \link{siamcat} object accordingly.
#' @param siamcat an object of class \link{siamcat}
#' @param filter string, name of the meta variable on which the selection
#'        should be done
#' @param allowed.set a vector of allowed values
#' @param allowed.range a range of allowed values
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'        for only information about progress and success, \code{2} for normal
#'        level of information and \code{3} for full debug information, defaults to \code{1}
#' @keywords SIAMCAT select.samples
#' @export
#' @return siamcat an object of class \link{siamcat}
#' @examples
#' siamcat.sel <- select.samples(siamcat, 'Age', allowed.range=c(20,80))
#' siamcat.sel <- select.samples(siamcat, 'Gender', allowed.set=c('M', 'F'), verbose=2)
select.samples  <- function(siamcat, filter, allowed.set = NULL, allowed.range = NULL, verbose=1){
  if(verbose>1) cat("+ starting select.samples\n")
  s.time <- proc.time()[3]

  if(verbose>2) cat("+++ checking allowed values\n")
  if(!is.null(allowed.range )) {
    allowed.range  <- gsub('\\[|\\]','', allowed.range)
    allowed.range  <- as.numeric(unlist(strsplit(allowed.range,',')))
    stopifnot(length(allowed.range) == 2)
    stopifnot(allowed.range [1] <= allowed.range[2])
  }

  if (!is.null(allowed.set)) {
    # parse set
    allowed.set  <- gsub('\\{|\\}','', allowed.set)
    allowed.set  <- as.numeric(unlist(strsplit(allowed.set,',')))
    allowed.set  <- sort(unique(allowed.set))
  }

  if (!xor(is.null(allowed.range ), is.null(allowed.set))) {
    stop('Neither allowed.range nor allowed.set (or both at the same time) have been provided, exiting!\n')
  } else {
    if (!is.null(allowed.range )) {
      if (verbose > 2) cat('+++ allowed.range  = [', paste(allowed.range , collapse=','), ']\n', sep='')
    } else {
      if (verbose > 2) cat('+++ allowed.set = {', paste(allowed.set, collapse=','), '}\n', sep='')
    }
  }

  if(!filter %in% colnames(siamcat@phyloseq@sam_data)) stop("! The filter name is not present in colnames of the siamcat@phyloseq@sam_data. Stopping.\n")
  filter.var    <- siamcat@phyloseq@sam_data[,filter]

  if (!is.null(allowed.range)) {
    s.idx <- !is.na(filter.var) & filter.var >= allowed.range [1] & filter.var <= allowed.range [2]
    if (verbose > 1) {
      cat('+++ removed ', sum(!s.idx), ' samples with ', filter, ' not in [',
      paste(allowed.range , collapse=', '), '] (retaining ', sum(s.idx), ')\n',
      sep='')
    }

  } else {
    s.idx <- !is.na(filter.var) & filter.var %in% allowed.set
    if (verbose > 1){
      cat('+++ removed ', sum(!s.idx), ' samples with ', filter, ' not in {',
      paste(allowed.set, collapse=', '), '} (retaining ', sum(s.idx), ')\n',
      sep='')
    }
  }

  s.names <- rownames(siamcat@phyloseq@sam_data)[s.idx]
  siamcat@phyloseq <- prune_samples(x = siamcat@phyloseq, samples = s.names)
  siamcat          <- filter.label(siamcat, verbose=verbose)
  e.time <- proc.time()[3]
  if(verbose>1) cat("+ finished select.samples in",e.time-s.time,"s\n")
  if(verbose==1)cat("Selecting samples finished\n")
  return(siamcat)
}
