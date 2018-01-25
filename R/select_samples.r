###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# RScript flavor
#
# written by Georg Zeller
# with additions by Nicolai Karcher and Konrad Zych
# EMBL Heidelberg 2012-2017
#
# version 0.2.0
# file last updated: 26.06.2017
# GNU GPL 3.0
###

#' @title Select samples based on siamcat@phyloseq@sam_datadata
#' @description This functions filters features, labels, and siamcat@phyloseq@sam_datadata based on a specific column in the siamcat@phyloseq@sam_datadata.
#' Provided with a column in the siamcat@phyloseq@sam_datadata and a range or a set of allowed values, the function will return the labels, siamcat@phyloseq@sam_datadata, and features for the samples matching the allowed range or set.
#' @param siamcat an object of class \link{siamcat}
#' @param filter name of the column from siamcat@phyloseq@sam_datadata on which the selection should be done
#' @param allowed.set a vector of allowed values
#' @param allowed.range a range of allowed values
#' @keywords SIAMCAT select.samples
#' @export
#' @return siamcat an object of class \link{siamcat}
#' 
select.samples  <- function(siamcat, filter, allowed.set = NULL, allowed.range = NULL){

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
    stop('Neither allowed.range  nor allowed.set (or both at the same time) have been provided, exiting!\n')
  } else {
    if (!is.null(allowed.range )) {
      cat('allowed.range  = [', paste(allowed.range , collapse=','), ']\n', sep='')
    } else {
      cat('allowed.set = {', paste(allowed.set, collapse=','), '}\n', sep='')
    }
  }
  cat('\n')

  if(!filter %in% colnames(siamcat@phyloseq@sam_data)) stop("The filter name is not present in colnames of the siamcat@phyloseq@sam_datadata. Stopping.\n")
  filter.var    <- siamcat@phyloseq@sam_data[,filter]

  if (!is.null(allowed.range )) {
    s.idx <- !is.na(filter.var) & filter.var >= allowed.range [1] & filter.var <= allowed.range [2]
    cat('Removed ', sum(!s.idx), ' samples with ', filter,
        ' not in [', paste(allowed.range , collapse=', '), '] (retaining ', sum(s.idx), ')\n', sep='')

  } else {
    s.idx <- !is.na(filter.var) & filter.var %in% allowed.set
    cat('Removed ', sum(!s.idx), ' samples with ', filter,
        ' not in {', paste(allowed.set, collapse=', '), '} (retaining ', sum(s.idx), ')\n', sep='')
  }

  s.names <- rownames(siamcat@phyloseq@sam_data)[s.idx]
  siamcat@phyloseq <- prune_samples(x = siamcat@phyloseq, samples = s.names)
  siamcat          <- filter.label(siamcat,)

  return(siamcat)
}
