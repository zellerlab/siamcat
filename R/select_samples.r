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

#' @title Select samples based on metadata
#' @description This functions filters features, labels, and metadata based on a specific column in the metadata.
#' Provided with a column in the metadata and a range or a set of allowed values, the function will return the labels, metadata, and features for the samples matching the allowed range or set.
#' @param meta metadata object
#' @param feat features object
#' @param label labels object
#' @param filter name of the column from metadata on which the selection should be done
#' @param allowed.set a vector of allowed values
#' @param allowed.range a range of allowed values
#' @keywords SIAMCAT select.samples
#' @export
#' @return list containing values after selection: \itemize{
#'  \item \code{$feat} = features;
#'  \item \code{$label} = labels;
#'  \item \code{$meta} = metadata
#' }
select.samples  <- function(meta, feat, label, filter, allowed.set = NULL, allowed.range = NULL){
  # TODO at the moment, this functions takes label$label, not an label-object. Should be corrected
  # parse interval
  if(!is.null(allowed.range )) {
    allowed.range  <- gsub('\\[|\\]','', allowed.range)
    allowed.range  <- as.numeric(unlist(strsplit(allowed.range,',')))
    stopifnot(length(allowed.range) == 2)
    stopifnot(allowed.range [1] <= allowed.range[2])
  }
  #allowed.set  <- opt$allowed.set # TODO should be removed, i guess?
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
  # TODO throws an error when called in Rstudio, should be removed i guess
  # if (substr(source.dir, nchar(source.dir), nchar(source.dir)) != '/') {
  #   source.dir   <- paste(source.dir, '/', sep='')
  # }
  if(!filter %in% colnames(meta)) stop("The filter name is not present in colnames of the metadata. Stopping.\n")
  filter.var    <- meta[,filter]

  if (!is.null(allowed.range )) {
    s.idx <- !is.na(filter.var) & filter.var >= allowed.range [1] & filter.var <= allowed.range [2]
    cat('Removed ', sum(!s.idx), ' samples with ', filter,
        ' not in [', paste(allowed.range , collapse=', '), '] (retaining ', sum(s.idx), ')\n', sep='')
  } else {
    s.idx <- !is.na(filter.var) & filter.var %in% allowed.set
    cat('Removed ', sum(!s.idx), ' samples with ', filter,
        ' not in {', paste(allowed.set, collapse=', '), '} (retaining ', sum(s.idx), ')\n', sep='')
  }
  results         <- list(NULL) # why?
  results$feat    <- feat[,s.idx]
  results$label   <- label[s.idx]
  results$meta    <- meta[s.idx,]
  # why invisible return?
  invisible(results)
}
