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

#' @title Perform unsupervised feature filtering.
#' @description This function may convert absolute abundances into relative abundances and then performs unsupervised feature filtering. Features can be filtered based on abundance or prevalence. Additionally, unmapped reads may be removed.
#' @param feat feature object
#' @param filter.method method used for filtering the features, can be one of these: \code{c("abundance", "cum.abundance", "prevalence")}
#' @param cutoff float, abundace or prevalence cutoff
#' @param recomp.prop boolean, should absolute abundances be converted into relative abundances?
#' @param rm.unmapped boolean, should unmapped reads be discarded?
#' @details Currently, there are three filtering methods implemented:
#' \itemize{
#'  \item \code{"abundance"} remove features for which the abundance is never above the threshold value (e.g. 0.5\%) in any of the samples
#'  \item \code{"cum.abundance"} remove features with very low abundance in all samples, i.e. features that are never among the most abundant entities that collectively make up \code{(1-cutoff)} of the reads in any samples
#'  \item \code{"prevalence"} remove features with low prevalence across samples, i.e. features that are undetected in more than \code{(1-cutoff)} proportion of samples
#' }
#' @keywords SIAMCAT filter.feat
#' @export
#' @return Returns the filtered feature matrix
filter.feat <- function(feat, filter.method=abundance, cutoff=0.001, recomp.prop=FALSE, rm.unmapped=TRUE){
  ### this statement does not have the purpose to calculate relative abundances on the fly and return them.
  ### Instead, it's purpose is to be able to calculate f.idx (specifying the indices of features which are to be kept)
  ### when feature list has already been transformed to relative abundances, but e.g. certain features have been removed manually.
  ## TODO check filter.method, add default value for cutoff, recomp.prop, and rm.unmapped?
  if (recomp.prop) {
    # recompute relative abundance values (proportions)
    ra.feat <- prop.table(feat, 2)
  } else {
    ra.feat <- feat
  }

  ### apply filters
  if (filter.method == 'abundance') {
    # remove features whose abundance is never above the threshold value (e.g. 0.5%) in any of the samples
    f.max <- apply(ra.feat, 1, max)
    f.idx <- which(f.max >= cutoff)
  } else if (filter.method == 'cum.abundance') {
    # remove features with very low abundance in all samples i.e. ones that are never among the most abundant
    # entities that collectively make up (1-cutoff) of the reads in any sample
    f.idx <- vector('numeric', 0)
    # sort features per sample and apply cumsum to identify how many collectively have weight K
    for (s in 1:ncol(ra.feat)) {
      srt   <- sort(ra.feat[,s], index.return=TRUE)
      cs    <- cumsum(srt$x)
      m     <- max(which(cs < cutoff))
      f.idx <- union(f.idx, srt$ix[-(1:m)])
    }
    # an index of those features that collectively make up more than 1-K of the read mass in any sample
    f.idx <- sort(f.idx)
  } else if (filter.method == 'prevalence') {
    # remove features with low prevalence across samples
    # i.e. ones that are 0 (undetected) in more than (1-cutoff) proportion of samples
    f.idx <- which(rowSums(ra.feat > 0) / ncol(ra.feat) > cutoff)
  } else {
    stop('unrecognized filter.method, exiting!\n')
  }

  cat('Removed ', nrow(feat)-length(f.idx), ' features whose values did not exceed ', cutoff,
      ' in any sample (retaining ', length(f.idx), ')\n', sep='')
  feat <- feat[f.idx,]

  ### postprocessing and output generation
  if (rm.unmapped) {
    # remove 'unmapped' feature
    unm.idx <- rownames(feat) == 'UNMAPPED' | rownames(feat) == 'unmapped' | rownames(feat) == '-1' | rownames(feat) == 'X.1' | rownames(feat) == 'UNCLASSIFIED' | rownames(feat) == 'unclassified' | rownames(feat) == 'UNASSIGNED' | rownames(feat) == 'unassigned'
    if (any(unm.idx)) {
      feat <- feat[!unm.idx,]
      cat('Removed ', sum(unm.idx), ' features corresponding to UNMAPPED reads',
          ' (retaining ', nrow(feat), ')\n', sep='')
    } else {
      cat('tried to remove unmapped reads, but could not find them. Continue anyway.')
    }
  }
  return(feat)
}
