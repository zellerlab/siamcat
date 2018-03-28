#!/usr/bin/Rscript
###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# R flavor
# EMBL Heidelberg 2012-2018
# GNU GPL 3.0
###

#' @title Perform unsupervised feature filtering.
#' @description This function performs unsupervised feature filtering. Features
#'        can be filtered based on abundance or prevalence. Additionally,
#'        unmapped reads may be removed.
#' @param siamcat an object of class \link{siamcat-class}
#' @param filter.method method used for filtering the features, can be one of
#'        these: \code{c("abundance", "cum.abundance", "prevalence")},
#'        defaults to \code{"abundance"}
#' @param cutoff float, abundace or prevalence cutoff, default to \code{0.001}
#' @param recomp.prop boolean, should relative abundances be recomputed?,
#'        defaults to \code{FALSE}
#' @param rm.unmapped boolean, should unmapped reads be discarded?, defaults to
#'        \code{TRUE}
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'        for only information about progress and success, \code{2} for normal
#'        level of information and \code{3} for full debug information, defaults to \code{1}
#' @keywords SIAMCAT filter.features
#' @details This function filters the features in a \link{siamcat-class} object in a
#'        unsupervised manner.
#'
#'        The different filter methods work in the following way: \itemize{
#'        \item \code{"abundace"} remove features whose abundance is never above
#'          the threshold value in any of the samples
#'        \item \code{"cum.abundance"} remove features with very low abundance
#'          in all samples i.e. ones that are never among the most abundant
#'          entities that collectively make up (1-cutoff) of the reads in any sample
#'        \item \code{"prevalence"} remove features with low prevalence across samples
#'          i.e. ones that are 0 (undetected) in more than (1-cutoff) proportion
#'          of samples.
#'        }
#' @export
#' @return siamcat an object of class \link{siamcat-class}
#' @examples
#'  # Example dataset
#'  data(siamcat_example)
#'  # since the whole pipeline has been run in the example data, the feature were
#'  # filtered already. Remove...
#'  siamcat_example@phyloseq@otu_table <- NULL
#'
#' # Simple examples
#' siamcat.filtered <- filter.features(siamcat_example, filter.method='abundance', cutoff=1e-03)
filter.features <- function(siamcat, filter.method="abundance", cutoff=0.001,
                            recomp.prop=FALSE, rm.unmapped=TRUE, verbose=1){
  ### this statement does not have the purpose to calculate relative abundances on the fly and return them.
  ### Instead, it's purpose is to be able to calculate f.idx (specifying the indices of features which are to be kept)
  ### when feature list has already been transformed to relative abundances, but e.g. certain features have been removed manually.
  ## TODO check filter.method, add default value for cutoff, recomp.prop, and rm.unmapped?
  # Check filter methods
  if(verbose>1) cat("+ starting filter.features\n")
  s.time <- proc.time()[3]

  if (!filter.method %in% c("abundance", "cum.abundace", "prevalence")){
    stop('! Unrecognized filter.method, exiting!\n')
  }

  if (verbose > 1) cat("+++ before filtering, the data has", nrow(siamcat@phyloseq@otu_table), "features\n")
  if (recomp.prop) {
    # recompute relative abundance values (proportions)
    ra.feat <- prop.table(siamcat@phyloseq@otu_table, 2)
    siamcat@phyloseq@otu_table <- otu_table(ra.feat,taxa_are_rows = T )
  } else {
    ra.feat <- siamcat@phyloseq@otu_table
  }

  ### apply filters
  if(verbose>2) cat("+++ applying",filter.method,"filter\n")
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
  }


  if(verbose>2) cat("+++ checking for unmapped reads\n")
  ### postprocessing and output generation
  if (rm.unmapped) {
    # remove 'unmapped' feature
    names.unmapped <- c('UNMAPPED', '-1', 'X.1', 'unmapped', 'UNCLASSIFIED', 'unclassified', 'UNASSIGNED', 'unassigned')
    unm.idx <- rownames(siamcat@phyloseq@otu_table) %in% names.unmapped
    if (any(unm.idx)) {
      f.idx <- f.idx[-which(f.idx%in%unm.idx)]
      if (verbose > 2) cat("+++ removing ", rownames(siamcat@phyloseq@otu_table)[unm.idx], ' as unmapped reads\n')
      if (verbose > 1) cat('+++ removed ', sum(unm.idx), ' features corresponding to UNMAPPED reads\n', sep='')
    } else {
      if (verbose > 1) cat('+++ tried to remove unmapped reads, but could not find them. Continue anyway.\n')
    }
  }
  if(verbose>2) cat("+++ applying prune_taxa\n")
  if (verbose > 1) cat('+++ removed ', nrow(siamcat@phyloseq@otu_table)-length(f.idx)-sum(unm.idx), ' features whose values did not exceed ',
    cutoff, ' in any sample (retaining ', length(f.idx), ')\n', sep='')
  f.names <- rownames(siamcat@phyloseq@otu_table)[f.idx]
  siamcat@phyloseq <- prune_taxa(x = siamcat@phyloseq, taxa = f.names)
  e.time <- proc.time()[3]
  if(verbose>1) cat("+ finished filter.features in",e.time-s.time,"s\n")
  if(verbose==1)cat("Features successfully filtered\n")
  return(siamcat)
}
