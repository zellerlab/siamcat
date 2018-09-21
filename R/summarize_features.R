#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Summarize features
#'
#' @description Summarize features
#'
#' @param siamcat object of class \link{siamcat-class}
#'
#' @param siamcat at which level to summarize (g__ = genus)
#'
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'     for only information about progress and success, \code{2} for normal
#'     level of information and \code{3} for full debug information,
#'     defaults to \code{1}
#'
#' @keywords internal
#'
#' @export
#'
#' @return object of class \link{siamcat-class} with a summarized feature table
#'
#' @examples
#'
#'     data(siamcat_example)
#'     # simple working example
#'     siamcat_evaluated <- summarize.features(siamcat_example)
#'
summarize.features <- function(siamcat, level = "g__",
                               keep_original_names=TRUE,
                               feature.type='original', verbose=1){

    if (verbose > 1) message("+ starting summarize.features")

    s.time <- proc.time()[3]

    if (verbose > 2) message("+++ summarizing on level: ",level)

    if (!feature.type %in% c('original', 'filtered', 'normalized')){
        stop("Unrecognised feature type, exiting...\n")
    }
    # get the right features
    if (feature.type == 'original'){
        feat <- get.orig_feat.matrix(siamcat)
    } else if (feature.type == 'filtered'){
        if (is.null(filt_feat(siamcat, verbose=0))){
            stop('Features have not yet been filtered, exiting...\n')
        }
        feat <- get.filt_feat.matrix(siamcat)
    } else if (feature.type == 'normalized'){
        if (is.null(norm_feat(siamcat, verbose=0))){
            stop('Features have not yet been normalized, exiting...\n')
        }
        feat <- get.norm_feat.matrix(siamcat)
    }
    # make sure that seperating characters are dots
    rownames(feat) <- make.names(rownames(feat))

    # search for bins
    prefix <- ifelse(keep_original_names==TRUE, '^.*', '')
    bins <- str_extract_all(rownames(feat),
        pattern=paste0(prefix, level, '[A-Za-z0-9]+\\.'))
    bins.unique <- unique(unlist(bins))

    # make empty matrix
    summarized.feat <- matrix(NA, nrow=length(bins.unique), ncol=ncol(feat),
                              dimnames=list(bins.unique, colnames(feat)))

    if (verbose > 0) pb <- txtProgressBar(max=length(bins.unique), style=3)
    for (x in bins.unique){
        summarized.feat[x,] <- colSums(feat[grep(x, rownames(feat)),,
            drop=FALSE])
        if (verbose > 0) setTxtProgressBar(pb, pb$getVal()+1)
    }

    if (verbose > 2) message("\n+++ summarized features table contains: ",
        nrow(summarized.feat)," features\n")

    e.time <- proc.time()[3]
    if (verbose > 1)
    message(paste("+ finished summarize.features in",
        formatC(e.time - s.time, digits = 3), "s" ))

    if (verbose == 1)
        message("\nSummarized features successfully.")
    features(siamcat, feature.type) <-
        otu_table(summarized.feat,taxa_are_rows = TRUE)

    return(siamcat)
}
