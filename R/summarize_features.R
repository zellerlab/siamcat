#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Summarize features
#'
#' @description Summarize features
#'
#' @usage summarize.features(siamcat, level = 'g__',
#'                     feature.type='original', verbose=1)
#'
#' @param siamcat object of class \link{siamcat-class}
#'
#' @param level at which level to summarize (g__ = genus)
#'
#' @param feature.type On which type of features should the function work? Can
#'   be either "original", "filtered", or "normalized". Please only change this
#'   paramter if you know what you are doing!
#'
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'     for only information about progress and success, \code{2} for normal
#'     level of information and \code{3} for full debug information,
#'     defaults to \code{1}
#'
#' @keywords internal
#'
#' @export summarize.features
#'
#' @return object of class \link{siamcat-class} with a summarized feature table
#'
#' @examples
#'
#' data("GlobalPatterns") ## phyloseq example data
#' feat <- otu_table(GlobalPatterns)[1:500,]
#' label <- create.label(meta=sample_data(GlobalPatterns),
#'     label = "SampleType",
#'     case = c("Freshwater",
#'         "Freshwater (creek)",
#'         "Ocean"))
#' # rename features
#' temp <- tax_table(GlobalPatterns)[1:500,]
#' test <- apply(temp, 1, FUN=function(vec){
#'     out <- ''
#'     for (i in seq_along(vec)){
#'         end <- ifelse(i == ncol(temp), '', ';')
#'         x <- colnames(temp)[i]
#'         x2 <- tolower(substr(x, 1, 1))
#'         out <- paste0(out, x2, '__', vec[i], end)
#'     }
#'     return(out)})
#' rownames(feat) <- test
#' # run the constructor function
#' siamcat <- siamcat(feat=feat, label=label)
#' siamcat <- summarize.features(siamcat, level='g__', verbose=3)
summarize.features <- function(siamcat, level = "g__",
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

    # TODO
    # checks and balances
    # check that all feature names are on the same taxonomic Level
    # check that desired level is in the names

    # search for bins
    tax.table <- str_split(rownames(feat),
        pattern='(\\.[a-z]__)|(^[a-z]__)',
        simplify = TRUE)[,-1]
    colnames(tax.table) <- str_extract_all(rownames(feat)[1], '[a-z]__')[[1]]

    idx <- match(level, colnames(tax.table))
    if (is.na(idx)){
        stop('Level ', level, ' not found in the feature names. Exiting...')
    }
    # rename unclassified features
    tax.table[tax.table[,idx] == '', idx] <- 'unclassified'
    bins.unique <- unique(tax.table[,idx])

    # make empty matrix
    summarized.feat <- matrix(NA, nrow=length(bins.unique), ncol=ncol(feat),
        dimnames=list(bins.unique, colnames(feat)))

    if (verbose > 0) pb <- txtProgressBar(max=length(bins.unique), style=3)
    for (x in bins.unique){
        summarized.feat[x,] <- colSums(feat[which(tax.table[,idx] == x),,
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

    if (feature.type == 'original'){
        orig_feat(siamcat) <- otu_table(summarized.feat,taxa_are_rows = TRUE)
    } else if (feature.type == 'filtered'){
        filt_feat(siamcat) <- new('filt_feat',
            filt.feat=otu_table(summarized.feat,taxa_are_rows = TRUE),
            filt.param=filt_params(siamcat))
    } else if (feature.type == 'normalized'){
        norm_feat(siamcat) <- new("norm_feat",
            norm.feat=otu_table(summarized.feat,taxa_are_rows = TRUE),
            norm.param=norm_params(siamcat))
    }

    return(siamcat)
}
