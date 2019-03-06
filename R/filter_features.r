#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Perform unsupervised feature filtering.
#'
#' @description This function performs unsupervised feature filtering. Features
#'     can be filtered based on abundance, prevalence, or on variance.
#'     Additionally, unmapped reads may be removed.
#'
#' @usage filter.features(siamcat, filter.method = "abundance",
#'     cutoff = 0.001, rm.unmapped = TRUE,
#'     feature.type='original', verbose = 1)
#'
#' @param siamcat an object of class \link{siamcat-class}
#'
#' @param filter.method string, method used for filtering the features, can be
#'     one of these: \code{c('abundance', 'cum.abundance', 'prevalence',
#'     'variance')}, defaults to \code{'abundance'}
#'
#' @param cutoff float, abundace, prevalence, or variance cutoff, defaults
#'     to \code{0.001} (see Details below)
#'
#' @param rm.unmapped boolean, should unmapped reads be discarded?, defaults to
#'     \code{TRUE}
#'
#' @param feature.type string, on which type of features should the function
#'     work? Can be either "original", "filtered", or "normalized". Please
#'     only change this paramter if you know what you are doing!
#'
#' @param verbose integer, control output: \code{0} for no output at all,
#'     \code{1} for only information about progress and success, \code{2} for
#'     normal level of information and \code{3} for full debug information,
#'     defaults to \code{1}
#'
#' @keywords SIAMCAT filter.features
#'
#' @details This function filters the features in a \link{siamcat-class}
#'     object in a unsupervised manner.
#'
#'     The different filter methods work in the following way: \itemize{
#'     \item \code{'abundace'} remove features whose maximum abundance is never
#'     above the threshold value in any of the samples
#'     \item \code{'cum.abundance'} remove features with very low abundance
#'     in all samples, i.e. those that are never among the most abundant
#'     entities that collectively make up (1-cutoff) of the reads in
#'     any sample
#'     \item \code{'prevalence'} remove features with low prevalence across
#'     samples, i.e. those that are undetected (relative abundance of 0)
#'     in more than \code{1 - cutoff} percent of samples.
#'     \item \code{'variance'} remove features with low variance across
#'     samples, i.e. those that have a variance lower than \code{cutoff}
#'     }
#'
#'     Features can also be filtered repeatedly with different methods, e.g.
#'     first using the maximum abundance filtering and then using prevalence
#'     filtering.
#'     However, if a filtering method has already been applied to the dataset,
#'     SIAMCAT will default back on the original features for filtering.
#' @export filter.features
#'
#' @return siamcat an object of class \link{siamcat-class}
#'
#' @examples
#' # Example dataset
#' data(siamcat_example)
#'
#' # Simple examples
#' siamcat_filtered <- filter.features(siamcat_example,
#'     filter.method='abundance',
#'     cutoff=1e-03)
#'
filter.features <- function(siamcat,
    filter.method = "abundance",
    cutoff = 0.001,
    rm.unmapped = TRUE,
    feature.type='original',
    verbose = 1) {

    if (verbose > 1) message("+ starting filter.features")
    s.time <- proc.time()[3]

    # checks
    if (!filter.method %in% c("abundance", "cum.abundance",
                              "prevalence", "variance")) {
        stop("Unrecognized filter.method, exiting!\n")
    }
    if (!feature.type %in% c('original', 'filtered', 'normalized')){
        stop("Unrecognised feature type, exiting...\n")
    }
    if (!is.logical(rm.unmapped)){
        stop("rm.unmapped should be logical, exiting...\n")
    }

    # get the right features
    if (feature.type=='original'){
        feat <- get.orig_feat.matrix(siamcat)
        param.set <- list(list(filter.method=filter.method,
                cutoff=cutoff, rm.unmapped=rm.unmapped,
                feature.type=feature.type))
    } else if (feature.type == 'filtered'){
        # if not yet there, stop
        if (is.null(filt_feat(siamcat, verbose=0))){
            stop("Features have not yet been filtered, exiting...\n")
        }
        feat <- get.filt_feat.matrix(siamcat)
        param.set <- filt_params(siamcat)
        param.set[[length(param.set)+1]] <-
            list(filter.method=filter.method,
                cutoff=cutoff, rm.unmapped=rm.unmapped,
                feature.type=feature.type)
    } else if (feature.type == 'normalized'){
        # if not yet there, stop
        if (is.null(norm_feat(siamcat, verbose=0))){
            stop("Features have not yet been normalized, exiting...\n")
        }
        if (is.null(filt_feat(siamcat, verbose=0))){
            param.set <- list(list(filter.method=filter.method,
                    cutoff=cutoff, rm.unmapped=rm.unmapped,
                    feature.type=feature.type, feature.type=feature.type))
        } else {
            param.set <- filt_params(siamcat)
            param.set[[length(param.set)+1]] <-
                list(filter.method=filter.method,
                    cutoff=cutoff, rm.unmapped=rm.unmapped,
                    feature.type=feature.type)
        }
        feat <- get.norm_feat.matrix(siamcat)
    }

    # check if there are NAs in the data
    if (any(is.na(feat))){
        stop("There are NAs in the feature matrix! Exiting...")
    }
    if (verbose > 1)
        message(paste("+++ before filtering, the data have",
            nrow(feat), "features"))

    ### apply filters
    if (verbose > 2)
        message(paste("+++ applying", filter.method, "filter"))
    if (filter.method == "abundance") {
        # remove features whose abundance is never above the threshold value
        # (e.g. 0.5%) in any of the samples
        f.max <- rowMaxs(feat)
        f.idx <- which(f.max >= cutoff)
    } else if (filter.method == "cum.abundance") {
        # remove features with very low abundance in all samples i.e. ones that
        # are never among the most abundant entities that collectively make up
        # (1-cutoff) of the reads in any sample
        f.idx <- vector("numeric", 0)
        # sort features per sample and apply cumsum to identify how many
        # collectively have weight K
        for (s in seq_len(ncol(feat))) {
            srt <- sort(feat[, s], index.return = TRUE)
            cs <- cumsum(srt$x)
            m <- max(which(cs < cutoff))
            f.idx <- union(f.idx, srt$ix[-(seq_len(m))])
        }
        # an index of those features that collectively make up more than 1-K of
        # the read mass in any sample
        f.idx <- sort(f.idx)
    } else if (filter.method == "prevalence") {
        # remove features with low prevalence across samples i.e. ones that are
        # 0 (undetected) in more than (1-cutoff)
        # proportion of samples
        f.idx <-
            which(rowSums(feat > 0) / ncol(feat) > cutoff)
    } else if (filter.method == "variance"){
        # remove features with very low variance
        f.var <- rowVars(feat)
        f.idx <- which(f.var >= cutoff)
    }


    if (verbose > 2)
        message("+++ checking for unmapped reads")
    ### postprocessing and output generation
    if (rm.unmapped) {
        # remove 'unmapped' feature
        names.unmapped <- c("UNMAPPED", "-1", "X.1", "unmapped",
            "UNCLASSIFIED", "unclassified", "UNASSIGNED", "unassigned")

        unm.idx <- rownames(feat) %in% names.unmapped

        if (any(unm.idx)) {
            f.idx <- f.idx[-which(f.idx %in% which(unm.idx))]
            if (verbose > 2)
                message(paste("+++ removing",
                    rownames(feat)[unm.idx], "as unmapped reads"))
            if (verbose > 1)
                message(paste("+++ removed", sum(unm.idx),
                    "features corresponding to UNMAPPED reads"))
        } else {
            if (verbose > 1)
                message(paste0("+++ tried to remove unmapped reads ",
                    "but could not find any. Continue anyway."))
        }
    }

    if (verbose > 1)
        message(paste0("+++ removed ",
            nrow(feat) - length(f.idx) - sum(unm.idx),
            " features whose values did not exceed ", cutoff,
            " in any sample (retaining ", length(f.idx), ")" ))

    f.names <- rownames(feat)[f.idx]
    if (length(f.idx) == 0){
        stop('No features retained after filtering!',
             ' Try changing your cutoff. Exiting...\n')
    }


    if (verbose > 2)
        message("+++ saving filtered features")

    filt_feat(siamcat) <- new("filt_feat", filt.feat=otu_table(feat[f.names,],
            taxa_are_rows=TRUE), filt.param=param.set)

    e.time <- proc.time()[3]
    if (verbose > 1)
        message(paste("+ finished filter.features in",
            formatC(e.time - s.time, digits = 3), "s"))

    if (verbose == 1)
        message("Features successfully filtered")

    return(siamcat)
}
