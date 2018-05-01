#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Perform unsupervised feature filtering.
#'
#' @description This function performs unsupervised feature filtering. Features
#'     can be filtered based on abundance or prevalence. Additionally,
#'     unmapped reads may be removed.
#'
#' @usage filter.features(siamcat, filter.method = "abundance",
#'     cutoff = 0.001, recomp.prop = FALSE, rm.unmapped = TRUE, verbose = 1)
#'
#' @param siamcat an object of class \link{siamcat-class}
#'
#' @param filter.method method used for filtering the features, can be one of
#'     these: \code{c('abundance', 'cum.abundance', 'prevalence')},
#'     defaults to \code{'abundance'}
#'
#' @param cutoff float, abundace or prevalence cutoff, default to \code{0.001}
#'
#' @param recomp.prop boolean, should relative abundances be recomputed?,
#'     defaults to \code{FALSE}
#'
#' @param rm.unmapped boolean, should unmapped reads be discarded?, defaults to
#'     \code{TRUE}
#'
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'     for only information about progress and success, \code{2} for normal
#'     level of information and \code{3} for full debug information,
#'     defaults to \code{1}
#'
#' @keywords SIAMCAT filter.features
#'
#' @details This function filters the features in a \link{siamcat-class}
#'     object in a unsupervised manner.
#'
#'     The different filter methods work in the following way: \itemize{
#'     \item \code{'abundace'} remove features whose abundance is never
#'     above the threshold value in any of the samples
#'     \item \code{'cum.abundance'} remove features with very low abundance
#'     in all samples i.e. ones that are never among the most abundant
#'     entities that collectively make up (1-cutoff) of the reads in
#'     any sample
#'     \item \code{'prevalence'} remove features with low prevalence across
#'     samples i.e. ones that are 0 (undetected) in more than (1-cutoff)
#'     proportion of samples.
#'     }
#' @export
#'
#' @return siamcat an object of class \link{siamcat-class}
#'
#' @examples
#'     # Example dataset
#'     data(siamcat_example)
#'     # since the whole pipeline has been run in the example data, the feature
#'     # were filtered already.
#'     siamcat_example <- reset.features(siamcat_example)
#'
#' # Simple examples
#' siamcat_filtered <- filter.features(siamcat_example,
#'     filter.method='abundance',
#'     cutoff=1e-03)
#'
filter.features <- function(siamcat,
    filter.method = "abundance",
    cutoff = 0.001,
    recomp.prop = FALSE,
    rm.unmapped = TRUE,
    verbose = 1) {
    ### this statement does not have the purpose to calculate relative
    ### abundances on the fly and return them.  Instead, it's purpose is to be
    ### able to calculate f.idx (specifying the indices of features which are
    ### to be kept) when feature list has already been transformed to relative
    ### abundances, but e.g. certain features have been removed manually.

    if (verbose > 1)
        message("+ starting filter.features")
    s.time <- proc.time()[3]

    if (!filter.method %in% c("abundance", "cum.abundace", "prevalence")) {
        stop("! Unrecognized filter.method, exiting!\n")
    }

    if (verbose > 1)
        message(paste(
            "+++ before filtering, the data has",
            nrow(features(siamcat)),
            "features"
        ))
    if (recomp.prop) {
        # recompute relative abundance values (proportions)
        ra.feat           <- prop.table(features(siamcat), 2)
        features(siamcat) <-
            otu_table(ra.feat, taxa_are_rows = TRUE)
    } else {
        ra.feat <- get.features.matrix(siamcat)
    }

    ### apply filters
    if (verbose > 2)
        message(paste("+++ applying", filter.method, "filter"))
    if (filter.method == "abundance") {
        # remove features whose abundance is never above the threshold value
        # (e.g. 0.5%) in any of the samples
        f.max <- rowMaxs(ra.feat)
        f.idx <- which(f.max >= cutoff)
    } else if (filter.method == "cum.abundance") {
        # remove features with very low abundance in all samples i.e. ones that
        # are never among the most abundant entities that collectively make up
        # (1-cutoff) of the reads in any sample
        f.idx <- vector("numeric", 0)
        # sort features per sample and apply cumsum to identify how many
        # collectively have weight K
        for (s in seq_len(ncol(ra.feat))) {
            srt <- sort(ra.feat[, s], index.return = TRUE)
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
            which(rowSums(ra.feat > 0) / ncol(ra.feat) > cutoff)
    }


    if (verbose > 2)
        message("+++ checking for unmapped reads")
    ### postprocessing and output generation
    if (rm.unmapped) {
        # remove 'unmapped' feature
        names.unmapped <- c(
            "UNMAPPED",
            "-1",
            "X.1",
            "unmapped",
            "UNCLASSIFIED",
            "unclassified",
            "UNASSIGNED",
            "unassigned"
        )
        unm.idx <- rownames(features(siamcat)) %in% names.unmapped
        if (any(unm.idx)) {
            f.idx <- f.idx[-which(f.idx %in% which(unm.idx))]
            if (verbose > 2)
                message(paste(
                    "+++ removing",
                    rownames(features(siamcat))[unm.idx],
                    "as unmapped reads"
                ))
            if (verbose > 1)
                message(paste(
                    "+++ removed",
                    sum(unm.idx),
                    "features corresponding to UNMAPPED reads"
                ))
        } else {
            if (verbose > 1)
                message(
                    "+++ tried to remove unmapped reads, but could not find
                    them. Continue anyway."
                )
        }
        }
    if (verbose > 2)
        message("+++ applying prune_taxa")
    if (verbose > 1)
        message(
            paste0(
                "+++ removed ",
                nrow(features(siamcat)) - length(f.idx) - sum(unm.idx),
                " features whose values did not exceed ",
                cutoff,
                " in any sample (retaining ",
                length(f.idx),
                ")"
            )
        )
    f.names <- rownames(features(siamcat))[f.idx]
    physeq(siamcat) <-
        prune_taxa(x = physeq(siamcat), taxa = f.names)
    e.time <- proc.time()[3]
    if (verbose > 1)
        message(paste(
            "+ finished filter.features in",
            formatC(e.time - s.time, digits = 3),
            "s"
        ))
    if (verbose == 1)
        message("Features successfully filtered")
    return(siamcat)
}
