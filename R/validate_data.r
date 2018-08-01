#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Validate samples in labels, features, and metadata
#' @description This function checks if labels are available for all samples in
#'     features. Additionally validates metadata, if available.
#' @param siamcat an object of class \link{siamcat-class}
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'     for only information about progress and success, \code{2} for normal
#'     level of information and \code{3} for full debug information,
#'     defaults to \code{1}
#' @keywords SIAMCAT validate.data
#' @export
#' @details This function validates the data by checking that labels are
#'     available for all samples in the feature matrix. Furthermore,
#'     the number of samples per class is checked to ensure a minimum
#'     number. If metadata is available, the overlap between labels and
#'     metadata is checked as well.
#' @return an object of class \link{siamcat-class} with validated data
#' @examples
#'
#'     data(siamcat_example)
#'     # simple working example
#'     siamcat_validated <- validate.data(siamcat_example)
#'
validate.data <- function(siamcat, verbose = 1) {
    if (verbose > 1)
        message("+ starting validate.data")
    label <- label(siamcat)
    feat  <- features(siamcat)
    meta  <- meta(siamcat)
    s.time <- proc.time()[3]

    # Check if features contain any missing values, i.e. no NAs
    if (any(is.na(feat))){
        stop('### The features contain missing data! Exiting...')
    }
    # check for compositional data
    if (any(colSums(feat) > 1.01) || any(colSums(feat) < 0.99)) {
        message('\t### Warning: The data does not seem to be compositional!')
    }

    # Check if labels are available for all samples in features
    if (verbose > 2) {
        message("+++ checking overlap between labels and features")
    }
    s.intersect <- intersect(names(label$label), colnames(feat))
    # check and re-order features
    s.removed <- ncol(feat) - length(s.intersect)
    features(siamcat) <- feat[,s.intersect]
    feat <- features(siamcat)
    if (verbose > 1 & s.removed != 0)
        message(paste0("+ Removed ", s.removed,
                       " samples from the feature matrix..."))
    # check and re-order labels
    s.removed <- length(label$label) - length(s.intersect)
    ids <- match(s.intersect, names(label$label))
    siamcat <- filter.label(siamcat, ids=ids, verbose=verbose)
    label <- label(siamcat)
    stopifnot(all(names(label$label) == colnames(feat)))
    if (verbose > 1 & s.removed != 0)
        message(paste0("+ Removed ", s.removed,
                       " samples from the label object..."))

    # Check for sample number in the different classes
    if (verbose > 2)
        message("+++ checking sample number per class")
    for (i in seq_along(label$info)) {
        if (sum(label$label == label$info[i]) <= 5) {
            stop(
                "Data set has only",
                sum(label$label == label$info[i]),
                "training examples
                of class",
                names(label$info)[i],
                " This is not enough for SIAMCAT to proceed"
            )
        }
        if (sum(label$label == label$info[i]) < 10) {
            message(
                paste(
                    "Data set has only",
                    sum(label$label == label$info[i]),
                    "training
                    examples of class",
                    names(label$info)[i],
                    " . Note that a dataset this small/skewed
                    is not necessarily
                    suitable for analysis in this pipeline."
                )
            )
        }
    }


    # if metadata is available, check for overlap in labels
    if (!is.null(meta)) {
        if (verbose > 2)
            message("+++ checking overlap between samples and metadata")
        if (!all(names(label$label) %in% rownames(meta))){
            stop('Metadata is not available for all samples! Exiting...')
        }

        s.intersect <- intersect(names(label$label), rownames(meta))
        # check and re-order metadata
        s.removed <- nrow(meta) - length(s.intersect)
        meta(siamcat) <- meta[s.intersect,]
        if (verbose > 1 & s.removed != 0)
            message(paste0("+ Removed ", s.removed,
                           " samples from the metadata..."))
        stopifnot(all(names(label$label) == rownames(meta(siamcat))))
    }

    orig_feat(siamcat) <- otu_table(physeq(siamcat))
    e.time <- proc.time()[3]
    if (verbose > 1)
        message(paste(
            "+ finished validate.data in",
            formatC(e.time - s.time, digits = 3),
            "s"
        ))
    if (verbose == 1)
        message("Data succesfully validated")
    return(siamcat)
}
