#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Validate samples in labels, features, and metadata
#' @description This function checks if labels are available for all samples in
#'        features. Additionally validates metadata, if available.
#' @param siamcat an object of class \link{siamcat-class}
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'        for only information about progress and success, \code{2} for normal
#'        level of information and \code{3} for full debug information,
#'        defaults to \code{1}
#' @keywords SIAMCAT validate.data
#' @export
#' @details This function validates the data by checking that labels are
#'        available for all samples in the feature matrix. Furthermore,
#'        the number of samples per class is checked to ensure a minimum
#'        number. If metadata is available, the overlap between labels and
#'        metadata is checked as well.
#' @return an object of class \link{siamcat-class} with validated data
#' @examples
#'
#'  data(siamcat_example)
#'  # simple working example
#'  siamcat_validated <- validate.data(siamcat_example)
#'
validate.data <- function(siamcat, verbose = 1) {
    if (verbose > 1)
        message("+ starting validate.data")
    s.time <- proc.time()[3]
    # Check if labels are available for all samples in features
    if (verbose > 2) {
        message("+++ checking if labels are available for all samples in features")
    }
    if (length(label(siamcat)@label) == ncol(features(siamcat))) {
        stopifnot(all(names(label(siamcat)) %in% colnames(features(siamcat))) && all(colnames(features(siamcat)) %in%
            names(label(siamcat))))
        # if of the same length, everything should match and be in the same order
        m <- match(names(label(siamcat)), colnames(features(siamcat)))
        features(siamcat) <- features(siamcat)[, m]
        stopifnot(all(names(label(siamcat)) == colnames(features(siamcat))))

    } else if (length(label(siamcat)) >= ncol(features(siamcat))) {
        # if there are more labels than samples in features, remove them in labels
        stopifnot(all(colnames(features(siamcat)) %in% names(label(siamcat))))

        ids <- colnames(features(siamcat))[colnames(features(siamcat)) %in% names(label(siamcat))]
        # create new labels object with reduced labels
        siamcat <- filter.label(siamcat, ids)

    } else if (length(label(siamcat)) <= ncol(features(siamcat))) {
        # if there are more samples in features, remove them and keep only the ones for which labels are available
        stopifnot(all(names(label(siamcat)) %in% colnames(features(siamcat))))

        message(paste("Warning: Removing", ncol(features(siamcat)) - length(label(siamcat)), "sample(s) for which no labels are available."))

        m <- match(names(label(siamcat)), colnames(features(siamcat)))
        features(siamcat) <- features(siamcat)[, m]

    }

    # Check for sample number in the different classes
    if (verbose > 2)
        message("+++ checking sample number per class")
    for (i in label(siamcat)@info$class.descr) {
        if (sum(label(siamcat) == i) <= 5) {
            stop("Data set has only", sum(label(siamcat) == i), "training examples of class", i, " This is not enough for SIAMCAT to proceed")
        }
        if (sum(label(siamcat) == i) < 10) {
            message(paste("Data set has only", sum(label(siamcat) == i), "training examples of class", i, " . Note that a dataset this small/skewed is not necessarily
        suitable for analysis in this pipeline."))
        }
    }


    # if meta(siamcat)data is available, check for overlap in labels
    if (!is.null(meta(siamcat))) {
        if (verbose > 2)
            message("+++ check for overlap between labels and metadata")
        if (length(label(siamcat)) == nrow(meta(siamcat))) {
            stopifnot(all(names(label(siamcat)) %in% rownames(meta(siamcat))) && all(rownames(meta(siamcat)) %in%
                names(label(siamcat))))
            m <- match(names(label(siamcat)), rownames(meta(siamcat)))
            meta(siamcat) <- meta(siamcat)[m, ]
            stopifnot(all(names(label(siamcat)) == rownames(meta(siamcat))))
        } else if (length(label(siamcat)) <= nrow(meta(siamcat))) {
            stopifnot(all(names(label(siamcat)) %in% rownames(meta(siamcat))))
            m <- match(names(label(siamcat)), rownames(meta(siamcat)))
            message(paste("Warning: Removing metadata information for", nrow(meta(siamcat)) - length(label(siamcat)),
                "superfluous samples."))
            meta(siamcat) <- meta(siamcat)[m, ]
        } else {
            stop("! Metadata is not available for all samples!")
        }
    }
    orig_feat(siamcat) <- otu_table(siamcat@phyloseq)
    e.time <- proc.time()[3]
    if (verbose > 1)
        message(paste("+ finished validate.data in", formatC(e.time - s.time, digits=3), "s"))
    if (verbose == 1)
        message("Data succesfully validated")
    return(siamcat)
}
