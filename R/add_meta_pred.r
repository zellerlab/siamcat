#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Add metadata as predictors
#' @description This function adds metadata to the feature matrix to be
#'        later used as predictors
#' @param siamcat object of class \link{siamcat-class}
#' @param pred.names vector of names of the variables within the metadata to be
#'        added to the feature matrix as predictors
#' @param std.meta boolean, should added metadata features be standardized?,
#'        defaults to \code{TRUE}
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'        for only information about progress and success, \code{2} for normal
#'        level of information and \code{3} for full debug information,
#'        defaults to \code{1}
#' @keywords SIAMCAT add.meta.pred
#' @export
#' @return an object of class \link{siamcat-class} with metadata added to the
#'        features
#' @examples
#'  data(siamcat_example)
#'  # Add the Age of the patients as potential predictor
#'  siamcat_age_added <- add.meta.pred(siamcat_example, pred.names=c('age'))
#'
#'  # Add Age, BMI, and Gender as potential predictors
#'  # Additionally, prevent standardization of the added features
#'  siamcat_meta_added <- add.meta.pred(siamcat_example, pred.names=c('age',
#'      'bmi', 'gender'), std.meta=FALSE)
add.meta.pred <- function(siamcat, pred.names = NULL, std.meta = TRUE, verbose = 1) {
    if (verbose > 1)
        message("+ starting add.meta.pred")
    s.time <- proc.time()[3]
    ### add metadata as predictors to the feature matrix
    cnt <- 0

    if (pred.names != "" && !is.null(pred.names)) {
        if (verbose > 2)
            message("+ starting to add metadata predictors")
        for (p in pred.names) {
            if (verbose > 2)
                message("+++ adding metadata predictor:", p)
            if (!p %in% colnames(meta(siamcat))) {
                stop("There is no metadata variable called ", p)
            }
            idx <- which(colnames(meta(siamcat)) == p)
            if (length(idx) != 1)
                stop(p, "matches multiple columns in the metada")

            m <- unlist(meta(siamcat)[, idx])

            if (!all(is.finite(m))) {
                na.cnt <- sum(!is.finite(m))
                if (verbose > 1)
                  message(paste("++++ filling in", na.cnt, "missing values by mean imputation"))
                mn <- mean(m, na.rm = TRUE)
                m[!is.finite(m)] <- mn
            }

            if (std.meta) {
                if (verbose > 1)
                  message(paste("++++ standardizing metadata feature", p))
                m.mean <- mean(m, na.rm = TRUE)
                m.sd <- sd(m, na.rm = TRUE)
                stopifnot(!m.sd == 0)
                m <- (m - m.mean)/m.sd
            }
            features.with.meta <- otu_table(rbind(features(siamcat),m),
                taxa_are_rows = TRUE)
            rownames(features.with.meta)[nrow(features.with.meta)] <- paste(
            "META_", toupper(p), sep = "")
            features(siamcat) <- features.with.meta

            cnt <- cnt + 1
        }
        if (verbose > 1)
            message(paste("+++ added", cnt, "meta-variables as predictor to the feature matrix"))
    } else {
        if (verbose > 0)
            message("+++ Not adding any of the meta-variables as predictor to the feature matrix")
    }
    e.time <- proc.time()[3]
    if (verbose > 1)
        message(paste("+ finished add.meta.pred in", formatC(e.time - s.time, digits = 3), "s"))
    if (verbose == 1)
        message("Adding metadata as predictor finished")
    return(siamcat)
}
