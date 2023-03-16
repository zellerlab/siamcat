#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Add metadata as predictors
#'
#' @description This function adds metadata to the feature matrix to be later
#' used as predictors
#'
#' @usage add.meta.pred(siamcat, pred.names, std.meta = TRUE, 
#' feature.type='normalized', verbose = 1)
#'
#' @param siamcat object of class \link{siamcat-class}
#'
#' @param pred.names vector of names of the variables within the metadata to 
#' be added to the feature matrix as predictors
#'
#' @param std.meta boolean, should added metadata features be standardized?,
#' defaults to \code{TRUE}
#'
#' @param feature.type string, on which type of features should the function 
#' work? Can be either \code{"original"}, \code{"filtered"}, or 
#' \code{"normalized"}. Please only change this paramter if you know what
#' you are doing!
#'
#' @param verbose integer, control output: \code{0} for no output at all, 
#' \code{1} for only information about progress and success, \code{2} for
#' normal level of information and \code{3} for full debug information,
#' defaults to \code{1}
#'
#' @keywords SIAMCAT add.meta.pred
#'
#' @export
#' 
#' @encoding UTF-8
#'
#' @details This functions adds one or several metadata variables to the set
#' of features, so that they can be included for model training.
#'
#' Usually, this function should be called before \link{train.model}.
#'
#' Numerical meta-variables are added as z-scores to the feature matrix unless
#' specified otherwise.
#'
#' Please be aware, that non-numerical metadata variables will be converted to
#' numerical values by using \code{as.numeric()} and could therefore lead to
#' errors. Thus, it makes sense to encode non-numerical metadata variables to
#' numerically before you start the SIAMCAT workflow.
#'
#' @return an object of class \link{siamcat-class} with metadata added to the
#' features
#'
#' @examples
#' data(siamcat_example)
#'
#' # Add the Age of the patients as potential predictor
#' siamcat_age_added <- add.meta.pred(siamcat_example, pred.names=c('Age'))
#'
#' # Add Age and BMI as potential predictors
#' # Additionally, prevent standardization of the added features
#' siamcat_meta_added <- add.meta.pred(siamcat_example, 
#'     pred.names=c('Age', 'BMI'), std.meta=FALSE)
add.meta.pred <- function(siamcat, pred.names, std.meta = TRUE,
    feature.type = 'normalized', verbose = 1) {

    if (verbose > 1)
    message("+ starting add.meta.pred")
    s.time <- proc.time()[3]

    if (is.null(meta(siamcat))) {
        stop('SIAMCAT object has no metadata. Exiting...')
    }
    if (!feature.type %in% c('original', 'filtered', 'normalized')){
        stop("Unrecognised feature type, exiting...\n")
    }
    if (feature.type != 'normalized'){
        warning('It is recommended to add a meta-predictor only',
            ' after feature normalization.')
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

    ### add metadata as predictors to the feature matrix
    cnt <- 0

    if (all(pred.names != "") && all(!is.null(pred.names))) {
        if (verbose > 2)
        message("+ starting to add metadata predictors")
        for (p in pred.names) {
            if (verbose > 2)
            message("+++ adding metadata predictor: ", p)
            if (!p %in% colnames(meta(siamcat))) {
                stop("There is no metadata variable called ", p)
            }
            if (paste0('META_', toupper(p)) %in% rownames(feat)){
                stop("This meta-variable has already been added. Exiting...")
            }
            idx <- which(colnames(meta(siamcat)) == p)
            if (length(idx) != 1)
            stop(p, "matches multiple columns in the metada")

            m <- unlist(meta(siamcat)[, idx])

            # check if the meta-variable is a factor or numeric
            if (is.factor(m)){
                warning(paste0('WARNING: meta-variable ', p,' is a factor and',
                ' not numeric...\n   The values will be converted to numerical',
                ' values with "as.numeric()"\n'))
                m <- as.numeric(m)
            }
            stopifnot(is.numeric(m))

            if (!all(is.finite(m))) {
                na.cnt <- sum(!is.finite(m))
                if (verbose > 1)
                    message(paste("++++ filling in", na.cnt,
                        "missing values by mean imputation"))
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
            feat <- rbind(feat, m)
            rownames(feat)[nrow(feat)] <- paste0('META_', toupper(p))
            cnt <- cnt + 1
        }

        if (feature.type == 'original'){
            orig_feat(siamcat) <- otu_table(feat, taxa_are_rows=TRUE)
        } else if (feature.type == 'filtered'){
            filt_feat(siamcat) <- list(
                filt.feat=feat,
                filt.param=filt_params(siamcat))
        } else if (feature.type == 'normalized'){
            norm_feat(siamcat) <- list(
                norm.feat=feat,
                norm.param=norm_params(siamcat))
        }

        if (verbose > 1)
        message(paste("+++ added", cnt,
            "meta-variables as predictor to the feature matrix"))
        } else {
            if (verbose > 0)
            message("+++ Not adding any of the meta-variables as predictor",
                " to the feature matrix")
        }
        e.time <- proc.time()[3]
        if (verbose > 1)
        message(paste("+ finished add.meta.pred in", formatC(e.time - s.time,
            digits = 3), "s"))
        if (verbose == 1)
        message("Adding metadata as predictor finished")
        return(siamcat)
}
