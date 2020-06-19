#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

# ##############################################################################
# label object

# check label object for validity
#'@keywords internal
check.label <- function(object){
    errors <- character()
    # check that all entries are there
    if (!all(names(object) == c('label', 'info', 'type'))){
        msg <- 'Label object does not contain all needed entries!'
        errors <- c(errors, msg)
    }
    # check that label is binary or test
    if (object$type != 'BINARY' & object$type != 'TEST'){
        msg <- 'Label object is neither binary nor a test label!'
        errors <- c(errors, msg)
    }
    # check that info and label match up
    if (!all(sort(unique(object$label)) == object$info)){
        msg <- 'label info does not match to label entries!'
        errors <- c(errors, msg)
    }
    # check that info has names
    if (is.null(names(object$info))){
        msg <- 'Label info does not contain group names!'
        errors <- c(errors, msg)
    }
    if (length(errors) == 0) NULL else errors
}

# ##############################################################################
# Filtered data/parameters

# check filtered data/parameters for validity
#'@keywords internal
check.filt.feat <- function(object){
    errors <- character()
    if (is.null(object$filt.param)){
        msg <- "Filtering parameters are missing!"
        errors <- c(errors, msg)
    }
    # check if all entries are lists
    if (!all(vapply(object$filt.param, class,
        FUN.VALUE = character(1)) == 'list')){
            msg <- "Filtering parameters are not in the right format!"
            errors <- c(errors, msg)
        }
    # check if each entry is of length 3
    if (!all(vapply(object$filt.param, length, FUN.VALUE = numeric(1)) == 4)){
        msg <- "Filtering parameters are not in the right format!"
        errors <- c(errors, msg)
    }
    # check if they have the correct entries
    if (!all(vapply(object$filt.param, names,
        FUN.VALUE = character(4)) ==
        c('filter.method', 'cutoff', 'rm.unmapped', 'feature.type'))){
            msg <- "Filtering parameters do not contain all needed entries!"
            errors <- c(errors, msg)
        }

    # check that taxa are rows == TRUE
    if (is.null(object$filt.feat)){
        msg <- "Filtered features are missing!"
        errors <- c(errors, msg)
    }
    if (length(errors) == 0) NULL else errors
}

# ##############################################################################
# Associations

# check associations for validity
#'@keywords internal
check.assoc <- function(object){
    errors <- character()
    if (is.null(object$assoc.param)){
        msg <- "Association parameters are missing!"
        errors <- c(errors, msg)
    }
    if (is.null(object$assoc.results)){
        msg <- "Association results are missing!"
        errors <- c(errors, msg)
    }
    # check that assoc.param contains all entries
    if (!all(names(object$assoc.param) == c('detect.lim', 'pr.cutoff',
        'probs.fc', 'mult.corr', 'alpha', 'feature.type'))){
            msg <- 'Association testing parameters do not contain all entries!'
            errors <- c(errors, msg)
        }
    # check that all entries are valid and in the expected ranges
    if (!all(vapply(object$assoc.param, class,
        FUN.VALUE=character(1)) == c('numeric', 'numeric', 'numeric',
            'character', 'numeric', 'character'))){
    msg<-'Association testing parameters do not contain the expected classes!'
    errors <- c(errors, msg)
    }
    # detect.lim
    if (object$assoc.param$detect.lim > 1 | object$assoc.param$detect.lim < 0){
        msg<-'Detection limit (pseudocount) is not valid (not between 1 and 0)!'
        errors <- c(errors, msg)
    }
    # pr.cutoff
    if (object$assoc.param$pr.cutoff > 1 | object$assoc.param$pr.cutoff < 0){
        msg<-'Prevalence cutoff is not valid (not between 1 and 0)!'
        errors <- c(errors, msg)
    }
    # probs.fc
    if (any(object$assoc.param$probs.fc > 1) |
        any(object$assoc.param$probs.fc < 0)){
        msg<-'Quantiles for FC calculation are not valid (not between 1 and 0)!'
        errors <- c(errors, msg)
    }
    # mult.corr
    if (!object$assoc.param$mult.corr %in%
        c('none', 'bonferroni', 'holm', 'fdr', 'bhy')){
            msg<-'Multiple testing correction method not valid!'
            errors <- c(errors, msg)
    }
    # alpha
    if (object$assoc.param$alpha > 1 | object$assoc.param$alpha < 0){
        msg<-'Significance level (alpha) is not valid (not between 1 and 0)!'
        errors <- c(errors, msg)
    }
    # feat.type
    if (!object$assoc.param$feature.type %in%
        c('filtered', 'original', 'normalized')){
        msg <- 'Feature type is not valid!
            (should be original, filtered, or normalized)'
        errors <- c(errors, msg)
    }
    # check that assoc.results contains all that it should
    if (!all(colnames(object$assoc.results)==c("fc", "p.val", "auc",
        "auc.ci.l", "auc.ci.h", "pr.shift", "pr.n", "pr.p", "bcol", "p.adj"))){
        msg <- 'Association results do not contain all needed entries!'
        errors <- c(errors, msg)
    }
    if (nrow(object$assoc.results) < 1){
        msg <- 'Association results are empty!'
        errors <- c(errors, msg)
    }
    if (length(errors) == 0) NULL else errors
}

# ##############################################################################
# Normalization data/parameters

# check normalization data/parameters for validity
#'@keywords internal
check.norm.feat <- function(object){
    errors <- character()
    if (is.null(object$norm.param)){
        msg <- "Normalization parameters are missing!"
        errors <- c(errors, msg)
    }
    norm.method <- object$norm.param$norm.method
    # check if norm method is there
    if (is.null(norm.method)){
        msg <- 'norm.method is NULL!'
        errors <- c(errors, msg)
    }
    if (!norm.method %in%
        c('rank.std','rank.unit','log.unit',
            'log.std','log.clr', 'std', 'pass')){
            mgs <- paste0('norm.method ', norm.method, 'not recognized!')
            errors <- c(errors, msg)
        }
    # check if retained feat is there
    if (is.null(object$norm.param$retained.feat)){
        msg <- 'retained.feat is missing!'
        errors <- c(errors, msg)
    }
    # for each norm method, check additional entries
    if (startsWith(norm.method, 'log')){
        if (is.null(object$norm.param$log.n0)){
            msg <- 'Detection limit (pseudocount) is missing!'
            errors <- c(errors, msg)
        }
        if (object$norm.param$log.n0 > 1 | object$norm.param$log.n0 < 0){
            msg<-'Detection limit (pseudocount) is not between 0 and 1!'
            errors <- c(errors, msg)
        }
    }
    # std
    if (endsWith(norm.method, 'std')){
        # feat.mean or feat.adj.sd
        if (is.null(object$norm.param$feat.mean) |
            is.null(object$norm.param$feat.adj.sd)){
            msg<-'Needed entries for std-normalization are missing!'
            errors <- c(errors, msg)
        }
    }
    # log.unit n.p, norm.margin, norm.fun, feat.norm.denom
    if (norm.method == 'log.unit'){
        if (is.null(object$norm.param$n.p)){
            msg<-'Vector norm is missing!'
            errors <- c(errors, msg)
        }
        if (is.null(object$norm.param$norm.fun)){
            msg<-'Normalization function is missing!'
            errors <- c(errors, msg)
        }
        if (is.null(object$norm.param$norm.margin)){
            msg<-'Normalization margin is missing!'
            errors <- c(errors, msg)
        }
        if (is.null(object$norm.param$feat.norm.denom)){
            msg<-'Feature specific denominator is missing!'
            errors <- c(errors, msg)
        }
    }
    if (is.null(object$norm.feat)){
        msg <- "Normalized features are missing!"
        errors <- c(errors, msg)
    }
    if (length(errors) == 0) NULL else errors
}

# ##############################################################################
# Data split

# check data split for validity
#'@keywords internal
check.data.split <- function(object){
    errors <- character()
    # check names
    if (!all(names(object) == c('training.folds', 'test.folds',
                                'num.resample', 'num.folds'))){
        msg <- 'Data split does not contain all needed entries!'
        errors <- c(errors, msg)
    }
    # check that num.resample and num.folds are numbers
    if (length(object$num.resample) != 1 |
        !is.numeric(object$num.resample)){
        msg <- 'num.resample should be numeric and of length 1!'
        errors <- c(errors, msg)
    }
    if (length(object$num.folds) != 1 |
        !is.numeric(object$num.folds)){
        msg <- 'num.folds should be numeric and of length 1!'
        errors <- c(errors, msg)
    }
    # check that training.folds is a list (of the right length)
    if (!is.list(object$training.folds)){
        msg <- 'training.folds should be a list!'
        errors <- c(errors, msg)
    }
    if (length(object$training.folds) != object$num.resample){
        msg <- 'training.folds should be of length num.resample!'
        errors <- c(errors, msg)
    }
    if (!all(vapply(object$training.folds, length,
            FUN.VALUE = numeric(1)) == object$num.folds)){
        msg <- 'All training.folds should be of length num.folds!'
        errors <- c(errors, msg)
    }
    # same for test.folds
    if (!is.list(object$test.folds)){
        msg <- 'test.folds should be a list!'
        errors <- c(errors, msg)
    }
    if (length(object$test.folds) != object$num.resample){
        msg <- 'test.folds should be of length num.resample!'
        errors <- c(errors, msg)
    }
    if (!all(vapply(object$test.folds, length,
            FUN.VALUE = numeric(1)) == object$num.folds)){
        msg <- 'All test.folds should be of length num.folds!'
        errors <- c(errors, msg)
    }
    # check that no samples are in training and test fold at the same time
    test.consistency <- any(c(vapply(seq_len(object$num.resample),
        FUN=function(x){vapply(seq_len(object$num.folds),
            FUN=function(y){
                any(object$test.folds[[x]][[y]] %in%
                    object$training.folds[[x]][[y]])
                }, FUN.VALUE = logical(1))
            }, FUN.VALUE=logical(object$num.folds))))
    if (test.consistency){
        msg <- 'Some samples are in both training and test folds!'
        errors <- c(errors, msg)
    }


    if (length(errors) == 0) NULL else errors
}

# ##############################################################################
# Model list

# check model list for validity
#'@keywords internal
check.model.list <- function(object){
    errors <- character()
    if (!all(c('model.type', 'feature.type', 'models') %in% names(object))){
        msg <- 'Model list does not contain all needed entries!'
        errors <- c(errors, msg)
    }
    # check model.type
    length.model.type <- length(object$model.type)
    if (length.model.type != 1){
        msg <- paste0('Model type is of length ', length.model.type, '.',
            ' Should be length 1!')
        errors <- c(errors, msg)
    }
    # check that all models in the list are mlr models
    if (any(vapply(object$models, class,
            FUN.VALUE = character(1)) != 'WrappedModel')){
        msg <- 'Models are supposed to be mlr-WrappedModels!'
        errors <- c(errors, msg)
    }
    # check feature type
    if (!object$feature.type %in% c('original', 'filtered', 'normalized')){
        msg <- paste0('Feature type ', object$feature.type,
            ' is not recognized!
            Should be either original, filtered, or normalized!')
        errors <- c(errors, msg)
    }
    if (length(errors) == 0) NULL else errors
}

# ##############################################################################
# Prediction matrix

# check.prediction.matrix function

# ##############################################################################
# Evaluation data

# check evaluation data for validity
#'@keywords internal
check.eval.data <- function(object){
    errors <- character()
    # check that all entries are there
    if (!all(c('roc', 'auroc', 'prc', 'auprc', 'ev') %in% names(object))){
        msg <- 'Not all needed entries are given!'
        errors <- c(errors, msg)
    }
    # check roc
    if (!is(object$roc,'roc')){
        msg <- 'Entry for roc is not an object of class roc (from pROC)!'
        errors <- c(errors, msg)
    }
    # check prc
    if (!is.list(object$prc)  |
        !all(names(object$prc) ==c('recall', 'precision')) |
        length(unique(vapply(object$prc, length, FUN.VALUE=numeric(1))))!=1){
        msg <- paste0('No valid entry for prc ',
            '(missing entries or no list with entries of equal length)!')
        errors <- c(errors, msg)
    }
    # check ev
    if (!is.list(object$ev) |
        !all(names(object$ev) == c("tp", "tn", "fp", "fn", "thresholds"))){
        msg <- 'Not a valid entry for ev (missing entries or no list)!'
        errors <- c(errors, msg)
    }
    # check that the lenghts of each entry in ev are the same
    if (length(unique(vapply(object$ev, length, FUN.VALUE = numeric(1)))) != 1){
        msg <- 'No concordance for the entries in ev (unequal length)!'
        errors <- c(errors, msg)
    }
    # check concordance between ev and prc
    if (length(object$prc$recall) != length(object$ev$thresholds)){
        msg <- 'No concordance for the entries in ev and prc (unequal length)!'
        errors <- c(errors, msg)
    }

    # for the case that there are multiple repeats
    if (!is.null(object$roc.all)){
        # check if all entries are there
        if (!all(c('roc.all', 'auroc.all', 'prc.all', 'auprc.all', 'ev.all')
            %in% names(object))){
                msg <- 'Not all needed entries are given!'
                errors <- c(errors, msg)
        }
        # test roc.all
        if (!all(vapply(object$roc.all, class,
            FUN.VALUE = character(1)) == 'roc')){
                msg <- 'roc.all entries are not objects of class roc!'
                errors <- c(errors, msg)
            }
        # test lenght concordance
        if (length(unique(
            vapply(object[grep('.all', names(object))], length,
            FUN.VALUE = integer(1))))!=1){
            msg<-'entries for individual repeats do not have concordant length!'
            errors <- c(errors, msg)
        }
        ### MORE CHECKS FOR EVAL_DATA WITH MULTIPLE REPEATS?
    }
    if (length(errors) == 0) NULL else errors
}

# ##############################################################################
# check siamcat object for validity
#'@keywords internal
check.siamcat <- function(object){
    errors <- character()
    temp <- check.label(label(object))
    if (!is.null(temp)){
        errors <- c(errors, temp)
    }
    if (!is.null(filt_feat(object, verbose=0))){
        temp <- check.filt.feat(filt_feat(object))
        if (!is.null(temp)){
            errors <- c(errors, temp)
        }
    }
    if (!is.null(associations(object, verbose=0))){
        temp <- check.assoc(object@associations)
        if (!is.null(temp)){
            errors <- c(errors, temp)
        }
    }
    if (!is.null(norm_feat(object, verbose=0))){
        temp <- check.norm.feat(norm_feat(object))
        if (!is.null(temp)){
            errors <- c(errors, temp)
        }
    }
    if (!is.null(data_split(object, verbose=0))){
        temp <- check.data.split(data_split(object))
        if (!is.null(temp)){
            errors <- c(errors, temp)
        }
    }
    if (!is.null(model_list(object, verbose=0))){
        temp <- check.model.list(model_list(object))
        if (!is.null(temp)){
            errors <- c(errors, temp)
        }
    }
    if (!is.null(eval_data(object, verbose=0))){
        temp <- check.eval.data(eval_data(object))
        if (!is.null(temp)){
            errors <- c(errors, temp)
        }
    }
    if (length(errors) == 0) TRUE else errors
}


# ##############################################################################
# SIAMCAT object

#' The S4 SIAMCAT class
#'
#' @name siamcat-class
#' @rdname siamcat-class
#'
#' @description The SIAMCAT class
#'
#' @details The S4 SIAMCAT class stores the results from the SIAMCAT
#' workflow in different slots. The different slots will be filled by
#' different functions (referenced in the description below).
#'
#' In order to contruct a SIAMCAT class object, please refer to the
#' documentation of the construction function \link{siamcat}.
#'
#' The SIAMCAT class is based on the \link[phyloseq]{phyloseq-class}. Therefore,
#' you can easily import a \code{phyloseq} object into SIAMCAT.
#'
#'
#' @slot phyloseq object of class \link[phyloseq]{phyloseq-class}
#'
#' @slot label list containing the label information for the samples and
#' some metadata about the label, created by \link{create.label} or when
#' creating the \link{siamcat-class} object by calling \link{siamcat}
#'
#' @slot filt_feat list containing the filtered features as matrix and
#' the list of filtering parameters, created by calling the
#' \link{filter.features} function
#'
#' @slot associations list containing the parameters for association
#' testing and the results of association testing with these parameters in
#' a dataframe, created by calling the \link{check.associations} function
#'
#' @slot norm_feat list containing the normalized features as matrix and
#' the list of normalziation parameters (for frozen normalization), created by
#' calling the \link{normalize.features} function
#'
#' @slot data_split list containing cross-validation instances, created by
#' calling the \link{create.data.split} function
#'
#' @slot model_list list containing the trained models, the type of model
#' that was trained, and on which kind of features it was trained, created by
#' calling the \link{train.model} function
#'
#' @slot pred_matrix matrix of predictions, created by calling the
#' \link{make.predictions} function
#'
#' @slot eval_data list containing different evaluation metrics, created by
#' calling the \link{evaluate.predictions} function
#'
#' @exportClass siamcat
setClass(
    "siamcat",
    representation(
        phyloseq = "phyloseq",
        label = "list",
        filt_feat = "list",
        associations = "list",
        norm_feat = "list",
        data_split = "list",
        model_list = "list",
        pred_matrix = "matrix",
        eval_data = "list"
    ),
    validity=check.siamcat
)
