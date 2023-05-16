#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Model training
#'
#' @description This function trains the a machine learning model on the
#' training data
#'
#' @usage train.model(siamcat, method = "lasso", measure = "classif.acc",
#' param.set = NULL, grid.size=11, min.nonzero=5, perform.fs = FALSE,
#' param.fs = list(no_features = 100, method = "AUC", direction="absolute"),
#' feature.type='normalized', verbose = 1)
#'
#' @param siamcat object of class \link{siamcat-class}
#'
#' @param method string, specifies the type of model to be trained, may be one
#' of these: \code{c('lasso', 'enet', 'ridge', 'lasso_ll', 'ridge_ll',
#' 'randomForest')}
#'
#' @param measure character, specifies the model selection criterion during
#' internal cross-validation, see \link[mlr3]{mlr_measures} for more details,
#' defaults to \code{'classif.acc'}
#'
#' @param param.set list, set of extra parameters for mlr, see below for
#' details, defaults to \code{NULL}
#'
#' @param grid.size integer, grid size for internal tuning (needed for some
#' machine learning methods, for example \code{lasso_ll}), defaults to
#' \code{11}
#'
#' @param min.nonzero integer number of minimum nonzero coefficients that
#' should be present in the model (only for \code{'lasso'}, \code{'ridge'},
#' and \code{'enet'}), defaults to \code{5}
#'
#' @param perform.fs boolean, should feature selection be performed? Defaults
#' to \code{FALSE}
#'
#' @param param.fs list, parameters for the feature selection, see Details,
#' defaults to \code{list(thres.fs=100, method.fs="AUC", direction='absolute')}
#'
#' @param feature.type string, on which type of features should the function
#' work? Can be either \code{"original"}, \code{"filtered"}, or
#' \code{"normalized"}. Please only change this paramter if you know what you
#' are doing!
#'
#' @param verbose integer, control output: \code{0} for no output at all,
#' \code{1} for only information about progress and success, \code{2} for
#' normal level of information and \code{3} for full debug information,
#' defaults to \code{1}
#'
#' @export
#'
#' @encoding UTF-8
#'
#' @keywords SIAMCAT plm.trainer
#'
#' @return object of class \link{siamcat-class} with added \code{model_list}
#'
#'
#'
#' @section Machine learning methods:
#' This functions performs the training of the machine learning model
#' and functions as an interface to the \code{mlr3}-package.
#'
#' The function expects a \link{siamcat-class}-object with a prepared
#' cross-validation (see \link{create.data.split}) in the
#' \code{data_split}-slot of the object. It then trains a model for each fold
#' of the data split.
#'
#' The different machine learning methods are implemented as Learners from the
#' \link{mlr3learners} package:
#' \itemize{
#' \item \code{'lasso'}, \code{'enet'}, and \code{'ridge'} use the
#' \code{'classif.cv_glmnet'} or \code{'regr.cv_glmnet'} Learners, which
#' interface to the \link{glmnet} package,
#' \item \code{'lasso_ll'} and \code{'ridge_ll'} use a custom Learner, which
#' is only available for classification tasks. The underlying package is the
#' \link{LiblineaR} packge.
#' \item \code{'randomForest'} is implemented via the \code{'classif.ranger'}
#' or \code{regr.ranger} Learners available trough the \link[ranger]{ranger}
#' package.}
#'
#' @section Hyperparameter tuning:
#' There is additional control over the machine learning procedure by
#' supplying information through the \code{param.set} parameter within the
#' function. We encourage you to check out the excellent
#' \href{https://mlr3book.mlr-org.com/optimization.html}{mlr documentation}
#' for more in-depth information.
#'
#' Here is a short overview which parameters you can supply in which form:
#' \itemize{
#' \item enet The \strong{alpha} parameter describes the mixture between
#' lasso and ridge penalty and is -per default- determined using internal
#' cross-validation (the default would be equivalent to
#' \code{param.set=list('alpha'=c(0,1))}). You can supply either the limits of
#' the hyperparameter exploration (e.g. with limits 0.2 and 0.8:
#' \code{param.set=list('alpha'=c(0.2,0.8))}) or you can supply a fixed alpha
#' value as well (\code{param.set=list('alpha'=0.5)}).
#' \item lasso_ll/ridge_ll You can supply both \strong{class.weights} and
#' the \strong{cost} parameter (cost of the constraints violation, see
#' \link[LiblineaR]{LiblineaR} for more info). The default values would be
#' equal to \code{param.set=list('class.weights'=c(5, 1),
#' 'cost'=c(-2, 3))}.
#' \item randomForest You can supply the two parameters \strong{num.trees}
#' (Number of trees to grow) and \strong{mtry} (Number of variables randomly
#' sampled as candidates at each split). See also
#' \link[ranger]{ranger} for more info. The default values
#' correspond to
#' \code{param.set=list('num.trees'=c(100, 1000), 'mtry'=
#' c(round(sqrt.mdim / 2), round(sqrt.mdim), round(sqrt.mdim * 2)))} with
#' \code{sqrt.mdim=sqrt(nrow(data))}.
#' }
#'
#' @section Feature selection:
#' If feature selection should be performed (for example for functional data
#' with a large number of features), the \code{param.fs} list should contain:
#' \itemize{ \item \code{no_features} - Number of features to be retained after
#' feature selection,
#' \item \code{method} - method for the feature selection, may be
#' \code{AUC}, \code{gFC}, or \code{Wilcoxon} for binary classification
#' problems or \code{spearman}, \code{pearson}, or \code{MI} (mutual
#' information) for regression problems
#' \item \code{direction} - indicates if the feature selection should be
#' performed in a single direction only. Can be either \itemize{
#' \item \code{absolute} -  select the top associated features (independent of
#' the sign of enrichment),
#' \item \code{positive}the top positively associated featured (enriched in
#' the case group for binary classification or enriched in higher values
#' for regression),
#' \item \code{negative} the top negatively associated features (inverse of
#' positive)} Direction will be ignored for \code{Wilcoxon} and \code{MI}.}
#'
#' @examples
#' data(siamcat_example)
#'
#' # simple working example
#' siamcat_example <- train.model(siamcat_example, method='lasso')
train.model <- function(siamcat, method = "lasso",
    measure = "classif.acc", param.set = NULL,
    grid.size=11, min.nonzero=5, perform.fs = FALSE,
    param.fs = list(no_features = 100, method = "AUC", direction="absolute"),
    feature.type='normalized', verbose = 1) {

    if (verbose > 1)
        message("+ starting train.model")
    s.time <- proc.time()[3]

    # check and get features
    if (!feature.type %in% c('original', 'filtered', 'normalized')){
        stop("Unrecognised feature type, exiting...\n")
    }
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

    # make sure the names fit
    rownames(feat) <- make.names(rownames(feat))

    # checks
    label <- label(siamcat)
    if (label$type == "TEST"){
        stop('Model can not be trained to SIAMCAT object with a TEST label.',
            ' Exiting...')
    }
    if (is.null(data_split(siamcat, verbose=0))){
        stop("SIAMCAT object needs a data split for model training! Exiting...")
    }
    data.split <- data_split(siamcat)

    # Create List to save models.
    models.list <- list()
    num.runs <- data.split$num.folds * data.split$num.resample

    ## Create Learner
    if (perform.fs){
        if (!'no_features' %in% names(param.fs)){
            stop("Parameter '' is missing in `param.fs`!")
        }
        lrn <- create.mlr.learner(method, param.fs$no_features, param.set,
                                    type=label$type)
    } else {
        lrn <- create.mlr.learner(method, nrow(feat), param.set,
                                    type=label$type)
    }

    get_logger("mlr3")$set_threshold('off')
    get_logger("bbotk")$set_threshold('off')
    # loop over the folds
    bar <- 0
    if (verbose > 1){
        msg <- paste("+ training", method, "models on", num.runs,
            "training sets")
        message(msg)
    }
    if (verbose > 1 & perform.fs){
        message('+ Performing feature selection ',
                'with following parameters:')
        for (i in seq_along(param.fs)) {
            msg <- paste0('    ', names(param.fs)[i], ' = ',
                ifelse(is.null(param.fs[[i]]), 'NULL', param.fs[[i]]))
            message(msg)
        }
    }

    if (verbose == 1 || verbose == 2)
        pb <- progress_bar$new(total = num.runs)

    for (fold in seq_len(data.split$num.folds)) {
        if (verbose > 2){
            msg <- paste("+++ training on cv fold:", fold)
            message(msg)
        }

        for (resampling in seq_len(data.split$num.resample)) {
            if (verbose > 2){
                msg <- paste("++++ repetition:", resampling)
                message(msg)
            }
            ## Prepare data
            fold.name <- paste0("cv_fold", as.character(fold), "_rep",
                as.character(resampling))
            fold.exm.idx <-
                match(data.split$training.folds[[resampling]][[fold]],
                    names(label$label))

            ### subselect examples for training
            label.fac <- label$label
            if (label$type == 'BINARY'){
                label.fac <- factor(label.fac,  levels = sort(label$info))
            }
            train.label <- label.fac[fold.exm.idx]
            data <-
                as.data.frame(t(feat)[fold.exm.idx,])
            stopifnot(nrow(data) == length(train.label))
            stopifnot(all(rownames(data) == names(train.label)))

            if (perform.fs){
                data <- perform.feature.selection(data, train.label,
                    param.fs, verbose)
            }

            data$label <- train.label


            if (label$type == 'BINARY' && data$label[1] != -1) {
                data <- data[c(which(data$label == -1)[1],
                    c(seq_len(nrow(data)))[-which(data$label == -1)[1]]), ]
            }


            ## create task
            if (label$type == 'BINARY'){
                task <- TaskClassif$new(id='classif', backend=data,
                                        target='label', positive="1")
            } else if (label$type == 'CONTINUOUS') {
                task <- TaskRegr$new(id='regr', backend=data, target='label')
            }
            lrn.fold <- lrn$clone(deep=TRUE)

            ## Train model
            any.tuner <- unlist(lapply(lrn$param_set$values, FUN=class))
            if (any(any.tuner=='TuneToken')){
                instance <- tune(tnr("grid_search", resolution=grid.size),
                    task = task,
                    learner = lrn.fold,
                    resampling = rsmp("cv", folds = 5),
                    measures = msr(measure))
                lrn.fold$param_set$values <- instance$result_learner_param_vals
            }
            model <- lrn.fold$train(task = task)
            if (method %in% c('lasso', 'enet', 'ridge')){
                model <- get.best.glmnet.lambda(model, measure,
                                                min.nonzero, task)
                s.idx <- which(model$model$lambda==model$param_set$values$s)
                feat.weights <- -model$model$glmnet.fit$beta[,s.idx]
            } else if (method %in% c('lasso_ll', 'ridge_ll')){
                feat.weights <- t(model$model$W)[,1]
                idx.intercept <- which(names(feat.weights) == 'Bias')
                feat.weights <- feat.weights[-idx.intercept]
            } else if (method == 'randomForest'){
                feat.weights <- model$importance()
            } else if (method == 'SVM'){
                feat.weights <- rep(NA_real_, ncol(model$model$SV))
                names(feat.weights) <- colnames(model$model$SV)
            } else if (method == 'LSVM'){
                feat.weights <- coef(model$model)
                feat.weights <- feat.weights[-1]
                names(feat.weights) <- colnames(model$model$SV)
            }

            # add feature weights to the model
            bar <- bar + 1

            models.list[[fold.name]] <- list(model=model, features=feat.weights)
            if (verbose == 1 || verbose == 2)
                pb$tick()
        }
    }

    model_list(siamcat) <- list(
        models = models.list,
        model.type = method,
        feature.type = feature.type)

    e.time <- proc.time()[3]

    if (verbose > 1){
        msg <- paste("+ finished train.model in", 
            formatC(e.time - s.time, digits = 3), "s")
        message(msg)
    }
    if (verbose == 1){
        msg <- paste("Trained", method, "models successfully.")
        message(msg)
    }

    return(siamcat)
}
