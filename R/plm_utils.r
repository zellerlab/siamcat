#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0


##### Internal function to create the mlr3 learner with all kinds of parameters
#' @keywords internal
create.mlr.learner <- function(method, nrow.data, param.set=NULL,
                                type='BINARY'){
    if (!method %in% c("lasso", "enet", "ridge", "lasso_ll",
            "ridge_ll", "randomForest", "SVM", "LSVM")){
        stop("Unsupported method!")
    }
    standard.param.set <- list(
        "cost" = c(-2, 3),
        "epsilon" = 1e-08,
        "ntree" = c(100, 1000),
        "mtry" = c(round(sqrt(nrow.data) / 2),
            round(sqrt(nrow.data)),
            round(sqrt(nrow.data) * 2)),
        "alpha" = c(0, 1),
        "class.weights" = c("-1"=5, "1"=1))

    if (is.null(param.set)){
        use.param.set <- standard.param.set
    } else {
        use.param.set <- param.set
        for (i in names(standard.param.set)){
            if (!i %in% names(param.set)){
                use.param.set[[i]] <- standard.param.set[[i]]
            }
        }
    }

    # lasso/enet/ridge
    if (method %in% c('lasso', 'enet', 'ridge')){
        if (type == 'BINARY'){
            learner <- lrn("classif.cv_glmnet")
            learner$predict_type <- 'prob'
        } else if (type == 'CONTINUOUS'){
            learner <- lrn("regr.cv_glmnet")
        }
        if (method == 'lasso'){
            learner$param_set$values$alpha <- 1
            if ('alpha' %in% names(param.set)){
                warning("Parameter 'alpha' will be ignored and set to 1!")
            }
        } else if (method == 'ridge'){
            learner$param_set$values$alpha <- 0
            if ('alpha' %in% names(param.set)){
                warning("Parameter 'alpha' will be ignored and set to 0!")
            }
        } else if (method == 'enet'){
            if (length(use.param.set$alpha)==1){
                learner$param_set$values$alpha <- use.param.set$alpha
            } else {
                learner$param_set$values$alpha <- to_tune(
                    lower=use.param.set$alpha[1],
                    upper=use.param.set$alpha[2])
            }
        }
        use.param.set$alpha <- NULL
    } else if (method=='randomForest'){
        if (type == 'BINARY'){
            learner <- lrn("classif.ranger", importance='permutation')
            learner$predict_type <- 'prob'
        } else if (type == 'CONTINUOUS'){
            learner <- lrn("regr.ranger", importance='impurity')
        }
        # number of trees
        if (length(use.param.set$ntree) == 1){
            learner$param_set$values$num.trees <- use.param.set$ntree
        } else if (length(use.param.set$ntree) == 2) {
            learner$param_set$values$num.trees <- to_tune(
                lower=use.param.set$ntree[1], upper=use.param.set$ntree[2])
        } else {
            learner$param_set$values$num.trees <- to_tune(use.param.set$ntree)
        }
        # mtry
        if (length(use.param.set$mtry) == 1){
            learner$param_set$values$mtry <- use.param.set$mtry
        } else if (length(use.param.set$mtry) == 2){
            learner$param_set$values$mtry <- to_tune(
                lower=use.param.set$mtry[1], upper=use.param.set$mtry[2])
        } else {
            learner$param_set$values$mtry <- to_tune(use.param.set$mtry)
        }
        use.param.set$ntree <- NULL
        use.param.set$mtry <- NULL
        if (!'alpha' %in% names(param.set)){
            use.param.set$alpha <- NULL
        }
        if (!'class.weights' %in% names(param.set)){
            use.param.set$class.weights <- NULL
        }

    } else if (method %in% c('lasso_ll', 'ridge_ll')){
        if (type == 'BINARY'){
            if (method == 'lasso_ll'){
                type <- 6
            } else {
                type <- 0
            }
          mlr_learners$add("classif.liblinear", 
                           LearnerClassifLiblineaR)
            learner <- lrn('classif.liblinear', type=type,
                wi=use.param.set$class.weights,
                epsilon=use.param.set$epsilon)
            learner$predict_type <- 'prob'
        } else if (type == 'CONTINUOUS'){
            stop("Methods not usable for regression tasks!")
        }
        if (length(use.param.set$cost) == 1){
            learner$param_set$values$cost <- use.param.set$cost
        } else if (length(use.param.set$cost) == 2){
            learner$param_set$values$cost <- to_tune(p_dbl(
                lower=use.param.set$cost[1],
                upper=use.param.set$cost[2], trafo=.f_exp))
        } else {
            learner$param_set$values$cost <- to_tune(use.param.set$cost)
        }
        use.param.set$class.weights <- NULL
        use.param.set$epsilon <- NULL
        use.param.set$cost <- NULL
    } else if (method == 'SVM'){
        if (type == 'BINARY'){
            learner <- lrn('classif.svm', kernel='radial')
            learner$predict_type <- 'prob'
            learner$param_set$values$type <- 'C-classification'
            learner$param_set$values$cost <- to_tune(p_dbl(
                lower=-5,
                upper=16, trafo=.f_exp))
            use.param.set$cost <- NULL
            use.param.set$class.weights <- NULL
        } else if (type == 'CONTINUOUS'){
            stop("Methods not usable for regression tasks!")
        }
    } else if (method == 'LSVM'){
        if (type == 'BINARY'){
            learner <- lrn('classif.svm', kernel='linear')
            learner$predict_type <- 'prob'
            learner$param_set$values$type <- 'C-classification'
            learner$param_set$values$cost <- to_tune(p_dbl(
                lower=-5,
                upper=16, trafo=.f_exp))
            use.param.set$cost <- NULL
            use.param.set$class.weights <- NULL
        } else if (type == 'CONTINUOUS'){
            stop("Methods not usable for regression tasks!")
        }
    }
    
    # try to set additional parameters, i hope that mlr catches errors here
    param.settable <- learner$param_set$ids()
    for (x in intersect(names(use.param.set), param.settable)){
        if (length(use.param.set[[x]])==1){
            learner$param_set$values[[x]] <- use.param.set[[x]]
        } else {
            learner$param_set$values[[x]] <- to_tune(use.param.set[[x]])
        }
    }

    return(learner)

}

# function to perform feature selection
#' @keywords internal
perform.feature.selection <- function(data, train.label, param.fs, verbose){

    stopifnot(all(c('method', 'no_features', 'direction') %in%names(param.fs)))

    # test method.fs
    if (is.factor(train.label)){
        allowed.methods <- c('Wilcoxon', 'AUC', 'gFC')
        if (!param.fs$method %in% allowed.methods) {
            stop('Unrecognised feature selection method. ',
                    'Must be one of those: {"',
                    paste(allowed.methods, collapse = '", "'), '"}')
        }
    } else if (is.numeric(train.label)){
        allowed.methods <- c('spearman', 'pearson', 'MI')
        if (!param.fs$method %in% allowed.methods) {
            stop('Unrecognised feature selection method. ',
                    'Must be one of those: {"',
                    paste(allowed.methods, collapse = '", "'), '"}')
        }
    }

    # assert the threshold
    stopifnot(param.fs$no_features > 10)
    stopifnot(param.fs$no_features < ncol(data))

    if (param.fs$method == 'Wilcoxon') {
        assoc <- vapply(data,
                        FUN=function(x, label){
                            d <- data.frame(x=x, y=label);
                            t <- wilcox.test(x~y, data=d)
                            return(t$p.val)
                        }, FUN.VALUE=double(1),
                        label=train.label)
        assoc <- sort(assoc)
    } else if (param.fs$method == 'AUC') {
        assoc <- vapply(data,
                        FUN=get.single.feat.AUC,
                        FUN.VALUE = double(1),
                        label=train.label,
                        pos=max(levels(train.label)),
                        neg=min(levels(train.label)))
        if (param.fs$direction == 'absolute'){
            assoc[assoc < 0.5] <- 1 - assoc[assoc < 0.5]
        } else if (param.fs$direction == 'negative'){
            assoc <- 1 - assoc
        }
        assoc <- assoc[assoc > 0.5]
        assoc <- sort(-assoc)
    } else if (param.fs$method == 'gFC') {
        assoc <- vapply(data,
                        FUN=get.quantile.FC,
                        FUN.VALUE = double(1),
                        label=train.label,
                        pos=max(levels(train.label)),
                        neg=min(levels(train.label)))
        if (param.fs$direction == 'absolute'){
            assoc <- abs(assoc)
        } else if (param.fs$direction == 'negative'){
            assoc <- -assoc
        }
        assoc <- assoc[assoc > 0]
        assoc <- sort(-assoc)
    } else if (param.fs$method %in% c('spearman', 'pearson')){
        assoc <- vapply(data, FUN=cor, FUN.VALUE = double(1), y=train.label,
                        method=param.fs$method)
        if (param.fs$direction == 'absolute'){
            assoc <- abs(assoc)
        } else if (param.fs$direction == 'negative'){
            assoc <- -assoc
        }
        assoc <- assoc[assoc > 0]
        assoc <- sort(-assoc)
    } else if (param.fs$method == 'MI'){
        assoc <- vapply(data, FUN=function(x){
            mutinformation(discretize(x, disc='equalwidth'),
                discretize(train.label, disc='equalwidth'))
        }, FUN.VALUE = double(1))
        assoc <- sort(-assoc)
    }

    data <- data[,names(assoc)[seq_len(param.fs$no_features)]]

    stopifnot(ncol(data) > 0)
    if (verbose > 2) {
        message(paste0('++ retaining ', ncol(data),
                ' features after selection based on ',
                param.fs$method, '; target number of features ',
                param.fs$no_features))
    }
    return(data)
}

#' @keywords internal
measureAUPRC <- function(probs, truth, negative, positive) {
    pr <- pr.curve(scores.class0 = probs[which(truth == positive)],
            scores.class1 = probs[which(truth == negative)])
    return(pr$auc.integral)
}

#' @keywords internal
get.single.feat.AUC <- function(x, label, pos, neg) {
    x.p <- x[label == pos]
    x.n <- x[label == neg]
    temp.auc <- roc(cases=x.p, controls=x.n, direction='<')$auc
    return(temp.auc)
}

#' @keywords internal
get.quantile.FC <- function(x, label, pos, neg){
    x.p <- x[label == pos]
    x.n <- x[label == neg]
    q.p <- quantile(x.p, probs=seq(.1, .9, length.out=9))
    q.n <- quantile(x.n, probs=seq(.1, .9, length.out=9))
    return(sum(q.p - q.n)/length(q.p))
}

#' @keywords internal
.f_exp <- function(x){10^x}

#' @keywords internal
get.best.glmnet.lambda <- function(model, measure, min.nonzero, task){
    idx <- which(model$model$nzero >= min.nonzero)
    new.model <- model$clone(deep = TRUE)
    perf <- vapply(idx, FUN=function(x){
        new.model$param_set$values$s <- new.model$model$lambda[x]
        pred <- new.model$predict(task)
        pred$score(msr(measure))
    }, FUN.VALUE = double(1))
    if (msr(measure)$minimize){
        f.idx <- which(perf==min(perf))[1]
    } else {
        f.idx <- which(perf==max(perf))[1]
    }
    new.model$param_set$values$s <- new.model$model$lambda[f.idx]
    return(new.model)
}
