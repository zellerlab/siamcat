#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Evaluate prediction results
#'
#' @description This function compares the predictions (from 
#' [make.predictions]) and true labels for all samples and evaluates 
#' the results. 
#'
#' @param siamcat object of class \link{siamcat-class}
#'
#' @param verbose integer, control output: \code{0} for no output at all, 
#' \code{1} for only information about progress and success, \code{2} for 
#' normal level of information and \code{3} for full debug information,
#' defaults to \code{1}
#'
#' @keywords SIAMCAT evaluate.predictions
#'
#' @section Binary classification problems:
#' This function calculates several metrics for the predictions in 
#' the \code{pred_matrix}-slot of the \link{siamcat-class}-object. 
#' The Area Under the Receiver Operating Characteristic (ROC) Curve (AU-ROC) 
#' and the Precision-Recall Curve will be evaluated and the results will be 
#' saved in the \code{eval_data}-slot of the supplied \link{siamcat-class}-
#' object. The \code{eval_data}-slot contains a list with several entries: 
#' \itemize{
#' \item \code{$roc} - average ROC-curve across repeats or a single ROC-curve 
#' on complete dataset (see \link[pROC]{roc});
#' \item \code{$auroc} - AUC value for the average ROC-curve;
#' \item \code{$prc} - list containing the positive predictive value 
#' (precision) and true positive rate (recall) values used to plot the mean 
#' PR curve;
#' \item \code{$auprc} - AUC value for the mean PR curve;
#' \item \code{$ev} - list containing for different decision thresholds the 
#' number of false positives, false negatives, true negatives, and true 
#' positives.}
#' For the case of repeated cross-validation, the function will additionally 
#' return \itemize{
#' \item \code{$roc.all} - list of roc objects (see \link[pROC]{roc}) 
#' for every repeat;
#' \item \code{$auroc.all} - vector of AUC values for the ROC curves 
#' for every repeat;
#' \item \code{$prc.all} - list of PR curves for every repeat;
#' \item \code{$auprc.all} - vector of AUC values for the PR curves 
#' for every repeat;
#' \item \code{$ev.all} - list of \code{ev} lists (see above) 
#' for every repeat.}
#'
#' @section Regression problems:
#' This function calculates several metrics for the evaluation of predictions
#' and will store the results in the \code{eval_data}-slot of the supplied 
#' \link{siamcat-class} objects. The \code{eval_data}-slot will contain:
#' \itemize{
#' \item \code{r2} - the mean R squared value across repeats or a single 
#' R-squared value on the complete dataset;
#' \item \code{mae} - them mean absolute error of the predictions;
#' \item \code{mse} - the mean squared error of the predictions.}
#' For the case of repeated cross-validation, the function will additionally 
#' compute all three of these measures for the individual cross-validation 
#' repeats and will store the results in the \code{eval_data} slot as 
#' \code{r2.all}, \code{mae.all}, and \code{mse.all}.
#'
#' @export
#'
#' @encoding UTF-8
#'
#' @return object of class \link{siamcat-class} with the slot 
#' \code{eval_data} filled
#'
#' @examples
#' data(siamcat_example)
#'
#' siamcat_evaluated <- evaluate.predictions(siamcat_example)
evaluate.predictions <- function(siamcat, verbose = 1) {
    if (verbose > 1)
        message("+ starting evaluate.predictions")
    s.time <- proc.time()[3]

    # check that predictions are present
    if (is.null(pred_matrix(siamcat, verbose=0))){
        stop('SIAMCAT object does not contain predictions. Exiting\n')
    }

    label  <- label(siamcat)
    if (label$type == 'TEST'){
        stop('SIAMCAT can not evaluate the predictions on a TEST',
            ' label. Exiting...')
    }

    if (label$type == 'BINARY'){
        r.object <- eval.binary(siamcat, s.time, verbose)
    } else if (label$type == 'CONTINUOUS'){
        r.object <- eval.regr(siamcat, s.time, verbose)
    }

    return(r.object)
}


#' @keywords internal
eval.regr <- function(siamcat, s.time, verbose=0){
    label  <- label(siamcat)
    m <- match(names(label$label), rownames(pred_matrix(siamcat)))

    pred <- pred_matrix(siamcat)[m, , drop = FALSE]
    stopifnot(all(names(label$label) == rownames(pred)))

    if (ncol(pred) > 1){
        if (verbose > 2)
            message("+ evaluating multiple predictions")

        # mean predictions
        pred.mean <- rowMeans(pred)
        ess <- sum((label$label - mean(label$label))^2)
        rss <- sum((pred.mean - label$label)^2)
        r2.mean <- 1 - rss/ess
        mae.mean <- mean(abs(pred.mean - label$label))
        mse.mean <- mean((abs(pred.mean - label$label))^2)

        # each repetition individually
        r2.all <- list()
        mae.all <- list()
        mse.all <- list()
        for (i in seq_len(ncol(pred))){
            rss <- sum((abs(pred[,i]) - label$label)^2)
            r2.all[[i]] <- 1 - rss/ess
            mae.all[[i]] <- mean(abs(pred[,i] - label$label))
            mse.all[[i]] <- mean((abs(pred[,i] - label$label))^2)
        }

        eval_data(siamcat) <- list(
            r2 = r2.mean, r2.all = r2.all,
            mae = mae.mean, mae.all = mae.all,
            mse = mse.mean, mse.all = mse.all)

    } else {
        if (verbose > 2)
            message("+ evaluating single prediction")
        rss <- sum((abs(pred[,1]) - label$label)^2)
        ess <- sum((label$label - mean(label$label))^2)
        r2 <- 1 - rss/ess
        mae <- mean(abs(pred[,1] - label$label))
        mse <- mean((abs(pred[,1] - label$label))^2)
        eval_data(siamcat) <- list(r2=r2, mae=mae, mse=mse)
    }
    e.time <- proc.time()[3]
    if (verbose > 1){
        msg <- paste(
            "+ finished evaluate.predictions in",
            formatC(e.time - s.time, digits = 3),
            "s")
        message(msg)
    }
    if (verbose == 1)
        message("Evaluated predictions successfully.")
    return(siamcat)
}

#' @keywords internal
eval.binary <- function(siamcat, s.time, verbose=0){
    label  <- label(siamcat)

    summ.stat <- "mean" # TODO make this a possible parameter?
    # TODO compare header to label make sure that label and prediction are in
    # the same order
    m <- match(names(label$label), rownames(pred_matrix(siamcat)))

    pred <- pred_matrix(siamcat)[m, , drop = FALSE]
    stopifnot(all(names(label$label) == rownames(pred)))

    # ##########################################################################
    # ROC curve
    if (verbose > 2)
        message("+ calculating ROC")
    auroc <- 0
    if (ncol(pred) > 1) {
        roc.all <- list()
        auroc.all <- vector("numeric", ncol(pred))
        for (c in seq_len(ncol(pred))) {
            roc.all[[c]] <- roc(response = label$label, predictor = pred[, c],
                        direction = '<', levels = label$info, ci = FALSE)
            auroc.all[c] <- roc.all[[c]]$auc
        }
        l.vec <- rep(label$label, ncol(pred))
    } else {
        l.vec <- label$label
    }

    # average data for plotting one mean prediction curve

    roc.mean <- roc(response = label$label,
                    predictor = apply(pred, 1, summ.stat),
                    ci = TRUE, of = "se",
                    sp = seq(0, 1, 0.05), direction = '<', levels = label$info)
    auroc <- roc.mean$auc

    # ##########################################################################
    # PR curve
    prc <- list()
    ev <- list()
    auprc <- 0
    if (ncol(pred) > 1) {
        auprc.all <- vector("numeric", ncol(pred))
        prc.all <- list()
        ev.all <- list()
        for (c in seq_len(ncol(pred))) {
            ev.all[[c]] <- evaluate.classifier(pred[, c], label$label, label,
                                            verbose = verbose)
            prc.all[[c]] <- evaluate.get.pr(ev.all[[c]], verbose = verbose)
            auprc.all[c] <- evaluate.calc.aupr(ev.all[[c]], verbose = verbose)
        }
        ev <- evaluate.classifier(apply(pred, 1, summ.stat), label$label, label)
    } else {
        ev <- evaluate.classifier(as.vector(pred), label$label, label,
                                verbose = verbose)
    }

    prc <- evaluate.get.pr(ev, verbose = verbose)
    auprc <- c(evaluate.calc.aupr(ev, verbose = verbose))

    if (ncol(pred) > 1) {
        if (verbose > 2)
            message("+ evaluating multiple predictions")
        eval_data(siamcat) <- list(
            roc= roc.mean, roc.all = roc.all,
            auroc = auroc, auroc.all = auroc.all,
            prc = prc, prc.all = prc.all,
            auprc = auprc, auprc.all = auprc.all,
            ev = ev, ev.all = ev.all
        )

    } else {
        if (verbose > 2)
            message("+ evaluating single prediction")
        eval_data(siamcat) <- list(
            roc=roc.mean, auroc=auroc, prc=prc, auprc=auprc, ev=ev
        )
    }
    e.time <- proc.time()[3]
    if (verbose > 1){
        msg <- paste(
            "+ finished evaluate.predictions in",
            formatC(e.time - s.time, digits = 3),
            "s")
        message(msg)
    }
    if (verbose == 1)
        message("Evaluated predictions successfully.")
    return(siamcat)
}

# evaluates the predictive performance of a classifier on a labeled data sets
# returns a list with vectors containing TP, FP, TN, FN for each threshold
# value on the predictions (where TP = true positives, FP = false positives,
# TN = true negatives, FN = false negatives)
#' @keywords internal
evaluate.classifier <-
    function(predictions, test.label, label, verbose = 0) {
        if (verbose > 2)
            message("+ starting evaluate.classifier")
        stopifnot(is.null(dim(test.label)))
        stopifnot(length(unique(test.label)) == 2)
        stopifnot(all(is.finite(predictions)))
        # calculate thresholds, one between each subsequent pair of sorted
        #prediction values this is ignorant to whether
        # predictions is in matrix or vector format (see below)
        thr <- predictions
        dim(thr) <- NULL
        thr <- sort(unique(thr))
        thr <-
            rev(c(min(thr) - 1, (thr[-1] + thr[-length(thr)]) / 2,
                max(thr) + 1))
        if (is.null(dim(predictions))) {
            # assuming that a single model was applied to predict the data set
            stopifnot(length(test.label) == length(predictions))
            # actual evaluations per threshold value
            tp <- vapply(
                thr,
                FUN = function(x) {
                    sum(test.label == max(label$info)
                        & predictions > x)
                },
                USE.NAMES = FALSE,
                FUN.VALUE = integer(1)
            )
            fp <- vapply(
                thr,
                FUN = function(x) {
                    sum(test.label == min(label$info)
                        & predictions > x)
                },
                USE.NAMES = FALSE,
                FUN.VALUE = integer(1)
            )
            tn <- vapply(
                thr,
                FUN = function(x) {
                    sum(test.label == min(label$info)
                        & predictions < x)
                },
                USE.NAMES = FALSE,
                FUN.VALUE = integer(1)
            )
            fn <- vapply(
                thr,
                FUN = function(x) {
                    sum(test.label == max(label$info)
                        & predictions < x)
                },
                USE.NAMES = FALSE,
                FUN.VALUE = integer(1)
            )
        } else {
            # assuming that several models were applied to predict the same data
            # and predictions of each model occupy one
            # column
            stopifnot(length(test.label) == nrow(predictions))
            tp <- t(vapply(
                thr,
                FUN = function(x) {
                    apply(
                        predictions,
                        2,
                        FUN = function(y) {
                            sum(test.label == max(label$info) & y > x)
                        }
                    )
                },
                USE.NAMES = FALSE,
                FUN.VALUE = integer(2)
            ))
            fp <- t(vapply(
                thr,
                FUN = function(x) {
                    apply(
                        predictions,
                        2,
                        FUN = function(y) {
                            sum(test.label == min(label$info) & y > x)
                        }
                    )
                },
                USE.NAMES = FALSE,
                FUN.VALUE = integer(2)
            ))
            tn <- t(vapply(
                thr,
                FUN = function(x) {
                    apply(
                        predictions,
                        2,
                        FUN = function(y) {
                            sum(test.label == min(label$info) & y < x)
                        }
                    )
                },
                USE.NAMES = FALSE,
                FUN.VALUE = integer(2)
            ))
            fn <- t(vapply(
                thr,
                FUN = function(x) {
                    apply(
                        predictions,
                        2,
                        FUN = function(y) {
                            sum(test.label == max(label$info) & y < x)
                        }
                    )
                },
                USE.NAMES = FALSE,
                FUN.VALUE = integer(2)
            ))
        }
        if (verbose > 2)
            message("+ finished evaluate.classifier")
        return(list(
            tp = tp,
            tn = tn,
            fp = fp,
            fn = fn,
            thresholds = thr
        ))
    }

# calculates the area under a curve using a trapezoid approximation
#' @keywords internal
evaluate.area.trapez <- function(x, y, verbose = 0) {
    if (verbose > 2)
        message("+ starting evaluate.area.trapez")
    if (x[1] > x[length(x)]) {
        x <- rev(x)
        y <- rev(y)
    }
    xd <- x[-1] - x[-length(x)]
    ym <- 0.5 * (y[-1] + y[-length(y)])
    if (verbose > 2)
        message("+ finished evaluate.area.trapez")
    return(xd %*% ym)
}

# returns a vector of x and y values for plotting a precision-recall curve
#' @keywords internal
evaluate.get.pr <- function(eval, verbose = 0) {
    if (verbose > 2)
        message("+ starting evaluate.get.pr")
    tpr <- eval$tp / (eval$tp + eval$fn)
    ppv <- eval$tp / (eval$tp + eval$fp)
    # at thresholds where the classifier makes no positive predictions at all,
    # we (somewhat arbitrarily) set its
    # precision to 1
    ppv[is.na(ppv)] <- 1
    if (verbose > 2)
        message("+ finished evaluate.get.pr")
    return(list(recall = tpr, precision = ppv))
}

# calculates the area under the precision-recall curve (over the interval
# [0, max.tpr], if specified)
#' @keywords internal
evaluate.calc.aupr <- function(eval,
    max.tpr = 1,
    verbose = 0) {
    if (verbose > 2)
        message("+ starting evaluate.calc.aupr")
    pr <- evaluate.get.pr(eval, verbose = verbose)
    idx <- pr$recall <= max.tpr
    if (verbose > 2)
        message("+ finished evaluate.calc.aupr")
    return(evaluate.area.trapez(pr$recall[idx], pr$precision[idx]))
}
