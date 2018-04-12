#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Evaluate prediction results
#' @description This function takes the correct labels and predictions for all
#'        samples and evaluates the results using the \itemize{
#'          \item Area Under the Receiver Operating Characteristic (ROC)
#'                Curve (AU-ROC)
#'          \item and the Precision-Recall Curve (PR)
#'        }
#'        as metric. Predictions can be supplied either for a single case or as
#'        matrix after resampling of the dataset.
#'
#'        Prediction results are usually produced with the function
#'        \link{make.predictions}.
#' @param siamcat object of class \link{siamcat-class}
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'        for only information about progress and success, \code{2} for normal
#'        level of information and \code{3} for full debug information,
#'        defaults to \code{1}
#' @keywords SIAMCAT evaluate.predictions
#' @details This functions calculates for the predictions in the
#'        \code{pred_matrix} -slot of the \link{siamcat-class}-object several
#'        metrices. The Area Under the Receiver Operating Characteristic (ROC)
#'        Curve (AU-ROC) and the Precision-Recall Curve will be evaluated and
#'        the results will be saved in the \code{eval_data}-slot of the
#'        supplied \link{siamcat-class}-object. The \code{eval_data}-slot
#'        contains a list with several entries: \itemize{
#'          \item \code{$roc.average} average ROC-curve across repeats or a
#'                single ROC-curve on complete dataset;
#'          \item \code{$auc.average} AUC value for the average ROC-curve;
#'          \item \code{$ev.list} list of \code{length(num.folds)}, containing
#'                for different decision thresholds the number of false
#'                positives, false negatives, true negatives, and
#'                true positives;
#'          \item \code{$pr.list} list of \code{length(num.folds)}, containing
#'                the positive predictive value (precision) and true positive
#'                rate (recall) values used to plot the PR curves.
#' }      For the case of repeated cross-validation, the function will
#'        additonally return \itemize{
#'          \item \code{$roc.all} list of roc objects (see \link[pROC]{roc})
#'                for every repeat;
#'          \item \code{$aucspr} vector of AUC values for the PR curves
#'                for every repeat;
#'          \item \code{$auc.all} vector of AUC values for the ROC curves
#'                for every repeat.
#'}
#' @export
#' @return object of class \link{siamcat-class} with the
#'         slot \code{eval_data} filled
#' @examples
#'
#'  data(siamcat_example)
#'  # simple working example
#'  siamcat_evaluated <- evaluate.predictions(siamcat_example)
#'
evaluate.predictions <- function(siamcat, verbose = 1) {
    if (verbose > 1)
        message("+ starting evaluate.predictions")
    label  <- get.label.list(siamcat)
    s.time <- proc.time()[3]
    # TODO compare header to label make sure that label and prediction are in the same order
    m <- match(names(label$label), rownames(pred_matrix(siamcat)))

    pred <- pred_matrix(siamcat)[m, , drop = FALSE]
    stopifnot(all(names(label$label) == rownames(pred)))

    # ROC curve
    if (verbose > 2)
        message("+ calculating ROC")
    auroc = 0
    if (ncol(pred) > 1) {
        rocc = list(NULL)
        aucs = vector("numeric", ncol(pred))
        for (c in seq_len(ncol(pred))) {
            rocc[c] = list(roc(response = label$label, predictor = pred[, c], ci = FALSE))
            aucs[c] = rocc[[c]]$auc
        }
        l.vec = rep(label$label, ncol(pred))
    } else {
        l.vec = label$label
    }
    # average data for plotting one mean prediction curve
    if (verbose > 2)
        message("+ calculating mean ROC")
    summ.stat = "mean"
    rocsumm = list(roc(response = label$label, predictor = apply(pred, 1, summ.stat), ci = TRUE, of = "se",
        sp = seq(0, 1, 0.05)))
    auroc = list(rocsumm[[1]]$auc)
    # precision recall curve
    pr = list(NULL)
    ev = list(NULL)
    if (ncol(pred) > 1) {
        aucspr = vector("numeric", dim(pred)[2])
        for (c in seq_len(ncol(pred))) {
            ev[c] = list(evaluate.classifier(pred[, c], label$label, label, verbose = verbose))
            pr[c] = list(evaluate.get.pr(ev[[c]], verbose = verbose))
            aucspr[c] = evaluate.calc.aupr(ev[[c]], verbose = verbose)
        }
        ev = append(ev, list(evaluate.classifier(apply(pred, 1, summ.stat), label$label, label)))
    } else {
        ev[1] = list(evaluate.classifier(as.vector(pred), label$label, label, verbose = verbose))
        pr[1] = list(evaluate.get.pr(ev[[1]]), verbose = verbose)
    }
    if (ncol(pred) > 1) {
        if (verbose > 2)
            message("+ evaluating multiple predictions")
        siamcat@eval_data <- list(roc.all = rocc, auc.all = aucs, roc.average = rocsumm, auc.average = auroc, ev.list = ev,
            pr.list = pr, aucspr = aucspr)
    } else {
        if (verbose > 2)
            message("+ evaluating single prediction")
        siamcat@eval_data <- list(roc.average = rocsumm, auc.average = auroc, ev.list = ev, pr.list = pr)
    }
    e.time <- proc.time()[3]
    if (verbose > 1)
        message(paste("+ finished evaluate.predictions in", formatC(e.time - s.time, digits=3), "s"))
    if (verbose == 1)
        message("Evaluated predictions successfully.")
    return(siamcat)
}

# evaluates the predictive performance of a classifier on a labeled data sets returns a list with vectors
# containing TP, FP, TN, FN for each threshold value on the predictions (where TP = true positives, FP = false
# positives, TN = true negatives, FN = false negatives)
#' @keywords internal
evaluate.classifier <- function(predictions, test.label, label, verbose = 0) {
    if (verbose > 2)
        message("+ starting evaluate.classifier")
    stopifnot(dim(test.label) == NULL)
    stopifnot(length(unique(test.label)) == 2)
    stopifnot(all(is.finite(predictions)))
    # calculate thresholds, one between each subsequent pair of sorted prediction values this is ignorant to whether
    # predictions is in matrix or vector format (see below)
    thr <- predictions
    dim(thr) <- NULL
    thr <- sort(unique(thr))
    thr <- rev(c(min(thr) - 1, (thr[-1] + thr[-length(thr)])/2, max(thr) + 1))
    if (is.null(dim(predictions))) {
        # assuming that a single model was applied to predict the data set
        stopifnot(length(test.label) == length(predictions))
        # actual evaluations per threshold value
        tp = vapply(thr, FUN=function(x){sum(test.label == label$positive.lab & predictions > x)}, USE.NAMES = FALSE, FUN.VALUE = integer(1))
        fp = vapply(thr, FUN=function(x){sum(test.label == label$negative.lab & predictions > x)}, USE.NAMES = FALSE, FUN.VALUE = integer(1))
        tn = vapply(thr, FUN=function(x){sum(test.label == label$negative.lab & predictions < x)}, USE.NAMES = FALSE, FUN.VALUE = integer(1))
        fn = vapply(thr, FUN=function(x){sum(test.label == label$positive.lab & predictions < x)}, USE.NAMES = FALSE, FUN.VALUE = integer(1))
    } else {
        # assuming that several models were applied to predict the same data and predictions of each model occupy one
        # column
        stopifnot(length(test.label) == nrow(predictions))
        tp = t(tp = vapply(thr, FUN=function(x){apply(predictions, 2, FUN=function(y){sum(test.label == label$positive.lab & y > x)})}, USE.NAMES = FALSE, FUN.VALUE = integer(2)))
        fp = t(tp = vapply(thr, FUN=function(x){apply(predictions, 2, FUN=function(y){sum(test.label == label$negative.lab & y > x)})}, USE.NAMES = FALSE, FUN.VALUE = integer(2)))
        tn = t(tp = vapply(thr, FUN=function(x){apply(predictions, 2, FUN=function(y){sum(test.label == label$negative.lab & y < x)})}, USE.NAMES = FALSE, FUN.VALUE = integer(2)))
        fn = t(tp = vapply(thr, FUN=function(x){apply(predictions, 2, FUN=function(y){sum(test.label == label$positive.lab & y < x)})}, USE.NAMES = FALSE, FUN.VALUE = integer(2)))
    }
    if (verbose > 2)
        message("+ finished evaluate.classifier")
    return(list(tp = tp, tn = tn, fp = fp, fn = fn, thresholds = thr))
}

# calculates the area under a curve using a trapezoid approximation
evaluate.area.trapez = function(x, y, verbose = 0) {
    if (verbose > 2)
        message("+ starting evaluate.area.trapez")
    if (x[1] > x[length(x)]) {
        x = rev(x)
        y = rev(y)
    }
    xd = x[-1] - x[-length(x)]
    ym = 0.5 * (y[-1] + y[-length(y)])
    if (verbose > 2)
        message("+ finished evaluate.area.trapez")
    return(xd %*% ym)
}

# returns a vector of x and y values for plotting a precision-recall curve
evaluate.get.pr = function(eval, verbose = 0) {
    if (verbose > 2)
        message("+ starting evaluate.get.pr")
    tpr = eval$tp/(eval$tp + eval$fn)
    ppv = eval$tp/(eval$tp + eval$fp)
    # at thresholds where the classifier makes no positive predictions at all, we (somewhat arbitrarily) set its
    # precision to 1
    ppv[is.na(ppv)] = 1
    if (verbose > 2)
        message("+ finished evaluate.get.pr")
    return(list(x = tpr, y = ppv))
}

# calculates the area under the precision-recall curve (over the interval [0, max.tpr], if specified)
evaluate.calc.aupr = function(eval, max.tpr = 1, verbose = 0) {
    if (verbose > 2)
        message("+ starting evaluate.calc.aupr")
    pr = evaluate.get.pr(eval, verbose = verbose)
    idx = pr$x <= max.tpr
    if (verbose > 2)
        message("+ finished evaluate.calc.aupr")
    return(evaluate.area.trapez(pr$x[idx], pr$y[idx]))
}
