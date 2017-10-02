###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# RScript flavor
#
# written by Georg Zeller
# with additions by Nicolai Karcher and Konrad Zych
# EMBL Heidelberg 2012-2017
#
# version 0.2.0
# file last updated: 26.06.2017
# GNU GPL 3.0
###

#' @title Evaluate prediction results
#' @description This function takes the correct labels and predictions for all samples and evaluates the results using the \itemize{
#'  \item Area Under the Receiver Operating Characteristic (ROC) Curve (AU-ROC)
#'  \item and the Precision-recall Curve (PR)
#' }
#' as metric. Predictions can be supplied either for a single case or as matrix after resampling of the dataset.
#'
#' Prediction results are usually produced with the function \link{plm.predictor}.
#' @param label label object
#' @param pred prediction for each sample by the model, should be a matrix with dimensions \code{length(label) x 1} or \code{length(label) x num.resample}
#' @keywords SIAMCAT eval.result
#' @export
#' @return list containing \itemize{
#'  \item \code{$roc.average} average ROC-curve across repeats or a single ROC-curve on complete dataset;
#'  \item \code{$auc.average} AUC value for the average ROC-curve;
#'  \item \code{$ev.list} list of \code{length(num.folds)}, containing for different decision thresholds the number of false positives, false negatives, true negatives, and true positives;
#'  \item \code{$pr.list} list of \code{length(num.folds)}, containing the positive predictive value (precision) and true positive rate (recall) values used to plot the PR curves;
#' }. If \code{prediction} had more than one column, i.e. if the models has been trained with several repeats, the function will additonally return \itemize{
#'  \item \code{$roc.all} list of roc objects (see \link[pROC]{roc}) for every repeat;
#'  \item \code{$aucspr} vector of AUC values for the PR curves for every repeat;
#'  \item \code{$auc.all} vector of AUC values for the ROC curves for every repeat
#'}
eval.result <- function(label, pred){

  # TODO compare header to label
  ### make sure that label and prediction are in the same order
  #stopifnot(all(names(label) %in% rownames(pred)) && all(rownames(pred) %in% names(label)))
  m    <- match(names(label$label), rownames(pred))
  #cat(m, '\n')
  pred <- pred[m,,drop=FALSE]
  stopifnot(all(names(label$label) == rownames(pred)))

  # ROC curve
  auroc = 0
  if (ncol(pred) > 1) {
    rocc = list(NULL)
    aucs = vector('numeric', ncol(pred))
    for (c in 1:ncol(pred)) {
      rocc[c] = list(roc(response=label$label, predictor=pred[,c], ci=FALSE))
      aucs[c] = rocc[[c]]$auc
    }
    l.vec = rep(label$label, ncol(pred))
  } else {
    l.vec = label$label
  }
  # average data for plotting one mean prediction curve
  summ.stat = 'mean'
  rocsumm = list(roc(response=label$label, predictor=apply(pred, 1, summ.stat),
                 ci=TRUE, of="se", sp=seq(0, 1, 0.05)))
  auroc = list(rocsumm[[1]]$auc)
  # precision recall curve
  pr = list(NULL)
  ev = list(NULL)
  if (ncol(pred) > 1) {
    aucspr = vector('numeric', dim(pred)[2])
    for (c in 1:ncol(pred)) {
      print(c)
      ev[c] = list(eval.classifier(pred[,c], label$label, label))
      pr[c] = list(get.pr(ev[[c]]))
      aucspr[c] = calc.aupr(ev[[c]])
    }
    ev = append(ev,list(eval.classifier(apply(pred, 1, summ.stat), label$label, label)))
  } else {
    ev[1] = list(eval.classifier(as.vector(pred), label$label, label))
    pr[1] = list(get.pr(ev[[1]]))
  }
  if (ncol(pred) > 1) {
    return(list("roc.all" = rocc,
                "auc.all"=aucs,
                "roc.average"=rocsumm,
                "auc.average"=auroc,
                "ev.list"=ev,
                "pr.list"=pr,
                "aucspr"=aucspr))
  } else {
    return(list("roc.average"=rocsumm,
                "auc.average"=auroc,
                "ev.list"=ev,
                "pr.list"=pr))
  }
}
