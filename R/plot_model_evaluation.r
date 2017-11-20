###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# R package flavor
#
# written by Georg Zeller
# with additions by Nicolai Karcher and Konrad Zych
# EMBL Heidelberg 2012-2017
#
# version 0.2.0
# file last updated: 26.06.2017
# GNU GPL-3.0
###

#' @title Model Evaluation Plot
#' @description Produces two plots for model evaluation. The first plot shows the Receiver Operating Receiver (ROC)-curves, the other the Precision-recall (PR)-curves for the different CV repetitions.
#' @param label a label object
#' @param pred matrix containing the model predictions for every CV repetition
#' @param eval.data list of evaluation results produced by \link{eval.results}
#' @param model.type type of the model, defaults to \code{"lasso"}
#' @keywords SIAMCAT evaluation.model.plot
#' @export
#' @return Does not return anything, but produces the model evaluation plot.
#' @examples
#' pdf(filename.pdf)
#' evaluation.model.plot(label, pred, eval.data, model.type)
#' dev.off()
evaluation.model.plot <- function(label, pred, eval.data, model.type='lasso'){
  plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab='False positive rate', ylab='True positive rate', type='n')
  title(paste('ROC curve for', model.type, 'model', sep=' '))
  abline(a=0, b=1, lty=3)
  if (dim(pred)[2] > 1) {
    aucs = vector('numeric', dim(pred)[2])
    for (c in 1:dim(pred)[2]) {
      roc.c = eval.data$roc.all[[c]]
      lines(1-roc.c$specificities, roc.c$sensitivities, col=gray(runif(1,0.2,0.8)))
      aucs[c] = eval.data$auc.all[c]
      cat('AU-ROC (resampled run ', c, '): ', format(aucs[c], digits=3), '\n', sep='')
    }
    l.vec = rep(label$label, dim(pred)[2])
  } else {
    l.vec = label$label
  }
  roc.summ = eval.data$roc.average[[1]]
  lines(1-roc.summ$specificities, roc.summ$sensitivities, col='black', lwd=2)
  auroc = eval.data$auc.average[1]
  # plot CI
  x = as.numeric(rownames(roc.summ$ci))
  yl = roc.summ$ci[,1]
  yu = roc.summ$ci[,3]
  polygon(1-c(x, rev(x)), c(yl, rev(yu)), col='#88888844' , border=NA)

  if (dim(pred)[2] > 1) {
    cat('Mean-pred. AU-ROC:', format(auroc, digits=3), '\n')
    cat('Averaged AU-ROC: ', format(mean(aucs), digits=3), ' (sd=', format(sd(aucs), digits=4), ')\n', sep='')
    text(0.7, 0.1, paste('Mean-prediction AUC:', format(auroc, digits=3)))
  } else {
    cat('AU-ROC:', format(auroc, digits=3), '\n')
    text(0.7, 0.1, paste('AUC:', format(auroc, digits=3)))
  }

  # precision recall curve
  plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab='Recall', ylab='Precision', type='n')
  title(paste('Precision-recall curve for', model.type, 'model', sep=' '))
  abline(h=mean(label$label==label$positive.lab), lty=3)

  if (dim(pred)[2] > 1) {
    aucspr = vector('numeric', dim(pred)[2])
    for (c in 1:dim(pred)[2]) {
      ev = eval.data$ev.list[[c]]
      pr = eval.data$pr.list[[c]]
      lines(pr$x, pr$y, col=gray(runif(1,0.2,0.8)))
      aucspr[c] = eval.data$aucspr[c]
      cat('AU-PRC (resampled run ', c, '): ', format(aucspr[c], digits=3), '\n', sep='')
    }
    ev = eval.data$ev.list[[length(eval.data$ev.list)]]
  } else {
    ev = eval.data$ev.list[[1]]
  }
  pr = get.pr(ev)
  lines(pr$x, pr$y, col='black', lwd=2)
  aupr = calc.aupr(ev)
  if (dim(pred)[2] > 1) {
    cat('Mean-pred. AU-PRC:', format(aupr, digits=3), '\n')
    cat('Averaged AU-PRC: ', format(mean(aucs), digits=3), ' (sd=', format(sd(aucs), digits=4), ')\n', sep='')
    text(0.7, 0.1, paste('Mean-prediction AUC:', format(aupr, digits=3)))
  } else {
    cat('AU-PRC:', format(aupr, digits=3), '\n')
    text(0.7, 0.1, paste('AUC:', format(aupr, digits=3)))
  }
}
