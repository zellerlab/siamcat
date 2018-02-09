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
#' @param siamcat object of class \link{siamcat-class}
#' @param fn.plot string, filename for the pdf-plot
#' @keywords SIAMCAT evaluation.model.plot
#' @export
#' @return Does not return anything, but produces the model evaluation plot.
evaluation.model.plot <- function(siamcat, fn.plot){

  pdf(fn.plot, onefile=TRUE)

  plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab='False positive rate', ylab='True positive rate', type='n')
  title(paste('ROC curve for the model', sep=' '))
  abline(a=0, b=1, lty=3)
  if (dim(siamcat@predMatrix)[2] > 1) {
    aucs = vector('numeric', dim(siamcat@predMatrix)[2])
    for (c in 1:dim(siamcat@predMatrix)[2]) {
      roc.c = siamcat@evalData$roc.all[[c]]
      lines(1-roc.c$specificities, roc.c$sensitivities, col=gray(runif(1,0.2,0.8)))
      aucs[c] = siamcat@evalData$auc.all[c]
      cat('AU-ROC (resampled run ', c, '): ', format(aucs[c], digits=3), '\n', sep='')
    }
    l.vec = rep(siamcat@label@label, dim(siamcat@predMatrix)[2])
  } else {
    l.vec = siamcat@label@label
  }
  roc.summ = siamcat@evalData$roc.average[[1]]
  lines(1-roc.summ$specificities, roc.summ$sensitivities, col='black', lwd=2)
  auroc = siamcat@evalData$auc.average[1]
  # plot CI
  x = as.numeric(rownames(roc.summ$ci))
  yl = roc.summ$ci[,1]
  yu = roc.summ$ci[,3]
  polygon(1-c(x, rev(x)), c(yl, rev(yu)), col='#88888844' , border=NA)

  if (dim(siamcat@predMatrix)[2] > 1) {
    cat('Mean-pred. AU-ROC:', format(auroc, digits=3), '\n')
    cat('Averaged AU-ROC: ', format(mean(aucs), digits=3), ' (sd=', format(sd(aucs), digits=4), ')\n', sep='')
    text(0.7, 0.1, paste('Mean-prediction AUC:', format(auroc, digits=3)))
  } else {
    cat('AU-ROC:', format(auroc, digits=3), '\n')
    text(0.7, 0.1, paste('AUC:', format(auroc, digits=3)))
  }

  # precision recall curve
  plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab='Recall', ylab='Precision', type='n')
  title(paste('Precision-recall curve for the model', sep=' '))
  abline(h=mean(siamcat@label@label==siamcat@label@positive.lab), lty=3)

  if (dim(siamcat@predMatrix)[2] > 1) {
    aucspr = vector('numeric', dim(siamcat@predMatrix)[2])
    for (c in 1:dim(siamcat@predMatrix)[2]) {
      ev = siamcat@evalData$ev.list[[c]]
      pr = siamcat@evalData$pr.list[[c]]
      lines(pr$x, pr$y, col=gray(runif(1,0.2,0.8)))
      aucspr[c] = siamcat@evalData$aucspr[c]
      cat('AU-PRC (resampled run ', c, '): ', format(aucspr[c], digits=3), '\n', sep='')
    }
    ev = siamcat@evalData$ev.list[[length(siamcat@evalData$ev.list)]]
  } else {
    ev = siamcat@evalData$ev.list[[1]]
  }
  pr = get.pr(ev)
  lines(pr$x, pr$y, col='black', lwd=2)
  aupr = calc.aupr(ev)
  if (dim(siamcat@predMatrix)[2] > 1) {
    cat('Mean-pred. AU-PRC:', format(aupr, digits=3), '\n')
    cat('Averaged AU-PRC: ', format(mean(aucs), digits=3), ' (sd=', format(sd(aucs), digits=4), ')\n', sep='')
    text(0.7, 0.1, paste('Mean-prediction AUC:', format(aupr, digits=3)))
  } else {
    cat('AU-PRC:', format(aupr, digits=3), '\n')
    text(0.7, 0.1, paste('AUC:', format(aupr, digits=3)))
  }
  tmp <- dev.off()
}
