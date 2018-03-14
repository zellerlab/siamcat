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
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'        for only information about progress and success, \code{2} for normal 
#'        level of information and \code{3} for full debug information, defaults to \code{1}
#' @keywords SIAMCAT model.evaluation.plot
#' @export
#' @return Does not return anything, but produces the model evaluation plot.
model.evaluation.plot <- function(siamcat, fn.plot, verbose=1){
  if(verbose>1) cat("+ starting model.evaluation.plot\n")
  s.time <- proc.time()[3]
  pdf(fn.plot, onefile=TRUE)
  
  if(verbose>2) cat("+ plotting ROC\n")
  plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab='False positive rate', ylab='True positive rate', type='n')
  title(paste('ROC curve for the model', sep=' '))
  abline(a=0, b=1, lty=3)
  if (dim(siamcat@predMatrix)[2] > 1) {
    aucs = vector('numeric', dim(siamcat@predMatrix)[2])
    for (c in 1:dim(siamcat@predMatrix)[2]) {
      roc.c = siamcat@eval_data$roc.all[[c]]
      lines(1-roc.c$specificities, roc.c$sensitivities, col=gray(runif(1,0.2,0.8)))
      aucs[c] = siamcat@eval_data$auc.all[c]
      if(verbose>2) cat('+++ AU-ROC (resampled run ', c, '): ', format(aucs[c], digits=3), '\n', sep='')
    }
    l.vec = rep(siamcat@label@label, dim(siamcat@predMatrix)[2])
  } else {
    l.vec = siamcat@label@label
  }
  roc.summ = siamcat@eval_data$roc.average[[1]]
  lines(1-roc.summ$specificities, roc.summ$sensitivities, col='black', lwd=2)
  auroc = siamcat@eval_data$auc.average[1]
  # plot CI
  x = as.numeric(rownames(roc.summ$ci))
  yl = roc.summ$ci[,1]
  yu = roc.summ$ci[,3]
  polygon(1-c(x, rev(x)), c(yl, rev(yu)), col='#88888844' , border=NA)

  if (dim(siamcat@predMatrix)[2] > 1) {
    if(verbose>1) cat('+ AU-ROC:\n+++ mean-prediction:', format(auroc, digits=3), 
                      '\n+++ averaged       :', format(mean(aucs), digits=3), 
                      '\n+++ sd             :', format(sd(aucs), digits=4), '\n', sep='')
    text(0.7, 0.1, paste('Mean-prediction AUC:', format(auroc, digits=3)))
  } else {
    if(verbose>1) cat('+ AU-ROC:', format(auroc, digits=3), '\n')
    text(0.7, 0.1, paste('AUC:', format(auroc, digits=3)))
  }

  # precision recall curve
  if(verbose>2) cat("+ plotting PRC\n")
  plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab='Recall', ylab='Precision', type='n')
  title(paste('Precision-recall curve for the model', sep=' '))
  abline(h=mean(siamcat@label@label==siamcat@label@positive.lab), lty=3)

  if (dim(siamcat@predMatrix)[2] > 1) {
    aucspr = vector('numeric', dim(siamcat@predMatrix)[2])
    for (c in 1:dim(siamcat@predMatrix)[2]) {
      ev = siamcat@eval_data$ev.list[[c]]
      pr = siamcat@eval_data$pr.list[[c]]
      lines(pr$x, pr$y, col=gray(runif(1,0.2,0.8)))
      aucspr[c] = siamcat@eval_data$aucspr[c]
      if(verbose>2) cat('+++ AU-PRC (resampled run ', c, '): ', format(aucspr[c], digits=3), '\n', sep='')
    }
    ev = siamcat@eval_data$ev.list[[length(siamcat@eval_data$ev.list)]]
  } else {
    ev = siamcat@eval_data$ev.list[[1]]
  }
  pr = evaluate.get.pr(ev,verbose=verbose)
  lines(pr$x, pr$y, col='black', lwd=2)
  aupr = evaluate.calc.aupr(ev,verbose=verbose)
  if (dim(siamcat@predMatrix)[2] > 1) {
    if(verbose>1) cat('+ AU-PRC:\n+++ mean-prediction:', format(aupr, digits=3), 
                      '\n+++ averaged       :', format(mean(aucspr), digits=3), 
                      '\n+++ sd             :', format(sd(aucspr), digits=4), '\n', sep='')
    text(0.7, 0.1, paste('Mean-prediction AU-PRC:', format(aupr, digits=3)))
  } else {
    if(verbose>1) cat('+ AU-PRC:', format(aupr, digits=3), '\n')
    text(0.7, 0.1, paste('AUC:', format(aupr, digits=3)))
  }
  tmp <- dev.off()
  e.time <- proc.time()[3]
  if(verbose>1) cat("+ finished model.evaluation.plot in",e.time-s.time,"s\n")
  if(verbose==1) cat("Plotted evaluation of predictions successfully to:",fn.plot,"\n")
}
