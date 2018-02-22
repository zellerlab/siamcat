###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# RScript flavor
#
# written by Georg Zeller
# with additions by Nicolai Karcher and Konrad Zych
# EMBL Heidelberg 2012-2017
#
# version 0.2.0
# file last updated: 01.06.2017
# GNU GPL 3.0
###

# evaluates the predictive performance of a classifier on a labeled data sets
# returns a list with vectors containing TP, FP, TN, FN for each threshold value on the predictions
# (where TP = true positives, FP = false positives, TN = true negatives, FN = false negatives)
eval.classifier <- function(predictions, test.label, label) {
  stopifnot(dim(test.label) == NULL)
  stopifnot(length(unique(test.label)) == 2)
  stopifnot(all(is.finite(predictions)))
  # calculate thresholds, one between each subsequent pair of sorted prediction values
  # this is ignorant to whether predictions is in matrix or vector format (see below)
  thr      <- predictions
  dim(thr) <- NULL
  thr      <- sort(unique(thr))
  thr      <- rev(c(min(thr)-1, (thr[-1]+thr[-length(thr)])/2, max(thr)+1))
  if (is.null(dim(predictions))) {
    # assuming that a single model was applied to predict the data set
    stopifnot(length(test.label) == length(predictions))
    # actual evaluations per threshold value
    tp = vector('numeric', length(thr))
    fp = vector('numeric', length(thr))
    tn = vector('numeric', length(thr))
    fn = vector('numeric', length(thr))
    for (i in 1:length(thr)) {
      tp[i] = sum(test.label==label@positive.lab & predictions>thr[i])
      fp[i] = sum(test.label==label@negative.lab & predictions>thr[i])
      tn[i] = sum(test.label==label@negative.lab & predictions<thr[i])
      fn[i] = sum(test.label==label@positive.lab & predictions<thr[i])
    }
  } else {
    # assuming that several models were applied to predict the same data and predictions of each model occupy one column
    stopifnot(length(test.label) == nrow(predictions))
    tp = matrix(0, nrow=length(thr), ncol=ncol(predictions))
    fp = matrix(0, nrow=length(thr), ncol=ncol(predictions))
    tn = matrix(0, nrow=length(thr), ncol=ncol(predictions))
    fn = matrix(0, nrow=length(thr), ncol=ncol(predictions))
    for (c in 1:ncol(predictions)) {
      for (r in 1:length(t)) {
        tp[r,c] = sum(test.label==label@positive.lab  & predictions[,c] > thr[r])
        fp[r,c] = sum(test.label==label@negative.lab  & predictions[,c] > thr[r])
        tn[r,c] = sum(test.label==label@negative.lab  & predictions[,c] < thr[r])
        fn[r,c] = sum(test.label==label@positive.lab  & predictions[,c] < thr[r])
      }
    }
  }
  return(list(tp=tp, tn=tn, fp=fp, fn=fn, thresholds=thr))
}

# returns a matrix of x and y values for plotting a receiver operating characteristic curve
# eval is a list produced by eval.classifier
get.roc = function(eval) {
  if (is.null(dim(eval$tp))) {
    stopifnot(!is.null(eval$fp))
    stopifnot(!is.null(eval$tn))
    stopifnot(!is.null(eval$fn))
    fpr = eval$fp / (eval$tn + eval$fp)
    tpr = eval$tp / (eval$tp + eval$fn)
  } else {
    stopifnot(ncol(eval$tp) == ncol(eval$fp))
    stopifnot(ncol(eval$tp) == ncol(eval$tn))
    stopifnot(ncol(eval$tp) == ncol(eval$fn))
    fpr = matrix(NA, nrow=nrow(eval$tp), ncol=ncol(eval$tp))
    tpr = matrix(NA, nrow=nrow(eval$tp), ncol=ncol(eval$tp))
    for (c in 1:ncol(eval$tp)) {
      fpr[,c] = eval$fp[,c] / (eval$tn[,c] + eval$fp[,c])
      tpr[,c] = eval$tp[,c] / (eval$tp[,c] + eval$fn[,c])
    }
  }
  return(list(x=fpr, y=tpr))
}

# plots a receiver operating characteristic curve
plot.roc.curve = function(eval) {
  roc = get.roc(eval)
  plot(roc$x, roc$y,
       xlim=c(0,1), ylim=c(0,1), type='l', # type='o' pch=16, cex=0.4,
       xlab='False positive rate', ylab='True positive rate')
}

# returns a vector of x and y values for plotting a precision-recall curve
get.pr = function(eval) {
  tpr = eval$tp / (eval$tp + eval$fn)
  ppv = eval$tp / (eval$tp + eval$fp)
  # at thresholds where the classifier makes no positive predictions at all,
  # we (somewhat arbitrarily) set its precision to 1
  ppv[is.na(ppv)] = 1.0
  return(list(x=tpr, y=ppv))
}

# plots a precision-recall curve
plot.pr.curve = function(eval) {
  pr = get.pr(eval)
  plot(pr$x, pr$y,
       xlim=c(0,1), ylim=c(0,1), type='l', # type='o' pch=16, cex=0.4,
       xlab='Recall (TPR)', ylab='Precision (PPV)')
}


# calculates the area under a curve using a trapezoid approximation
area.trapez = function(x, y) {
  if (x[1] > x[length(x)]) {
    x = rev(x)
    y = rev(y)
  }
  xd = x[-1] - x[-length(x)]
  ym = 0.5 * (y[-1] + y[-length(y)])
  return(xd %*% ym)
}

# calculates the area under the ROC curve (over the interval [0, max.fpr], if specified)
# get.roc is a function which takes the list returned by eval.classifier and returns a list (x,y) containing TP  and FP values for each
# threshhold set in eval.classifier.
calc.auroc = function(eval, max.fpr=1.0) {
  roc = get.roc(eval)
  idx = roc$x <= max.fpr
  return(area.trapez(roc$x[idx], roc$y[idx]))
}

# calculates the area under the precision-recall curve (over the interval [0, max.tpr], if specified)
calc.aupr = function(eval, max.tpr=1.0) {
  pr = get.pr(eval)
  idx = pr$x <= max.tpr
  return(area.trapez(pr$x[idx], pr$y[idx]))
}

# returns the positions of local minima in a given (ordered) vector
# local minima are defined as having ext adjacent smaller values in both directions
local.min = function(y, ext) {
  m = list()
  for (i in 1:length(y)) {
    env = c((i-ext):(i-1), (i+1):(i+ext))
    env = env[env > 0 & env <= length(y)]
    if (y[i] > max(y[env])) {
      m = c(m,i)
    }
  }
  return(unlist(m))
}

plot.model.eval <- function(fn.plot, pred, eval.data, model.type){
  pdf(fn.plot)
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
    l.vec = rep(label, dim(pred)[2])
  } else {
    l.vec = label
  }
  roc.summ = eval.data$roc.average[[1]]
  lines(1-roc.summ$specificities, roc.summ$sensitivities, col='black', lwd=2)
  auroc = eval.data$auc.average[1]
  # plot CI
  x = as.numeric(rownames(roc.summ$ci))
  yl = roc.summ$ci[,1]
  yu = roc.summ$ci[,3]
  polygon(1-c(x, rev(x)), c(yl, rev(yu)), col='#88888844', border=NA)

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
  abline(h=mean(label==PL), lty=3)

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

  tmp = dev.off()
}
