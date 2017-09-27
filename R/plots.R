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

#' @title model interpretation plot
#' @export
interpretor.model.plot <- function(feat, label, meta, model, pred, color.scheme, consens.thres){
  num.models   <- ncol(model$W)
  ### some color preprocessing
  if (color.scheme == 'matlab') {
    color.scheme = matlab.like(100)
  } else {
    # TODO check for valid param!
    color.scheme = rev(colorRampPalette(brewer.pal(11, color.scheme))(100))
  }
  col.p = color.scheme[length(color.scheme)-4]
  col.n = color.scheme[1+4]


  ### preprocess models: discard hyperpar terms, (optionally normalize) and handle ones that are completely 0
  # keep.idx = grep('Bias', rownames(model$W), invert=TRUE)
  keep.idx = grep('optimal hyperpar', rownames(model$W), invert=TRUE)
  model$W = model$W[keep.idx,]
  sum.w = colSums(abs(model$W))
  sum.w[sum.w == 0] = 1
  if (norm.models) {
    for (m in 1:dim(model$W)[2]) {
      model$W[,m] = model$W[,m] / sum.w[m]
    }
    sum.w = rep(1, dim(model$W)[2])
  }
  sel.idx = which(rowSums(model$W != 0) / dim(model$W)[2] >= consens.thres)
  sel.W = t(model$W[sel.idx,])

  # sort by mean relative model weight
  rel.model.weights = t(model$W) / rowSums(abs(t(model$W)))
  median.sorted.models <- sort(apply(rel.model.weights[,sel.idx], 2, median), decreasing=TRUE, index.return=TRUE)
  # Restrict amount of features to be plotted (Print 25 features with most positive and most negative feature weights, respectively)
  if (length(sel.idx) > 50){
    print("Restricting amount of features to be plotted to 50")
    median.sorted.models$x <- c(head(median.sorted.models$x,n = 25), tail(median.sorted.models$x, n = 25))
    median.sorted.models$ix <- c(head(median.sorted.models$ix,n = 25), tail(median.sorted.models$ix, n = 25))
  }
  sel.idx = sel.idx[median.sorted.models$ix]

  rel.model.weights = rel.model.weights[,sel.idx]
  #print(sel.idx)
  #print(rel.model.weights)
  sel.W = t(model$W[sel.idx,])
  num.sel.f = length(sel.idx)

  cat('Generating plot for a model with', num.sel.f, 'selected features\n')

  ### aggregate predictions of several models if more than one is given
  agg.pred = matrix(data=NA, nrow=num.models, ncol=length(label$label))
  for (s in 1:length(label$label)) {
    subj = names(label$label)[s]
    idx = which(rownames(pred) == subj)
    agg.pred[,s] = pred[idx,]
  }
  mean.agg.pred = agg.pred
  if (dim(agg.pred)[1] > 1) {
    mean.agg.pred = colMeans(agg.pred)
  }
  #cat(mean.agg.pred, '\n')

  ### idx to sort samples according to their class membership and prediction score
  srt.idx = sort(label$label+mean.agg.pred, index.return=TRUE)$ix


  ### start plotting model properties


  ### plot layout
  sel.f.cex = max(0.3, 0.8 - 0.01*num.sel.f)
  lmat = rbind(c(8, 1, 5, 2), c(3, 4,0, 6), c(0, 7, 0, 0))
  h_t = 0.10
  h_m = ifelse(is.null(meta), 0.8, max(0.5, 0.7-0.01*dim(meta)[2]))
  h_b = ifelse(is.null(meta), 0.1, 0.1+0.02*dim(meta)[2])
  cat('Layout height values: ', h_t, ', ', h_m, ', ', h_b, '\n', sep='')
  layout(lmat, widths=c(0.14, 0.58, 0.1, 0.14), heights=c(h_t, h_m, h_b))
  par(oma=c(3, 4, 3, 4))

  ### header row

  # field 1 will be plotted later together with feature heatmap

  # field 2: header for feature heatmap
  par(mar=c(0, 4.1, 3.1, 5.1))
  hm.label = label$label[srt.idx]
  plot(NULL, type='n', xlim=c(0,length(hm.label)), xaxs='i', xaxt='n',
       ylim=c(-0.5,0.5), yaxs='i', yaxt='n', xlab='', ylab='', bty='n')
  ul = unique(hm.label)
  for (l in 1:length(ul)) {
    idx = which(ul[l] == hm.label)
    lines(c(idx[1]-0.8, idx[length(idx)]-0.2), c(0, 0))
    lines(c(idx[1]-0.8, idx[1]-0.8), c(-0.2, 0))
    lines(c(idx[length(idx)]-0.2, idx[length(idx)]-0.2), c(-0.2, 0))
    h = (idx[1] + idx[length(idx)]) / 2
    t = gsub('_', ' ', names(label$info$class.descr)[label$info$class.descr==ul[l]])
    t = paste(t, ' (n=', length(idx), ')', sep='')
    mtext(t, side=3, line=-0.5, at=h, cex=0.7, adj=0.5)
  }
  mtext('Metagenomic features', side=3, line=2, at=length(hm.label)/2, cex=1, adj=0.5)

  # field 3: model header
  par(mar=c(0, 6.1, 3.1, 1.1))
  plot(NULL, type='n', xlim=c(-0.1,0.1), xaxt='n', xlab='',
       ylim=c(-0.1,0.1), yaxt='n', ylab='', bty='n')
  mtext('Linear model', side=3, line=2, at=0.04, cex=1, adj=0.5)
  mtext(paste('(|W| = ', num.sel.f, ')', sep=''), side=3, line=1, at=0.04, cex=0.7, adj=0.5)


  ### first data display row


  med = apply(rel.model.weights, 2, median)
  low.qt = apply(rel.model.weights, 2, quantile)[2,]
  upp.qt = apply(rel.model.weights, 2, quantile)[4,]
  # field 4: barplot of effect size associated with each feature
  par(mar=c(0.1, 1.1, 0, 1.1))
  mi = min(-med-(abs(low.qt-upp.qt)))
  mx = max(-med+(abs(low.qt-upp.qt)))
  barplot(-med, horiz = TRUE, width=1, space=0, yaxs='i',
          xlim=c(-max(abs(mi),abs(mx)), max(abs(mi), max(mx))), ylim=c(0, num.sel.f), xlab='', ylab='', yaxt='n')
  # to change the background color of the plot, invoking barplot twice seems easiest
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col='gray90', border=NA)
  plot.coords <- barplot(-med, col='gray30', border='white', horiz=TRUE, width=1, space=0, yaxs='i', xlab='', ylab='', yaxt='n', add=TRUE)
  draw.error.bar(plot.coords = plot.coords, -low.qt, -upp.qt)
  box(lwd=1)
  mtext('median relative feat. weight', side=1, line=2, at=0, cex=0.7, adj=0.5)
  # robustness indicated as percentage of models including a given feature (to the right of the barplot)
  for (f in 1:num.sel.f) {
    t = paste(format(100*rowSums(model$W[sel.idx[f],] != 0) / num.models, digits=1, scientific=FALSE), '%', sep='')
    mtext(t, side=4, line=2.5, at=(f-0.5), cex=sel.f.cex, las=2, adj=1)
    if (f == floor(num.sel.f/2)){
      mtext(gsub('_', ' ', names(label$info$class.descr)[label$info$class.descr==ul[1]]),
            side = 2,
            at = f,
            line = -2
      )
      mtext(gsub('_', ' ', names(label$info$class.descr)[label$info$class.descr==ul[2]]),
            side = 4,
            at = f,
            line = -2
      )
    }
  }
  mtext('effect size', side=3, line=1, at=(mx/2), cex=0.7, adj=1)
  mtext('robustness', side=3, line=1, at=mx, cex=0.7, adj=0)

  cat('  finished plotting feature weights.\n')


  # field 5: feature heatmap with feature names to the right

  par(mar=c(0.1, 4.1, 0, 5.1))
  if (heatmap.type == 'zscore'){
    # data is  transposed and transformed to feature z-scores for display
    img.data = t(feat[sel.idx, srt.idx])
    m = apply(img.data, 2, mean)
    s = apply(img.data, 2, sd)
    for (c in 1:dim(img.data)[2]) {
      img.data[,c] = (img.data[,c] - m[c]) / s[c]
    }
    zlim = z.score.lim
  } else if (heatmap.type == 'fc') {
    # extract only those features which are specified by sel.idx. Necessary since
    # the indices specified in sel.idx were generated on feat (and not on origin.feat)
    m = match(rownames(feat)[sel.idx], rownames(origin.feat))
    origin.feat.sel = origin.feat[m,]

    feat.ct.median = apply(origin.feat.sel[,label$n.idx], 1, median)
    img.data = log10(origin.feat.sel + detect.lim) - log10(feat.ct.median + detect.lim)
    # reorder columns
    img.data = img.data[,srt.idx]
    if (any(is.na(m))) {
      # this can be the case if meta-variables have been added as predictors
      # (then there's no corresponding original feature)
      idx = which(is.na(m))
      img.data[idx,] = feat[sel.idx[idx],]
      rownames(img.data)[idx] = rownames(feat)[sel.idx[idx]]
    }
    # transpose for heatmap plot
    img.data = t(img.data)
    img.data = img.data
    zlim = fc.lim
  } else {
    stop('unknown heatmap.type: ', heatmap.type)
  }
  # truncate extreme values for heatmap visualization
  img.data[img.data < zlim[1]] = zlim[1]
  img.data[img.data > zlim[2]] = zlim[2]
  image(img.data, zlim=zlim, col=color.scheme, xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
  for (f in 1:num.sel.f) {
    mtext(colnames(sel.W)[f], side=4, line=1, at=(f-1)/(num.sel.f-1), cex=sel.f.cex, las=2,
          col=ifelse(med[f]>0, col.n, col.p))
  }
  box(lwd=1)
  cat('  finished plotting feature heatmap.\n')

  # additonally add to header row a corresponding color bar
  par(mar=c(3.1, 1.1, 1.1, 1.1))
  # TODO would be nice to overplot a histogram
  #h = as.vector(img.data)
  #h = h[h > z.score.lim[1] & h < z.score.lim[2]]
  #h = hist(h, 100, plot=FALSE)$counts
  # end TODO
  barplot(as.matrix(rep(1,100)), col = color.scheme, horiz=TRUE, border=0, ylab='', axes=FALSE)
  if (heatmap.type == 'fc') {
    key.ticks = seq(round(min(img.data), digits = 1), round(max(img.data), digits = 1), length.out=7)
    key.label = 'Feature fold change over controls'
  } else if (heatmap.type == 'zscore') {
    key.ticks = seq(z.score.lim[1], z.score.lim[2], length.out=7)
    key.label = 'Feature z-score'
  }
  axis(side=1, at=seq(0, 100, length.out=7), labels=key.ticks)
  mtext(key.label, side=3, line=0.5, at=50, cex=0.7, adj=0.5)



  # field 6: boxplot displaying the poportion of weight per model that is actually shown
  par(mar=c(0.1, 6.1, 0, 1.1))
  boxplot(rowSums(abs(sel.W)) / sum.w, ylim=c(0,1))
  mtext('proportion of', side=1, line=1, at=1, adj=0.5, cex=0.7)
  mtext('weight shown', side=1, line=2, at=1, adj=0.5, cex=0.7)
  cat('  finished plotting proportion of model weight shown.\n')


  # empty field (left)

  # field 7 (middle): heatmap showing predictions and metadata (if given)
  par(mar=c(1.1, 4.1, 0.3, 5.1))
  img.data = as.matrix(mean.agg.pred[srt.idx])
  colnames(img.data) = 'Classification result'
  if (!is.null(meta)) {
    img.data = cbind(meta[srt.idx, dim(meta)[2]:1], img.data)
  }
  ### transform any categorial column into a numeric one
  for (m in 1:dim(img.data)[2]) {
    img.data[,m] = (img.data[,m] - min(img.data[,m], na.rm=TRUE))
    if (max(img.data[,m], na.rm=TRUE) != 0) {
      img.data[,m] = img.data[,m] / max(img.data[,m], na.rm=TRUE)
    }
  }

  grays = rev(gray(seq(0, 1, length.out=length(color.scheme))))
  image(as.matrix(img.data), col=grays, xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
  box(lwd=1)

  meta.cex = max(0.3, 0.7 - 0.01*dim(img.data)[2])
  for (m in 1:dim(img.data)[2]) {
    t = colnames(img.data)[m]
    t = gsub('\\.|_', ' ', t)
    mtext(t, side=4, line=1, at=(m-1)/(dim(img.data)[2]-1), cex=meta.cex, las=2)
  }
  # mark missing values
  for (m in 1:dim(img.data)[2]) {
    idx = which(is.na(img.data[,m]))
    for (i in idx) {
      x = (i-1) / (dim(img.data)[1]-1)
      y = (m-1) / (dim(img.data)[2]-1)
      text(x, y, 'NA', col='red', cex=0.4)
    }
  }
  cat('  finished plotting classification result and additional metadata.\n')

  # empty field (right)

  # field 8  : header for  feature weight barplot
  par(mar=c(0, 1.1, 3.1, 1.1))
  plot(NULL, type='n', xlim=c(-0.1,0.1), xaxt='n', xlab='',
       ylim=c(-0.1,0.1), yaxt='n', ylab='', bty='n')
  mtext('Feature Weights', side=3, line=2, at=0.04, cex=1, adj=0.5)


}

#' @title plor model evaluation plot
#' @export
evaluation.model.plot <- function(fn.plot, label, pred, eval.data, model.type){
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

# ### xv is vector containing values to draw a barplot, z and y determine upper and lower boundary of barplot, respectively.
draw.error.bar <- function(plot.coords, z, y){
  g <- (max(plot.coords)-min(plot.coords))/(3*length(plot.coords))
  for (i in 1:length(plot.coords)) {
    lines(c(z[i],y[i]),c(plot.coords[i], plot.coords[i]))
    lines(c(z[i],z[i]),c(plot.coords[i]+g, plot.coords[i]-g))
    lines(c(y[i],y[i]),c(plot.coords[i]+g, plot.coords[i]-g))
  }
}
