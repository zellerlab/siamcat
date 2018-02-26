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

#' @title Model Interpretation Plot
#' @description Produces a plot for model interpretation, displaying feature weights, robustness of feature weights, and features scores across patients.
#' @param siamcat object of class \link{siamcat-class}
#' @param fn.plot string, filename for the pdf-plot
#' @param color.scheme color scheme for the heatmap, defaults to =\code{"BrBG"}
#' @param consens.thres minimal ratio of models incorporating a feature in order to include it into the heatmap, defaults to \code{0.5}
#'        Note that for \code{"randomForest"} models, this cutoff specifies the minimum median Gini coefficient for a feature to be included and should therefore be much lower, e.g. \code{0.01}
#' @param heatmap.type type of the heatmap, can be either \code{"fc"} or \code{"zscore"}, defaults to \code{"zscore"}
#' @param norm.models boolean, should the feature weights be normalized across models?, defaults to \code{FALSE}
#' @param limits vector, cutoff for extreme values in the heatmap, defaults to \code{c(-3, 3)}
#' @param detect.lim float, pseudocount to be added before log10-transformation of features, defaults to \code{1e-08}
#' @param max.show integer, maximum number of features to be shown in the model interpretation plot, defaults to 50
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'        for only information about progress and success, \code{2} for normal 
#'        level of information and \code{3} for full debug information, defaults to \code{1}
#' @keywords SIAMCAT interpretor.model.plot
#' @return Does not return anything, but produces the model interpretion plot.
#' @export
#' @details Produces a plot consisting of \itemize{
#'  \item a barplot showing the feature weights and their robustness (i.e. in what proportion of models have they been incorporated)
#'  \item a heatmap showing the z-scores of the metagenomic features across patients
#'  \item another heatmap displaying the metadata categories (if applicable)
#'  \item a boxplot displaying the poportion of weight per model that is actually shown
#' for the features that are incorporated into more than \code{consens.thres} percent of the models.
#'}
interpretor.model.plot <- function(siamcat, fn.plot,
  color.scheme='BrBG', consens.thres=0.5,
  heatmap.type=c('zscore',"fc"), norm.models=FALSE,
  limits=c(-3,3), detect.lim=1e-08, max.show=50, verbose=1){
  if(verbose>1) cat("+ starting evaluation.model.plot\n")
  s.time <- proc.time()[3]
  # ############################################################################
  ### some color pre-processing
  if(verbose>2) cat("+++ preprocessing color scheme\n")
  if (!color.scheme %in% row.names(brewer.pal.info)){
    warning("Not a valid RColorBrewer palette name, defaulting to BrBG...\n
    See brewer.pal.info for more information about RColorBrewer palettes...")
    color.scheme <- 'BrBG'
  }
  color.scheme <- rev(colorRampPalette(brewer.pal(brewer.pal.info[color.scheme,'maxcolors'], color.scheme))(100))

  # ############################################################################
  # get model type from model
  if(verbose>2) cat("+++ retrieving model type\n")
  W.mat       <- get.weights.matrix(siamcat@modelList@models, verbose=verbose)
  all.weights <- W.mat[union(row.names(siamcat@phyloseq@otu_table), grep('META', row.names(W.mat), value = TRUE)),] # remove possible intercept parameters, but keep possible meta data included in the model
  rel.weights <- apply(all.weights, 2, function(x){x/sum(abs(x))})
  # ############################################################################
  ### preprocess models
  if(verbose>2) cat("+++ preprocessing models\n")
  sel.idx <- select.features(weights=all.weights,
                             model.type=siamcat@modelList@model.type,
                             consens.thres=consens.thres,
                             label=siamcat@label,
                             norm.models=norm.models,
                             max.show=max.show, verbose=verbose)
  num.sel.f <- length(sel.idx)
  # ############################################################################
  # aggreate predictions and sort patients by score
  ### aggregate predictions of several models if more than one is given
  mean.agg.pred <- rowMeans(siamcat@predMatrix)
  ### idx to sort samples according to their class membership and prediction score
  srt.idx <- sort(siamcat@label@label + mean.agg.pred, index.return=TRUE)$ix
  # ############################################################################
  # prepare heatmap
  if(verbose>2) cat("+++ preparing heatmap\n")
  if (heatmap.type == 'zscore'){
    feat <- matrix(siamcat@phyloseq@otu_table,nrow=nrow(siamcat@phyloseq@otu_table), ncol=ncol(siamcat@phyloseq@otu_table),
                 dimnames = list(rownames(siamcat@phyloseq@otu_table), colnames(siamcat@phyloseq@otu_table)))
    img.data <- prepare.heatmap.zscore(heatmap.data=feat[sel.idx, srt.idx],
                                       limits=limits)
  } else if (heatmap.type == 'fc') {
    feat <- matrix(siamcat@orig_feat,nrow=nrow(siamcat@orig_feat), ncol=ncol(siamcat@orig_feat),
                 dimnames = list(rownames(siamcat@orig_feat), colnames(siamcat@orig_feat)))
    if (is.null(detect.lim)){
      warning("WARNING: Pseudo-count before log-transformation not supplied! Estimating it as 5% percentile...\n")
      detect.lim <- quantile(feat[feat!=0], 0.05)
    }
    img.data <- prepare.heatmap.fc(heatmap.data=feat[, srt.idx],
                                   sel.feat=names(sel.idx),
                                   limits=limits, meta=siamcat@phyloseq@sam_data,
                                   label=siamcat@label, detect.lim=detect.lim)
  } else {
    stop('! unknown heatmap.type: ', heatmap.type)
  }

  # ############################################################################
  ### start plotting model properties
  if(verbose>2) cat("+++ plotting model properties\n")
  pdf(fn.plot, paper='special', height=8.27, width=11.69, onefile=TRUE)

  ### plot layout
  sel.f.cex <- max(0.3, 0.8 - 0.01*num.sel.f)
  lmat <- rbind(c(1, 2, 3, 4),
               c(5, 6, 0, 7),
               c(0, 8, 0, 0))
  h_t <- 0.10
  h_m <- ifelse(is.null(siamcat@phyloseq@sam_data), 0.8, max(0.5, 0.7-0.01*dim(siamcat@phyloseq@sam_data)[2]))
  h_b <- 1 - h_t - h_m # ifelse(is.null(meta), 0.1, 0.1+0.02*dim(meta)[2])
  cat('Layout height values: ', h_t, ', ', h_m, ', ', h_b, '\n', sep='')
  layout(lmat, widths=c(0.14, 0.58, 0.1, 0.14), heights=c(h_t, h_m, h_b))
  par(oma=c(3, 4, 3, 4))

  ### header row
  # ############################################################################
  # Title of Feature Weights plot
  if(verbose>2) cat("+++ plotting titles\n")
  par(mar=c(0, 1.1, 3.1, 1.1))
  plot(NULL, type='n', xlim=c(-0.1,0.1), xaxt='n', xlab='',
       ylim=c(-0.1,0.1), yaxt='n', ylab='', bty='n')
  mtext('Feature Weights', side=3, line=2, at=0.04, cex=1, adj=0.5)

  # ############################################################################
  # Title of heatmap and brackets for classes
  par(mar=c(0, 4.1, 3.1, 5.1))
  hm.label <- siamcat@label@label[srt.idx]
  plot(NULL, type='n', xlim=c(0,length(hm.label)), xaxs='i', xaxt='n',
       ylim=c(-0.5,0.5), yaxs='i', yaxt='n', xlab='', ylab='', bty='n')
  ul <- unique(hm.label)
  for (l in 1:length(ul)) {
    idx <- which(ul[l] == hm.label)
    lines(c(idx[1]-0.8, idx[length(idx)]-0.2), c(0, 0))
    lines(c(idx[1]-0.8, idx[1]-0.8), c(-0.2, 0))
    lines(c(idx[length(idx)]-0.2, idx[length(idx)]-0.2), c(-0.2, 0))
    h <- (idx[1] + idx[length(idx)]) / 2
    t <- gsub('_', ' ', names(siamcat@label@info$class.descr)[siamcat@label@info$class.descr==ul[l]])
    t <- paste(t, ' (n=', length(idx), ')', sep='')
    mtext(t, side=3, line=-0.5, at=h, cex=0.7, adj=0.5)
  }
  mtext('Metagenomic Features', side=3, line=2, at=length(hm.label)/2, cex=1, adj=0.5)

  # ############################################################################
  # Heatmap legend
  # LOW PRIORITY #
  # TODO would be nice to overplot a histogram
  #h = as.vector(img.data)
  #h = h[h > z.score.lim[1] & h < z.score.lim[2]]
  #h = hist(h, 100, plot=FALSE)$counts
  # end TODO
  # LOW PRIORITY #
  if(verbose>2) cat("+++ plotting legend\n")
  par(mar=c(3.1, 1.1, 1.1, 1.1))
  barplot(as.matrix(rep(1,100)), col = color.scheme,
          horiz=TRUE, border=0, ylab='', axes=FALSE)
  if (heatmap.type == 'fc') {
    key.ticks <- seq(round(min(img.data, na.rm=TRUE), digits = 1),
                    round(max(img.data, na.rm=TRUE), digits = 1), length.out=7)
    key.label <- 'Feature fold change over controls'
  } else if (heatmap.type == 'zscore') {
    key.ticks <- seq(limits[1], limits[2], length.out=7)
    key.label <- 'Feature z-score'
  }
  axis(side=1, at=seq(0, 100, length.out=7), labels=key.ticks)
  mtext(key.label, side=3, line=0.5, at=50, cex=0.7, adj=0.5)

  # ############################################################################
  # Model header (model sensitive)
  par(mar=c(0, 6.1, 3.1, 1.1))
  plot(NULL, type='n', xlim=c(-0.1,0.1), xaxt='n', xlab='',
       ylim=c(-0.1,0.1), yaxt='n', ylab='', bty='n')
  mtext(paste0(siamcat@modelList@model.type, ' model'), side=3, line=2, at=0.04, cex=.7, adj=0.5)
  mtext(paste('(|W| = ', num.sel.f, ')', sep=''), side=3, line=1, at=0.04, cex=0.7, adj=0.5)

  # ############################################################################
  # Feature weights ( model sensitive)
  if(verbose>2) cat("+++ plotting feature weights \n")
  plot.feature.weights(rel.weights=rel.weights, sel.idx=sel.idx,
                       mod.type=siamcat@modelList@model.type, label=siamcat@label)

  # ############################################################################
  # Heatmap
  if(verbose>2) cat("+++ plotting heatmap \n")
  if (siamcat@modelList@model.type!='RandomForest'){
    plot.heatmap(image.data=img.data,
                 limits=limits,
                 color.scheme=color.scheme,
                 effect.size=apply(rel.weights[sel.idx,], 1, median),verbose=verbose)
  } else {
    auroc.effect <-  sapply(colnames(img.data),
        FUN=function(f){roc(predictor=img.data[,f],
                        response=label@label, direction='<')$auc})
    bin.auroc.effect <- ifelse(auroc.effect >= 0.5, 1, 0)
    plot.heatmap(image.data=img.data,
                 limits=limits,
                 color.scheme=color.scheme,
                 effect.size=NULL, verbose=verbose)
  }

  # ############################################################################
  # Proportion of weights shown
  # if (model.type != 'RandomForest'){
  if(verbose>2) cat("+++ plotting proportion of weights shown \n")
    plot.proportion.of.weights(selected.weights=all.weights[sel.idx,],
                               all.weights=all.weights,verbose=verbose)
  # } else {
    # cat(dim(all.weights[sel.idx,]), '\n')
    # cat(dim(all.weights), '\n')
    # plot.percentage.of.features(selected.weights=all.weights[sel.idx,],
                               # all.weights=all.weights)
  # }

  # ############################################################################
  # Metadata and prediction
  if(verbose>2) cat("+++ plotting metadata and predictions\n")
  plot.pred.and.meta(prediction=mean.agg.pred[srt.idx],
                     label=siamcat@label,
                     meta=siamcat@phyloseq@sam_data[srt.idx,],verbose=verbose)

  tmp <- dev.off()
  e.time <- proc.time()[3]
  if(verbose>1) cat("+ finished check.associations in",e.time-s.time,"s\n")
  if(verbose==1) cat("Plotted associations between features and label successfully to:",fn.plot,"\n")
}

plot.feature.weights <- function(rel.weights, sel.idx, mod.type, label, verbose=0){
  if(verbose>2) cat("+ plot.feature.weights\n")
  if (mod.type != 'RandomForest'){
    relative.weights <- rel.weights[sel.idx,]
  } else {
    relative.weights <- -rel.weights[sel.idx,]
  }
  med = apply(relative.weights, 1, median)
  low.qt = apply(relative.weights, 1, quantile)[2,] # 25%
  upp.qt = apply(relative.weights, 1, quantile)[4,] # 75%

  if (mod.type != 'RandomForest'){
    par(mar=c(0.1, 1.1, 0, 1.1))
    mi = min(-med-(abs(low.qt-upp.qt)))
    mx = max(-med+(abs(low.qt-upp.qt)))

    barplot(-med, horiz = TRUE, width=1, space=0, yaxs='i', col='gray30',
            xlim=c(-max(abs(mi),abs(mx)), max(abs(mi), max(mx))),
            ylim=c(0, dim(relative.weights)[1]), xlab='', ylab='', yaxt='n')
    x<-par("usr")
    rect(x[1],x[3],x[2],x[4],col="lightgray")
    # grid lines
    grid(NULL, NA, lty=2, col='gray99')
    # plot the barplot again due to the idiotic way of how to change the background color of a plot in r
    barplot(-med, horiz = TRUE, width=1, space=0, yaxs='i', col='gray30',
            xlim=c(-max(abs(mi),abs(mx)), max(abs(mi), max(mx))),
            ylim=c(0, dim(relative.weights)[1]), xlab='', ylab='', yaxt='n', add=TRUE)

    # error bars
    arrows(y0=c(1:length(med))-.5, x0=-upp.qt, x1=-low.qt, angle=90, code=3, length=.04)
    mtext('median relative feat. weight', side=1, line=2, at=0, cex=0.7, adj=0.5)
    # robustness indicated as percentage of models including a given feature (to the right of the barplot)
    for (f in 1:length(sel.idx)) {
      t = paste(format(100*sum(rel.weights[sel.idx[f],] != 0) / dim(rel.weights)[2],
                digits=1, scientific=FALSE), '%', sep='')
      mtext(t, side=4, line=2.5, at=(f-0.5),
            cex=max(0.3, 0.8 - 0.01*length(sel.idx)), las=2, adj=1)
    }

    # label cancer/healthy
    mtext(gsub('_', ' ', label@n.lab), side = 2,
          at = floor(length(sel.idx)/2), line = -2)
    mtext(gsub('_', ' ', label@p.lab), side = 4,
          at = floor(length(sel.idx)/2), line = -2)

    mtext('effect size', side=3, line=1, at=(mx/2), cex=0.7, adj=1)
    mtext('robustness', side=3, line=1, at=mx, cex=0.7, adj=0)
  } else {
    par(mar=c(0.1, 1.1, 0, 1.1))
    mx = max(-med+(abs(low.qt-upp.qt)))

    barplot(-med, horiz = TRUE, width=1, space=0, yaxs='i', col='gray30',
            xlim=c(0, mx),
            ylim=c(0, dim(relative.weights)[1]), xlab='', ylab='', yaxt='n')
    x<-par("usr")
    rect(x[1],x[3],x[2],x[4],col="lightgray")
    # grid lines
    grid(NULL, NA, lty=2, col='gray99')
    # plot the barplot again due to the idiotic way of how to change the background color of a plot in r
    barplot(-med, horiz = TRUE, width=1, space=0, yaxs='i', col='gray30',
            xlim=c(0, mx),
            ylim=c(0, dim(relative.weights)[1]), xlab='', ylab='', yaxt='n', add=TRUE)

    # error bars
    arrows(y0=c(1:length(med))-.5, x0=-upp.qt, x1=-low.qt, angle=90, code=3, length=.04)
    # labels
    mtext('median relative Gini coefficient', side=1, line=2, at=mx/2, cex=0.7, adj=0.5)
    mtext('effect size', side=3, line=1, at=(mx/2), cex=0.7, adj=.5)
  }

  box(lwd=1)
  if(verbose>2) cat('+ finished plot.feature.weights\n')
}

plot.pred.and.meta <- function(prediction, label, meta=NULL, verbose=0){
  if(verbose>2) cat("+ plot.pred.and.meta\n")
  par(mar=c(1.1, 4.1, 0.3, 5.1))
  img.data = as.matrix(prediction)
  colnames(img.data) = 'Classification result'
  if (!is.null(meta)) {
    img.data = cbind(meta[,dim(meta)[2]:1], img.data)
    ### transform any categorial column into a numeric one
    for (m in 1:dim(img.data)[2]) {
      img.data[,m] = (img.data[,m] - min(img.data[,m], na.rm=TRUE))
      if (max(img.data[,m], na.rm=TRUE) != 0) {
        img.data[,m] = img.data[,m] / max(img.data[,m], na.rm=TRUE)
      }
    }
  }

  grays = rev(gray(seq(0, 1, length.out=100)))
  image(as.matrix(img.data), col=grays,
        xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
  box(lwd=1)
  # add deliminator between the different classes
  abline(v=length(which(label@label == label@negative.lab))/length(label@label),
         col='red')

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
  if(verbose>2) cat('+ finished plot.pred.and.meta\n')
}

plot.proportion.of.weights <- function(selected.weights, all.weights,verbose=0){
  if(verbose>2) cat("+ plot.proportion.of.weights\n")
  par(mar=c(0.1, 6.1, 0, 1.1))
  boxplot(colSums(abs(selected.weights)) / colSums(abs(all.weights)), ylim=c(0,1))
  mtext('proportion of', side=1, line=1, at=1, adj=0.5, cex=0.7)
  mtext('weight shown', side=1, line=2, at=1, adj=0.5, cex=0.7)
  if(verbose>2) cat('+ finished plotting plot.proportion.of.weights\n')
}

plot.percentage.of.features <- function(selected.weights, all.weights, verbose=0){
  if(verbose>2) cat("+ plot.percentage.of.features\n")
  par(mar=c(0.1, 6.1, 0, 1.1))
  boxplot(dim(selected.weights)[1] / colSums(all.weights != 0), ylim=c(0,1))
  mtext('percentage of', side=1, line=1, at=1, adj=0.5, cex=0.7)
  mtext('features shown', side=1, line=2, at=1, adj=0.5, cex=0.7)
  if(verbose>2) cat('+ finished plot.percentage.of.features\n')
}

plot.heatmap <- function(image.data, limits, color.scheme, effect.size, verbose=0){
  if(verbose>2) cat("+ plot.heatmap\n")
  par(mar=c(0.1, 4.1, 0, 5.1))
  image(image.data, zlim=limits, col=color.scheme,
        xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
  if (!is.null(effect.size)){
    for (f in 1:dim(image.data)[2]){
      mtext(colnames(image.data)[f], side=4, line=1,
            at=(f-1)/(dim(image.data)[2]-1), las=2,
            cex=max(0.3, 0.8 - 0.01*dim(image.data)[2]),
            col=ifelse(effect.size[f]>0, color.scheme[1+4], color.scheme[length(color.scheme)-4]))
    }
  } else {
    for (f in 1:dim(image.data)[2]){
      mtext(colnames(image.data)[f], side=4, line=1,
            at=(f-1)/(dim(image.data)[2]-1), las=2,
            cex=max(0.3, 0.8 - 0.01*dim(image.data)[2]),
            col='black')
    }
  }
  box(lwd=1)
  if(verbose>2) cat('+ finished plot.heatmap\n')
}

#not really working atm
prepare.heatmap.fc <- function(heatmap.data, limits, sel.feat, meta=NULL, label, detect.lim, verbose=0){
  if(verbose>2) cat("+ prepare.heatmap.fc\n")
  if (!any(grepl('META', sel.feat))){
    img.data <- apply(heatmap.data[sel.feat,], 1, FUN=function(x, label, detect.lim){
          log10(x + detect.lim) - log10(median(as.numeric(x[label@n.idx])) + detect.lim)
    }, label=label, detect.lim=detect.lim)
  } else {
    img.data <- matrix(NA, nrow=length(sel.feat), ncol=ncol(heatmap.data))
    row.names(img.data) <- sel.feat
    for (f in sel.feat){
      print(f)
      if (!grepl('META', f)){
        median.ctr   <- median(as.numeric(heatmap.data[f,label@n.idx]))
        img.data[f,] <- log10(heatmap.data[f,] + detect.lim) - log10(median.ctr + detect.lim)
      } else {
        meta.data <- meta[,grep(strsplit(f, '_')[[1]][2], colnames(meta), ignore.case=TRUE, value=TRUE)]
        # transform metadata to zscores
        meta.data <- apply(meta.data,c(1,2),as.numeric)
        meta.data <- (meta.data - mean(meta.data, na.rm=TRUE))/sd(meta.data, na.rm=TRUE)
        img.data[f,] <- meta.data[colnames(heatmap.data)]
      }
    }
    img.data <- t(img.data)
  }
  img.data[img.data < limits[1]] = limits[1]
  img.data[img.data > limits[2]] = limits[2]
  if(verbose>2) cat('+ finished plot.heatmap\n')
  return(img.data)
}

prepare.heatmap.zscore <- function(heatmap.data, limits, verbose=0){
    cat("+ prepare.heatmap.zscore\n")
    # data is transposed and transformed to feature z-scores for display
    img.data <- apply(heatmap.data, 1, FUN=function(x){(x-mean(x))/sd(x)})
    img.data[img.data < limits[1]] <- limits[1]
    img.data[img.data > limits[2]] <- limits[2]
    if(verbose>2) cat('+ finished plot.heatmap\n')
    return(img.data)
}

select.features <- function(weights, model.type, consens.thres, norm.models, label, max.show, verbose=0){
  cat("+ select.features\n")
  # for linear models, select those that have been selected more than consens.thres percent of the models
  if (model.type != 'RandomForest'){
    # normalize by overall model size
    if (norm.models) {weights <- apply(weights, 2, function(x){x/sum(abs(x))})}
    sel.idx = which(rowSums(weights != 0) / dim(weights)[2] >= consens.thres)
    # normalize by model size and order features by relative model weight
    weights.norm <- apply(weights, 2, function(x){x/sum(abs(x))})
    med.weights <- apply(weights.norm, 1, median)
    median.sorted.features <- sort(med.weights[sel.idx], decreasing=TRUE, index.return=TRUE)
    # restrict to plot at maximum fifty features
    if (length(sel.idx) > max.show){
      warning("WARNING: restricting amount of features to be plotted to 50\n")
      median.sorted.features.abs <- sort(abs(med.weights),
                                         decreasing=TRUE, index.return=TRUE)
      idx <- head(median.sorted.features.abs$ix, n=max.show)
      median.sorted.features <- sort(med.weights[idx],
                                     decreasing=TRUE, index.return=TRUE)
      sel.idx <- idx[median.sorted.features$ix]
    } else {
      sel.idx = sel.idx[median.sorted.features$ix]
    }
  } else {
    # for Random Forest, caluclate relative median feature weights
    # and sort by auroc as effect measure
    weights <- apply(weights, 2, function(x){x/sum(abs(x))})
    median.sorted.features <- sort(apply(weights, 1, median),
            decreasing=FALSE, index.return=TRUE)
    # take the feature with median higher than consens.threshold
    sel.idx <- median.sorted.features$ix[which(median.sorted.features$x >= consens.thres)]

    if (length(sel.idx) > max.show){
      sel.idx <- tail(sel.idx, n=max.show)
    }
    # sort again by auroc
    # auroc.effect <-  sapply(row.names(feat)[sel.idx],
    #    FUN=function(f, feat, label){roc(predictor=feat[f,],
    #                    response=label@label, direction='<')$auc}, feat=feat, label=label)
    #
    # if (length(sel.idx) > 50){
    #   cat("Restricting amount of features to be plotted to 50\n")
    #   auroc.transformed <- auroc.effect
    #   auroc.transformed[auroc.transformed <= .5] <- 1 - auroc.transformed[auroc.transformed <= .5]
    #   idx <- head(sort(auroc.transformed, decreasing=TRUE, index.return=TRUE)$ix, n=50)
    #   auroc.effect <- auroc.effect[idx]
    # }
    # sel.idx <- sel.idx[sort(auroc.effect, decreasing=TRUE, index.return=TRUE)$ix]
  }

  if(verbose>2) cat('+++ generating plot for a model with', length(sel.idx), 'selected features\n')
  if(verbose>2) cat('+ finished select.features\n')
  return(sel.idx)
}

get.weights.matrix <- function(models.list, verbose=0){
  if(verbose>2) cat("+ get.weights.matrix\n")
  W.mat <- as.numeric(models.list[[1]]$feat.weights)
  for(i in 2:length(models.list)){
    W.mat <- cbind(W.mat,as.numeric(models.list[[i]]$feat.weights))
  }
  rownames(W.mat) <- models.list[[1]]$features
  colnames(W.mat) <- paste('M', 1:ncol(W.mat), sep='_')
  if(verbose>2) cat('+ finished get.weights.matrix\n')
  return(W.mat)
}
