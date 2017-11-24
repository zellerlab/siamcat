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

#' @title Check and visualize associations between features and classes
#' @description This function calculates for each feature the median log10 fold
#' change between the different classes found in labels.
#'
#' Significance of the differences is computed for each feature using a wilcoxon
#' test followed by multiple hypothesis testing correction.
#'
#' Additionally, the Area under the Receiver Operating Characteristic Curve
#' '(AU-ROC) is computed for the features found to be associated with the two
#' different classes at a user-specified significance level \code{alpha}.
#'
#' Finally, the function produces a plot of the top \code{max.show} associated
#' features, showing the distribution of the log10-transformed abundances for
#' both classes, the median fold change, and the AU-ROC.
#' @param feat feature object
#' @param label label object
#' @param fn.plot filename for the pdf-plot
#' @param color.scheme valid R color scheme, defaults to \code{'RdYlBu'}
#' @param alpha float, significance level, defaults to \code{0.05}
#' @param mult.corr multiple hypothesis correction method, see \code{\link[stats]{p.adjust}}, defaults to \code{"fdr"}
#' @param sort.by string, sort features by p-value (\code{"pv"}) or by fold change (\code{"fc"}), defaults to \code{"pv"}
#' @param detect.lim float, pseudocount to be added before log-transormation of the data, defaults to \code{1e-08}
#' @param max.show integer, how many associated features should be shown, defaults to \code{50}
#' @param plot.type string, specify how the abundance should be plotted, must be one of these: \code{c("bean", "box", "quantile.box", "quantile.rect")}, defaults to \code{"bean"}
#' @return Does not return anything, but produces an association plot
#' @keywords SIAMCAT check.associations
#' @export
check.associations <- function(feat, label, fn.plot, color.scheme="RdYlBu",
<<<<<<< HEAD
                               alpha=0.05, min.fc=0, mult.corr="fdr", sort.by="pv",
                               detect.lim=1e-08, max.show=50, plot.type="bean"){
=======
                               alpha=0.05, mult.corr="fdr", sort.by="pv",
                               detect.lim=10^-8, max.show=50, plot.type="bean"){
>>>>>>> 3a7a4f2007aaa43d62c94673c1c9e1f05dc225f7


  ### some color pre-processing
  if (!color.scheme %in% row.names(brewer.pal.info)){
    warning("Not a valid RColorBrewer palette name, defaulting to RdYlBu...\n
    See brewer.pal.info for more information about RColorBrewer palettes...")
    color.scheme <- 'RdYlBu'
  }
  color.scheme <- rev(colorRampPalette(brewer.pal(brewer.pal.info[color.scheme,'maxcolors'], color.scheme))(100))
  
  col.p <- color.scheme[length(color.scheme)-4]
  col.n <- color.scheme[1+4]

  ### Define set of vectors that have the indeces and "description" of all positively and negatively labeled training examples.
  p.val <- vector('numeric', nrow(feat))
  fc    <- vector('numeric', nrow(feat))

  
  ### Calculate wilcoxon and FC for each feature
  for (i in 1:nrow(feat)) {
    fc[i]    <- median(log10(feat[i,label$p.idx] + detect.lim)) - median(log10(feat[i,label$n.idx] + detect.lim))
    p.val[i] <- wilcox.test(feat[i,label$n.idx], feat[i,label$p.idx], exact = FALSE)$p.value
  }

  ### Apply multi-hypothesis testing correction
  if(!tolower(mult.corr) %in% c('none','bonferroni','holm','fdr','bhy')) {
    stop("Unknown multiple testing correction method:', mult.corr,' Stopping!\n
          Must of one of c('none','bonferroni','holm','fdr','bhy')\n")
  }
  if (mult.corr == 'none') {
    warning('No multiple hypothesis testing performed...')
    p.adj <- p.val
  } else {
    p.adj <- p.adjust(p.val, method=tolower(mult.corr))
  }
  cat('Found', sum(p.adj < alpha, na.rm=TRUE), 'significant associations at a significance level <', alpha, '\n')

  idx <- which(p.adj < alpha)

  # if (min.fc > 0) {
  #   idx <- which(p.adj < alpha & abs(fc) > min.fc)
  #   cat('Found', length(idx), 'significant associations with absolute log10 fold change >', min.fc, '\n')
  # }
  ### Stop if no significant features were found
  if (length(idx) == 0){
    stop('No significant associations found. Stopping...\n')
  }

  ### Sort features
  if (sort.by == 'fc') {
    idx <- idx[order(fc[idx], decreasing=FALSE)]
  } else if (sort.by == 'pv') {
    idx <- idx[order(p.adj[idx], decreasing=TRUE)]
  } else {
    cat('Unknown sorting option:', sort.by, 'order by p-value...\n')
    idx <- idx[order(p.adj[idx], decreasing=TRUE)]
  }
  # Really needed?
  # for (i in idx) {
  #   cat(sprintf('%-40s', rownames(feat)[i]), 'p-value:', format(p.adj[i], digits=4), '\n')
  # }
  # # truncated the list for the following plots
  if (length(idx) > max.show) {
    idx <- idx[(length(idx)-max.show+1):length(idx)]
    cat('Truncating the list of significant associations to the top', max.show, '\n')
  }


  ### compute single-feature AUCs
  cat('\nCalculating the area under the ROC curve for each significantly associated feature\n')
  aucs <- vector('numeric', nrow(feat))
  for (i in idx) {
    f       <- feat[i,]
    ev      <- eval.classifier(f, label$label, label)
    aucs[i] <- calc.auroc(ev)
    if (aucs[i] < 0.5) {
      aucs[i] <- 1-aucs[i]
    }
  }

  # for (i in idx) {
  #   cat(sprintf('%-40s', rownames(feat)[i]), aucs[i], '\n')
  # }

    ### generate plots with significant associations between features and labels
    pdf(fn.plot, paper='special', height=8.27, width=11.69) # format: A4 landscape

    lmat  <- cbind(1,2,3,4)
    layout(lmat, widths=c(0.6,0.075,0.2,0.2))

    x <- log10(as.matrix(feat[idx, label$p.idx, drop=FALSE]) + detect.lim)
    y <- log10(as.matrix(feat[idx, label$n.idx, drop=FALSE]) + detect.lim)

    col <- c(paste(col.n, '77', sep=''), paste(col.p, '77', sep=''), 'gray')
    if (plot.type == 'box') {
      par(mar=c(5.1, 25.1, 4.1, 0))
      box.colors <- rep(c(col[1],col[2]),nrow(x))
      plot.data <- data.frame()
      for (i in 1:nrow(x)){
        temp <- as.data.frame(rbind(cbind(x[i,],rep(paste(label$n.lab, rownames(x)[i]), length(x[i,]))), cbind(y[i,], rep(paste(label$p.lab, rownames(x)[i]), length((y[i,]))))))
        temp[,1] <- as.numeric(as.character(temp[,1]))
        plot.data <- rbind(plot.data, temp)
        if (i == nrow(x)) {
          plot(NULL, xlab='', ylab='',xaxs='i', yaxs='i', axes=FALSE,
               xlim=c(min(plot.data[,1]-0.2), max(plot.data[,1]) + 1), ylim=c(+0.5, length(idx)*2+0.5), type='n')
          boxplot(plot.data[,1] ~ plot.data[,ncol(plot.data)],horizontal=TRUE,
                  names = c(""), show.names = FALSE, col = box.colors, axes = FALSE, outcol = c(col[1], col[2]), add = TRUE)
          mn          <- as.integer(c(min(plot.data[,1])))
          mx          <- as.integer(c(max(plot.data[,1])))
          ticks       <- mn:mx
          for (v in ticks) {
            abline(v=v, lty=3, col='lightgrey')
          }
          tick.labels <- formatC(10^ticks, format='E', digits=0)
          axis(side=1, at=ticks, labels=tick.labels, cex.axis=0.7)
          ### function label.plot.horizontal has been written in utils.r.
          label.plot.horizontal(x, y, rownames(feat)[idx], x.suff=paste(' (', label$p.lab, ')', sep=''),
                                y.suff=paste(' (', label$n.lab, ')', sep=''), outer.diff = 2, inner.diff.x = 0, inner.diff.y = -1)
        }
      }
    }
    else if (plot.type == "quantile.box"){
      par(mar=c(5.1, 25.1, 4.1, 0))
      plot.data.range(x, y, rownames(feat)[idx], x.col=col[2], y.col=col[1])
      label.plot.horizontal(x, y, rownames(feat)[idx], x.suff=paste(' (', label$p.lab, ')', sep=''),
                            y.suff=paste(' (', label$n.lab, ')', sep=''), outer.diff = 1, inner.diff.x = 0.15, inner.diff.y = -0.15)
    }

    else if (plot.type == "quantile.rect"){
      par(mar=c(5.1, 25.1, 4.1, 0))
      quantiles.vector <- c(0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9)
      x.q = apply(x, 1, function (x) quantile(x, quantiles.vector, na.rm=TRUE, names=FALSE))
      x.medians = apply(x,1,function (x) median(x))
      y.q = apply(y, 1, function (y) quantile(y, quantiles.vector, na.rm=TRUE, names=FALSE))
      y.medians = apply(y,1,function (y) median(y))

      p.m = min(c(min(x, na.rm=TRUE), min(y, na.rm=TRUE)))
      plot(rep(p.m, dim(x)[1]), 1:dim(x)[1],
           xlab='', ylab='', yaxs='i', axes=FALSE,
           xlim=c(min(x,y), max(x,y+2)), ylim=c(0, dim(x)[1]), frame.plot=FALSE, type='n')
      for (v in seq(p.m,0,1)) {
        abline(v=v, lty=3, col='lightgrey')
      }

      tck = floor(p.m):0
      axis(1, tck, formatC(10^tck, format='E', digits=0), las=1, cex.axis=0.7)
      for (i in 1:(nrow(x.q)/2)){
        if (i == 1) {
          rect(x.q[i,], 0.5:dim(x)[1], x.q[nrow(x.q)+1-i,], (0.5:dim(x)[1])+0.3, col = c("white"), border = c("black"), lwd = 0.9)
          rect(y.q[i,], 0.5:dim(y)[1], y.q[nrow(y.q)+1-i,], (0.5:dim(y)[1])-0.3, col = c("white"), border = c("black"), lwd = 0.9)

        }
        else {
          rect(x.q[i,], 0.5:dim(x)[1], x.q[nrow(x.q)+1-i,], (0.5:dim(x)[1])+0.3, col = col[2], border = c("black"), lwd = 0.9)
          rect(y.q[i,], 0.5:dim(y)[1], y.q[nrow(y.q)+1-i,], (0.5:dim(y)[1])-0.3, col = col[1], border = c("black"), lwd = 0.9)
        }
      }
      points(x.medians, y=(0.5:dim(x)[1])+0.15, pch=18, cex = min(35/nrow(x),4))
      points(y.medians, y=(0.5:dim(y)[1])-0.15, pch=18, cex = min(35/nrow(x),4))
      mtext('Quantiles', 3, line=0, at=1, adj = 1.675, padj = 0.45, las=1, cex=0.7)
      ### create.tints.rgb is in utils.r
      red.tints  <- create.tints.rgb(col2rgb(col[2])/255, nr.tints=5, tint.steps = 1/5)
      blue.tints <- create.tints.rgb(col2rgb(col[1])/255, nr.tints=5, tint.steps = 1/5)
      legend(-1.75, nrow(x), legend = c("80%","60%","40%","20%","median","","","","",""),
             bty='n', cex=1, fill=c(
               rgb(matrix(red.tints[,5], ncol = 3)),
               rgb(matrix(red.tints[,3], ncol = 3)),
               rgb(matrix(red.tints[,2], ncol = 3)),
               rgb(matrix(red.tints[,1], ncol = 3)),
               0,
               rgb(matrix(blue.tints[,5], ncol = 3)),
               rgb(matrix(blue.tints[,3], ncol = 3)),
               rgb(matrix(blue.tints[,2], ncol = 3)),
               rgb(matrix(blue.tints[,1], ncol = 3)),
               0),
             lty    <- c(0,0,0,0,0,0,0,0,0,0),
             lwd    <- c(1.3,1.3,1.3,1.3,2,1.3,1.3,1.3,1.3,1.3), ncol = 2,
             border <- c("black", "black","black","black","white","black","black","black","black","white"))
      legend(-1.675, nrow(x), legend = c("","","","",""),
             bty='n', lty = c(0,0,0,0,0),
             # cap legend size for diamond (should look symmetric to other symbols)
             pch = 18, cex = 1, pt.cex = c(0,0,0,0, min(35/nrow(x), 2.25)))

      label.plot.horizontal(x, y, rownames(feat)[idx], x.suff=paste(' (', label$p.lab, ')', sep=''),
                            y.suff=paste(' (', label$n.lab, ')', sep=''), outer.diff = 1, inner.diff.x = -0.3, inner.diff.y = -0.6)
    }
    else if (plot.type == "bean"){
      par(mar=c(5.1, 25.1, 4.1, 0))
      bean.data <- data.frame()
      for (i in 1:nrow(x)){
        temp      <- as.data.frame(rbind(cbind(x[i, ], rep(paste(label$n.lab, rownames(x)[i]), length(x[i, ]))),
                                    cbind(y[i, ], rep(paste(label$p.lab, rownames(x)[i]), length((y[i, ]))))))
        temp[,1]  <- as.numeric(as.character(temp[,1]))
        bean.data <- rbind(bean.data, temp)
        if (i == nrow(x)){
          plot(NULL, xlab='', ylab='',xaxs='i', yaxs='i', axes=FALSE,
               xlim = c(as.integer(min(x))-1.5,as.integer(max(x))+1), ylim=c(0.45, length(idx)+0.6), type='n')
          beanplot(bean.data[, 1] ~ bean.data[, ncol(bean.data)], side = "both", bw="nrd0", col = list(col[1],
                                                                                                       col[2]), horizontal = TRUE, names = c(""), show.names = FALSE, beanlines = "median", maxstripline = 0.2, what = c(FALSE,TRUE,TRUE,FALSE),
                   axes = FALSE, add = TRUE )
          mn    <- as.integer(c(min(bean.data[,1])-1.5))
          mx    <- as.integer(c(max(bean.data[,1])+1))
          ticks <- mn:mx
          for (v in ticks) {
            abline(v=v, lty=3, col='lightgrey')
          }
          tick.labels <- formatC(10^ticks, format='E', digits=0)
          axis(side=1, at=ticks, labels=tick.labels, cex.axis=0.7)
          label.plot.horizontal(x, y, rownames(feat)[idx], x.suff=paste(' (', label$p.lab, ')', sep=''),
                                y.suff=paste(' (', label$n.lab, ')', sep=''), outer.diff = 1, inner.diff.x = 0.15, inner.diff.y = -0.15)
        }
      }
    }
    else {
      print("plot type has not been specified properly; continue with quantileplot")
      plot.type <- "quantile.box"
      par(mar=c(5.1, 25.1, 4.1, 0))
      plot.data.range(x, y, rownames(feat)[idx], x.col=col[2], y.col=col[1])
      label.plot.horizontal(x, y, rownames(feat)[idx], x.suff=paste(' (', label$p.lab, ')', sep=''),
                            y.suff=paste(' (', label$n.lab, ')', sep=''), outer.diff = 1, inner.diff.x = 0.15, inner.diff.y = -0.15)
    }

    p.val.annot <- formatC(p.adj[idx], format='E', digits=2)
    if (sum(p.adj < alpha, na.rm=TRUE) <= max.show) {
      title(main='Differentially abundant features', xlab='Abundance (log10-scale)')
    } else {
      title(main=paste('Differentially abundant features\ntruncated to the top', max.show),
            xlab='Abundance (log10-scale)')
    }
    par(mar=c(5.1,0,4.1, 0))
    for (i in 1:length(p.val.annot)) {
      if (plot.type == 'box'){
        mtext(p.val.annot[i], 4, line=2, at=(2*i)-0.5, las=1, cex=min(0.7, 1-(length(idx)/100)))
      }
      else if (plot.type == "quantile.rect"){
        mtext(p.val.annot[i], 4, line=2, at=i-0.5, las=1, cex=min(0.7, 1-(length(idx)/100)))
      }
      else {
        mtext(p.val.annot[i], 4, line=2, at=i, las=1, cex=min(0.7, 1-(length(idx)/100)))
      }
    }
    plot(NULL, xlab='', ylab='',xaxs='i', yaxs='i', axes=FALSE,
         type='n', xlim=c(0,10), ylim=c(0,length(p.val.annot)+0.5))
    title(main='Adj. p-value')

    # plot fold changes
    par(mar=c(5.1, 2.1, 4.1, 2.1))
    bcol  <- ifelse(fc[idx] > 0, col[2], col[1])
    mn    <- floor(min(fc[idx]))
    mx    <- ceiling(max(fc[idx]))
    mx    <- max(abs(mn), abs(mx))
    if (!is.finite(mx)) {
      mx    <- 10
    }
    mn    <- -mx
    plot(NULL, xlab='', ylab='', xaxs='i', yaxs='i', axes=FALSE,
         xlim=c(mn, mx), ylim=c(0.2, length(idx)+0.2), type='n')
    barplot(fc[idx], horiz=TRUE, width=0.6, space=2/3, col=bcol, axes=FALSE, add=TRUE)

    ticks <- mn:mx
    for (v in ticks) {
      abline(v=v, lty=3, col='lightgrey')
    }
    tick.labels <- formatC(10^ticks, format='E', digits=0)
    axis(side=1, at=ticks, labels=tick.labels, cex.axis=0.7)
    title(main='Fold change', xlab='FC (log10-scale)')

    # plot single-feature AUCs
    par(mar=c(5.1, 1.1, 4.1, 3.1))
    plot(NULL, xlab='', ylab='',xaxs='i', yaxs='i', axes=FALSE,
         xlim=c(0.5,1), ylim=c(0.5, length(idx)+0.5), type='n')
    ticks       <- seq(0.5, 1.0, length.out=6)
    for (v in ticks) {
      abline(v=v, lty=3, col='lightgrey')
    }
    for (b in 1:length(idx)) {
      i <- idx[b]
      points(aucs[i], b, pch=18, col=bcol[b])
      points(aucs[i], b, pch=5, col='black', cex=0.9)
    }
    axis(side=1, at=ticks, cex.axis=0.7)
    title(main='Feature AUCs', xlab='AU-ROC')
    # close pdf device
    tmp <- dev.off()


}

### label.plot.horizontal() takes as input lists of (significantly) differentially abundant bacterial features and plots their names
### on the left side of a figur, parallel to each associated plot. inner.diff.x, inner.diff.y and outer.diff are numerical values that can be
### used to tweak the position of the text lines relatively to their plot. Specifically, inner.diff.y and inner.diff.x will shift the text
###  alongside the y-axis. outer.diff on the other hand is used as a multiplication factor which changes the distance between each different
### feature example combination globally.
label.plot.horizontal <- function(x, y, labels = NULL, x.suff, y.suff, inner.diff.x = NULL, inner.diff.y = NULL, outer.diff = NULL){
  stopifnot(length(labels) == dim(c)[1])
  if (!is.null(y) && !is.null(x.suff) && !is.null(y.suff)) {
    for (i in 1:dim(x)[1]) {
      mtext(paste(labels[i], x.suff), 2, line=1, at=i*outer.diff+inner.diff.x, las=1, cex=min(0.7, 1-(nrow(x)/70)))
      mtext(y.suff, 2, line=1, at=i*outer.diff+inner.diff.y, las=1, cex=min(0.7, 1-(nrow(x)/70)))
    }
  } else {
    for (i in 1:dim(x)[1]) {
      mtext(labels[i], 2, line=1, at=i*outer.diff+inner.diff.x, las=1, cex=min(0.7, 1-(nrow(x)/50)))
    }
  }
}

##### function to create different tints of a color based on the color's rgb specifications. Each column specifies one tint.
### Make sure that you specify your rgb values as rgb/255! Also, DO NOT call rgb() on your color vector!
create.tints.rgb <- function(color.rgb, nr.tints, tint.steps = 1/nr.tints) {
  tints <- matrix(rep(0,(3*(nr.tints))),nrow=3, ncol=nr.tints)
  for (i in 1:nr.tints){
    tints[1,i] = color.rgb[1] + (1 - color.rgb[1]) * (tint.steps*i)
    tints[2,i] = color.rgb[2] + (1 - color.rgb[2]) * (tint.steps*i)
    tints[3,i] = color.rgb[3] + (1 - color.rgb[3]) * (tint.steps*i)
  }
  return (tints)
}

### TODO docu!
# # # #' @export
plot.data.range <- function(x, y=NULL, x.col='black', y.col='black', labels=NULL, x.suff=NULL, y.suff=NULL) {
  if (is.null(y)) {
    p.m = min(x, na.rm=TRUE)
  } else {
    stopifnot(dim(x)[1] == dim(y)[1])
    p.m = min(c(min(x, na.rm=TRUE), min(y, na.rm=TRUE)))
  }
  plot(rep(p.m, dim(x)[1]), 1:dim(x)[1],
       xlab='', ylab='', yaxs='i', axes=FALSE,
       xlim=c(p.m, 0), ylim=c(0.5, dim(x)[1]+0.5), frame.plot=FALSE, type='n')
  for (v in seq(p.m,-1,1)) {
    abline(v=v, lty=3, col='lightgrey')
  }

  tck = floor(p.m):0
  axis(1, tck, formatC(10^tck, format='E', digits=0), las=1, cex.axis=0.7)

  x.q = apply(x, 1, function (x) quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm=TRUE, names=FALSE))
  if (is.null(y)) {
    # inter-quartile range
    rect(x.q[2,], (1:dim(x)[1])-0.2, x.q[4,], (1:dim(x)[1])+0.2)
    # 90% interval
    segments(x.q[1,], 1:dim(x)[1], x.q[2,], 1:dim(x)[1])
    segments(x.q[4,], 1:dim(x)[1], x.q[5,], 1:dim(x)[1])
    segments(x.q[1,], y0=(1:dim(x)[1])-0.15, y1=(1:dim(x)[1])+0.15)
    segments(x.q[5,], y0=(1:dim(x)[1])-0.15, y1=(1:dim(x)[1])+0.15)
    # median
    segments(x.q[3,], y0=(1:dim(x)[1])-0.2, y1=(1:dim(x)[1])+0.2, lwd=2)
    # scatter plot on top
    for (i in 1:dim(x)[1]) {
      if (nchar(x.col) > 7) {
        # adjust alpha channel by reducing transparency
        a = substr(x.col,nchar(x.col)-1, nchar(x.col))
        a = 1 - (1 - as.numeric(paste('0x', a, sep=''))/255)/2
        x.col = gsub('..$', toupper(as.hexmode(round(a*255))), x.col)
      }
      points(x[i,], rep(i, dim(x)[2])+rnorm(ncol(x),sd=0.05), pch=16, cex=0.6, col=x.col)
    }
  } else {
    y.q = apply(y, 1, function (x) quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm=TRUE, names=FALSE))
    # inter-quartile range
    rect(x.q[2,], 1:dim(x)[1], x.q[4,], (1:dim(x)[1])+0.3, col=x.col)
    rect(y.q[2,], (1:dim(y)[1])-0.3, y.q[4,], 1:dim(y)[1], col=y.col)
    # 90% interval
    segments(x.q[1,], 1:dim(x)[1], x.q[5,], 1:dim(x)[1])#, col=x.col)
    segments(y.q[1,], 1:dim(x)[1], y.q[5,], 1:dim(x)[1])#, col=x.col)
    segments(x.q[1,], y0=1:dim(x)[1], y1=(1:dim(x)[1])+0.2)
    segments(y.q[1,], y0=(1:dim(x)[1])-0.2, y1=1:dim(x)[1])
    segments(x.q[5,], y0=1:dim(x)[1], y1=(1:dim(x)[1])+0.2)
    segments(y.q[5,], y0=(1:dim(x)[1])-0.2, y1=1:dim(x)[1])
    # median
    segments(x.q[3,], y0=1:dim(x)[1], y1=(1:dim(x)[1])+0.3, lwd=3)#, col=x.col)
    segments(y.q[3,], y0=(1:dim(x)[1])-0.3, y1=1:dim(x)[1], lwd=3)#, col=y.col)
    # scatter plot on top
    for (i in 1:dim(x)[1]) {
      if (nchar(x.col) > 7) {
        # adjust alpha channel by reducing transparency
        a = substr(x.col,nchar(x.col)-1, nchar(x.col))
        a = 1 - (1 - as.numeric(paste('0x', a, sep=''))/255)/2
        x.col = gsub('..$', toupper(as.hexmode(round(a*255))), x.col)
      }
      if (nchar(y.col) > 7) {
        # adjust alpha channel by reducing transparency
        a = substr(y.col,nchar(y.col)-1, nchar(y.col))
        a = 1 - (1 - as.numeric(paste('0x', a, sep=''))/255)/2
        y.col = gsub('..$', toupper(as.hexmode(round(a*255))), y.col)
      }
      points(x[i,], rep(i+0.15, ncol(x))+rnorm(ncol(x),sd=0.03), pch=16, cex=0.6, col=x.col)
      points(y[i,], rep(i-0.15, ncol(y))+rnorm(ncol(y),sd=0.03), pch=16, cex=0.6, col=y.col)
    }
  }
}
