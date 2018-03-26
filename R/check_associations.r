#!/usr/bin/Rscript
###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# R flavor
# EMBL Heidelberg 2012-2018
# GNU GPL 3.0
###

#' @title Check and visualize associations between features and classes
#' @description This function calculates for each feature a pseudo-fold change
#' (geometrical mean of the difference between quantiles) between the different
#' classes found in labels.
#'
#' Significance of the differences is computed for each feature using a Wilcoxon
#' test followed by multiple hypothesis testing correction.
#'
#' Additionally, the Area Under the Receiver Operating Characteristic Curve
#' (AU-ROC) and a prevalence shift are computed for the features found to be
#' associated with the two different classes at a user-specified
#' significance level \code{alpha}.
#'
#' Finally, the function produces a plot of the top \code{max.show} associated
#' features, showing the distribution of the log10-transformed abundances for
#' both classes, and user-selected panels for the effect (AU-ROC, Prevalence
#' Shift, and Fold Change)
#' @param siamcat object of class \link{siamcat-class}
#' @param fn.plot filename for the pdf-plot
#' @param color.scheme valid R color scheme or vector of valid R colors (must
#'        be of the same length as the number of classes), defaults to \code{'RdYlBu'}
#' @param alpha float, significance level, defaults to \code{0.05}
#' @param mult.corr multiple hypothesis correction method, see \code{\link[stats]{p.adjust}},
#'        defaults to \code{"fdr"}
#' @param sort.by string, sort features by p-value (\code{"p.val"}), by fold change
#'        (\code{"fc"}) or by prevalence shift (\code{"pr.shift"}), defaults to
#'        \code{"fc"}
#' @param detect.lim float, pseudocount to be added before log-transformation of
#'        the data, defaults to \code{NULL} (estimated as the 5% quantile)
#' @param pr.cutoff float, cutoff for the prevalence computation, defaults to
#'        \code{1e-06}
#' @param max.show integer, how many associated features should be shown,
#'        defaults to \code{50}
#' @param plot.type string, specify how the abundance should be plotted, must be
#'        one of these: \code{c("bean", "box", "quantile.box", "quantile.rect")},
#'        defaults to \code{"quantile.box"}
#' @param panels vector, name of the panels to be plotted next to the log10-
#'        transformed abundances, possible entries are \code{c("fc", "auroc",
#'        "prevalence")}, defaults to \code{c("fc", "auroc")}
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'        for only information about progress and success, \code{2} for normal
#'        level of information and \code{3} for full debug information, defaults to \code{1}
#' @return Does not return anything, but produces an association plot
#' @keywords SIAMCAT check.associations
#' @export
#' @examples
#'  # Simple working example
#'  check.associations(siamcat, './assoc_plot.pdf')
#'
#'  # Plot associations as bean plot
#'  check.associations(siamcat, './assoc_plot_bean.pdf', plot.type='bean')
#'
#'  # Plot assocations as box plot
#'  # Additionally, sort by p-value instead of by fold change
#'  check.associations(siamcat, './assoc_plot_fc.pdf', plot.type='box', sort.by='p.val')
#'
#'  # Custom colors
#'  check.associations(siamcat, './assoc_plot_blue_yellow.pdf', plot.type='box',
#'    color.scheme=c('cornflowerblue', '#ffc125'))
check.associations <- function(siamcat, fn.plot, color.scheme="RdYlBu",
                               alpha=0.05, mult.corr="fdr", sort.by="fc",
                               detect.lim=NULL, pr.cutoff=10^-6, max.show=50,
                               plot.type="quantile.box", panels=c("fc", "auroc"), verbose=1){
  # check panel and plot.type parameter
  if(verbose>1) cat("+ starting check.associations\n")
  s.time <- proc.time()[3]

  if (!all(panels %in% c("fc", "auroc", "prevalence"))){
    stop("Unknown panel-type selected!")
  }
  if (length(panels) > 3){
    warning("Plot layout is not suited for more than 3 panels. Continuing with first three panels.")
    panels <- panels[1:3]
  }
  if ((!plot.type %in% c("bean", "box", "quantile.box", "quantile.rect")) || length(plot.type) != 1){
    warning("Plot type has not been specified properly! Continue with quantile.box.")
    plot.type <- "quantile.box"
  }
  # either give n_classes colors or color palette
  col <- check.color.scheme(color.scheme, siamcat@label)

  feat <- matrix(siamcat@phyloseq@otu_table,nrow=nrow(siamcat@phyloseq@otu_table), ncol=ncol(siamcat@phyloseq@otu_table),
                 dimnames = list(rownames(siamcat@phyloseq@otu_table), colnames(siamcat@phyloseq@otu_table)))
  ### Calculate different effect sizes
  if(verbose>2) cat("+++ analysing features\n")
  result.list <- analyse.binary.marker(feat=feat, label=siamcat@label, detect.lim=detect.lim, colors=col,
                                        pr.cutoff=pr.cutoff, mult.corr=mult.corr, alpha=alpha,
                                        max.show=max.show, sort.by=sort.by, probs.fc=seq(.1, .9, .05),verbose=verbose)

  ### TODO: remove at some point
  p.val     <- result.list$p.val
  p.adj     <- result.list$p.adj
  fc        <- result.list$fc
  aucs      <- result.list$aucs
  pr.shift  <- result.list$pr.shift
  feat.red  <- result.list$feat.red
  truncated <- result.list$truncated
  bcols     <- result.list$bcol
  detect.lim <- result.list$detect.lim

  ##############################################################################
  ### generate plots with significant associations between features and labels

  # make plot matrix dependent on panels parameters
  if(verbose>2) cat("+++ preparing plotting layout\n")
  if (length(panels) == 3){
    layout.mat <- cbind(2,1, t(seq(3, length.out=length(panels))))
    widths <- c(0.5, 0.1, rep(0.4/3, length(panels)))
  } else {
    layout.mat <- cbind(2,1, t(seq(3, length.out=length(panels))))
    widths <- c(0.5, 0.1, rep(0.2, length(panels)))
  }
  pdf(fn.plot, paper='special', height=8.27, width=11.69) # format: A4 landscape

  layout(mat=layout.mat, widths=widths)

  ##############################################################################
  # PANEL 2: P-VALUES
  # print p-values in second panel of the plot
  associations.pvals.plot(p.vals=p.adj, alpha=alpha,verbose=verbose)

  ##############################################################################
  # PANEL 1: DATA
  # prepare margins
  associations.margins.plot(species_names = row.names(feat.red),verbose=verbose)

  # get data
  data.n <- log10(as.matrix(feat.red[, siamcat@label@n.idx, drop=FALSE]) + detect.lim)
  data.p <- log10(as.matrix(feat.red[, siamcat@label@p.idx, drop=FALSE]) + detect.lim)

  if(verbose>2) cat("+++ plotting results\n")
  if (plot.type == "bean"){
    associations.bin.plot(data.n, data.p, siamcat@label, col=col, verbose=verbose)
  } else if (plot.type == "box"){
    associations.box.plot(data.n, data.p, siamcat@label, col=col, verbose=verbose)
  } else if (plot.type == "quantile.box"){
    associations.quantile.box.plot(data.p, data.n, siamcat@label, col=col, verbose=verbose)
  } else if (plot.type == "quantile.rect"){
    associations.quantile.rect.plot(data.p, data.n, siamcat@label, col=col, verbose=verbose)
  }

  # plot title
  if (!truncated) {
    title(main='Differentially abundant features',
          xlab='Abundance (log10-scale)')
  } else {
    title(main=paste('Differentially abundant features\nshowing top', max.show, 'features'),
          xlab='Abundance (log10-scale)')
  }

  ##############################################################################
  # OTHER PANELS
  for (p in panels){
    if (p == "fc"){
      associations.fcs.plot(fc.all=fc, binary.cols=bcols,verbose=verbose)
    } else if (p == "prevalence"){
      associations.pr.shift.plot(pr.shifts=pr.shift, col=col,verbose=verbose)
    } else if (p == "auroc"){
      associations.aucs.plot(aucs=aucs, binary.cols=bcols,verbose=verbose)
    }
  }

  # close pdf device
  tmp <- dev.off()
  e.time <- proc.time()[3]
  if(verbose>1) cat("+ finished check.associations in",e.time-s.time,"s\n")
  if(verbose==1) cat("Plotted associations between features and label successfully to:",fn.plot,"\n")
}


### one function for each type of plot
# bean plot
associations.bin.plot <- function(data1, data2, label, col, verbose=1){
  if(verbose>2) cat("+ starting associations.bin.plot\n")
  # create data.frame in format for beanplot
  bean.data <- data.frame()
  for (i in 1:nrow(data2)){
    temp      <- as.data.frame(rbind(cbind(data1[i, ], rep(paste(label@p.lab, rownames(data1)[i]), length(data1[i, ]))),
                                     cbind(data2[i, ], rep(paste(label@n.lab, rownames(data2)[i]), length((data2[i, ]))))))
    temp[,1]  <- as.numeric(as.character(temp[,1]))
    bean.data <- rbind(bean.data, temp)
  }

  plot(NULL, xlab='', ylab='',xaxs='i', yaxs='i', axes=FALSE,
       xlim = c(as.integer(min(data2))-1.5,as.integer(max(data2))+1),
       ylim=c(0.45, nrow(data1)+0.6), type='n')

  beanplot(bean.data[, 1] ~ bean.data[, ncol(bean.data)],
           side = "both", bw="nrd0", col = list(col[1], col[2]),
           horizontal = TRUE, names = c(""), show.names = FALSE,
           beanlines = "median", maxstripline = 0.2,
           what = c(FALSE,TRUE,TRUE,FALSE), axes = FALSE, add = TRUE )
  mn    <- as.integer(c(min(bean.data[,1])-1.5))
  mx    <- as.integer(c(max(bean.data[,1])+1))
  ticks <- mn:mx
  for (v in ticks) {
    abline(v=v, lty=3, col='lightgrey')
  }
  tick.labels <- formatC(10^ticks, format='E', digits=0)
  axis(side=1, at=ticks, labels=tick.labels, cex.axis=0.7)
  legend('topright', legend=c(label@p.lab, label@n.lab), fill=rev(col), bty='n')
  associations.labels.plot(row.names(data1), plot.type='bean', verbose=verbose)
  if(verbose>2) cat("+ finished associations.bin.plot\n")
}

# box plot
associations.box.plot <- function(data1, data2, label, col, verbose=1){
  if(verbose>2) cat("+ starting associations.box.plot\n")
  box.colors <- rep(c(col[1],col[2]),nrow(data1))

  plot.data <- data.frame()
  for (i in 1:nrow(data1)){
    temp <- as.data.frame(rbind(cbind(data2[i,],rep(paste(label@n.lab, rownames(data2)[i]), length(data2[i,]))),
                                cbind(data1[i,], rep(paste(label@p.lab, rownames(data1)[i]), length((data1[i,]))))))
    temp[,1] <- as.numeric(as.character(temp[,1]))
    plot.data <- rbind(plot.data, temp)
  }

  plot(NULL, xlab='', ylab='',xaxs='i', yaxs='i', axes=FALSE,
       xlim=c(min(plot.data[,1]-0.2), max(plot.data[,1]) + 1),
       ylim=c(+0.5, nrow(data1)*2+0.5), type='n')

  boxplot(plot.data[,1] ~ plot.data[,ncol(plot.data)],
          horizontal=TRUE, names = c(""), show.names = FALSE,
          col = box.colors, axes = FALSE,
          outcol = c(col[1], col[2]), add = TRUE)

  mn          <- as.integer(c(min(plot.data[,1])))
  mx          <- as.integer(c(max(plot.data[,1])))
  ticks       <- mn:mx
  for (v in ticks) {
    abline(v=v, lty=3, col='lightgrey')
  }
  tick.labels <- formatC(10^ticks, format='E', digits=0)
  axis(side=1, at=ticks, labels=tick.labels, cex.axis=0.7)
  legend('topright', legend=c(label@p.lab, label@n.lab), fill=rev(col), bty='n')
  associations.labels.plot(row.names(data1), plot.type='box', verbose=verbose)
  if(verbose>2) cat("+ finished associations.box.plot\n")
}

# quantile.box plot
associations.quantile.box.plot <- function(data1, data2, label, col, verbose=1){
  if(verbose>2) cat("+ starting associations.quantile.box.plot\n")
  x.col <- col[2]
  y.col <- col[1]

  p.m = min(c(min(data1, na.rm=TRUE), min(data2, na.rm=TRUE)))
  plot(rep(p.m, dim(data1)[1]), 1:dim(data1)[1],
       xlab='', ylab='', yaxs='i', axes=FALSE,
       xlim=c(p.m, 0), ylim=c(0.5, dim(data1)[1]+0.5), frame.plot=FALSE, type='n')
  for (v in seq(p.m,-1,1)) {
    abline(v=v, lty=3, col='lightgrey')
  }

  tck = floor(p.m):0
  axis(1, tck, formatC(10^tck, format='E', digits=0), las=1, cex.axis=0.7)

  x.q = apply(data1, 1, function (x) quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95),
                                              na.rm=TRUE, names=FALSE))
  y.q = apply(data2, 1, function (x) quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95),
                                              na.rm=TRUE, names=FALSE))

  # inter-quartile range
  rect(x.q[2,], 1:dim(data1)[1], x.q[4,], (1:dim(data1)[1])+0.3, col=x.col)
  rect(y.q[2,], (1:dim(data2)[1])-0.3, y.q[4,], 1:dim(data2)[1], col=y.col)

  # 90% interval
  segments(x.q[1,], 1:dim(data1)[1], x.q[5,], 1:dim(data1)[1])#, col=x.col)
  segments(y.q[1,], 1:dim(data1)[1], y.q[5,], 1:dim(data1)[1])#, col=x.col)
  segments(x.q[1,], y0=1:dim(data1)[1], y1=(1:dim(data1)[1])+0.2)
  segments(y.q[1,], y0=(1:dim(data1)[1])-0.2, y1=1:dim(data1)[1])
  segments(x.q[5,], y0=1:dim(data1)[1], y1=(1:dim(data1)[1])+0.2)
  segments(y.q[5,], y0=(1:dim(data1)[1])-0.2, y1=1:dim(data1)[1])

  # median
  segments(x.q[3,], y0=1:dim(data1)[1], y1=(1:dim(data1)[1])+0.3, lwd=3)#, col=x.col)
  segments(y.q[3,], y0=(1:dim(data1)[1])-0.3, y1=1:dim(data1)[1], lwd=3)#, col=y.col)

  # scatter plot on top
  for (i in 1:dim(data1)[1]) {
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

    points(data1[i,], rep(i+0.15, ncol(data1))+rnorm(ncol(data1),sd=0.03), pch=16, cex=0.6, col=x.col)
    points(data2[i,], rep(i-0.15, ncol(data2))+rnorm(ncol(data2),sd=0.03), pch=16, cex=0.6, col=y.col)
  }
  legend('topright', legend=c(label@p.lab, label@n.lab), fill=rev(col), bty='n')
  associations.labels.plot(row.names(data1), plot.type='quantile.box',verbose=verbose)
  if(verbose>2) cat("+ finished associations.quantile.box.plot\n")
}

# quantile.rect plot
associations.quantile.rect.plot <- function(data1, data2, label, col, verbose=1){
  if(verbose>2) cat("+ starting associations.quantile.rect.plot\n")
  quantiles.vector <- c(0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9)

  x.q = apply(data1, 1, function (x) quantile(x, quantiles.vector, na.rm=TRUE, names=FALSE))
  x.medians = apply(data1, 1, function (x) median(x))

  y.q = apply(data2, 1, function (x) quantile(x, quantiles.vector, na.rm=TRUE, names=FALSE))
  y.medians = apply(data2, 1, function (x) median(x))

  p.m = min(c(min(data1, na.rm=TRUE), min(data2, na.rm=TRUE)))

  plot(rep(p.m, dim(data1)[1]), 1:dim(data1)[1],
       xlab='', ylab='', yaxs='i', axes=FALSE,
       xlim=c(min(data1,data2), max(data1,data2+2)),
       ylim=c(0, dim(data1)[1]), frame.plot=FALSE, type='n')
  for (v in seq(p.m,0,1)) {
    abline(v=v, lty=3, col='lightgrey')
  }

  tck = floor(p.m):0
  axis(1, tck, formatC(10^tck, format='E', digits=0), las=1, cex.axis=0.7)
  # create different tints of the colours
  colors.p <- rev(sapply(seq(0,1, length.out = 4), FUN=function(x){rgb(matrix(col2rgb(col[2])/255 + (1 - col2rgb(col[2])/255)*x, ncol=3))}))
  colors.n <- rev(sapply(seq(0,1, length.out = 4), FUN=function(x){rgb(matrix(col2rgb(col[1])/255 + (1 - col2rgb(col[1])/255)*x, ncol=3))}))
  for (i in 1:(nrow(x.q)/2)){
    if (i == 1) {
      rect(x.q[i,], 0.5:dim(data1)[1], x.q[nrow(x.q)+1-i,], (0.5:dim(data1)[1])+0.3, col = c("white"), border = c("black"), lwd = 0.9)
      rect(y.q[i,], 0.5:dim(data2)[1], y.q[nrow(y.q)+1-i,], (0.5:dim(data2)[1])-0.3, col = c("white"), border = c("black"), lwd = 0.9)
    } else {
      rect(x.q[i,], 0.5:dim(data1)[1], x.q[nrow(x.q)+1-i,], (0.5:dim(data1)[1])+0.3, col = colors.p[i], border = c("black"), lwd = 0.9)
      rect(y.q[i,], 0.5:dim(data2)[1], y.q[nrow(y.q)+1-i,], (0.5:dim(data2)[1])-0.3, col = colors.n[i], border = c("black"), lwd = 0.9)
    }
  }

  points(x.medians, y=(0.5:dim(data1)[1])+0.15, pch=18, cex = min(35/nrow(data1),4))
  points(y.medians, y=(0.5:dim(data2)[1])-0.15, pch=18, cex = min(35/nrow(data2),4))

  mtext('Quantiles', 3, line=0, at=1, adj = 1.675, padj = 0.45, las=1, cex=0.7)
  legend(-1.75, nrow(data1), legend = c("40%-60%","30%-70%", "20%-80%","10%-90%","median","","","","",""), bty='n', cex=1,
         fill=c(rev(colors.p), 'white', rev(colors.n), 'white'), lwd <- 1.3, ncol = 2,
         border = c("black", "black", "black", "black", "white","black","black","black","black", "white"))
  legend(-1.675, nrow(data1), legend = c("","","","",""),
         bty='n', lty = c(0,0,0,0,0),
         # cap legend size for diamond (should look symmetric to other symbols)
         pch = 18, cex = 1, pt.cex = c(0,0,0,0, min(35/nrow(data1), 2.25)))
  associations.labels.plot(row.names(data1), plot.type='quantile.rect',verbose=verbose)
  if(verbose>2) cat("+ finished associations.quantile.rect.plot\n")
}

### Prepare margins for the first plots
#     make left margin as big as the longest label or maximally 20.1 lines
associations.margins.plot <- function(species_names, p.label, verbose=1){
  if(verbose>2) cat("+ starting associations.margins.plot\n")
  cex.org <- par()$cex
  par(mar=c(5.1, 18, 4.1, 1.1), cex=1)
  temp = par()$mai
  cex.labels <- min(.7,(((par()$pin[2]/length(species_names))*.6)/max(strheight(species_names, units = 'inches'))))
  max_name <- max(strwidth(species_names, units = 'inches', cex=cex.labels)) + temp[4]
  temp[2] <- min(temp[2], max_name)
  par(mai=temp, cex=cex.org)
  if(verbose>2) cat("+ finished associations.margins.plot\n")
}

### Plot single feature AUCs in single panel
associations.aucs.plot <- function(aucs, binary.cols,verbose=1){
  if(verbose>2) cat("+ starting associations.aucs.plot\n")
  # set margins
  par(mar=c(5.1, 0, 4.1, 1.6))
  # plot background
  plot(NULL, xlab='', ylab='',xaxs='i', yaxs='i', axes=FALSE,
       xlim=c(0,1), ylim=c(0.5, nrow(aucs)+0.5), type='n')
  ticks       <- seq(0, 1.0, length.out=5)
  tick.labels <- formatC(ticks, digits=2)
  # plot gridlines
  for (v in ticks) {
    abline(v=v, lty=3, col='lightgrey')
  }
  # make thicker line at .5
  abline(v=.5, lty=1, col='lightgrey')
  # plot single feature aucs
  for (i in 1:nrow(aucs)) {
    segments(x0=aucs[i, 2], x1=aucs[i, 3], y0=i, col='lightgrey', lwd=1.5)
    points(aucs[i, 1], i, pch=18, col=binary.cols[i])
    points(aucs[i, 1], i, pch=5, col='black', cex=0.9)
  }

  # Title and axis label
  axis(side=1, at=ticks, labels=tick.labels, cex.axis=0.7)
  title(main='Feature AUCs', xlab='AU-ROC')
  if(verbose>2) cat("+ finished associations.aucs.plot\n")
}

### Plot fold changes in single panel
associations.fcs.plot <- function(fc.all, binary.cols,  verbose=1){
  if(verbose>2) cat("+ starting associations.fcs.plot\n")
  # margins
  par(mar=c(5.1, 0, 4.1, 1.6))
  # get minimum and maximum fcs
  mx <- max(ceiling(abs(range(fc.all, na.rm=TRUE, finite=TRUE))))
  mn    <- -mx
  # plot background
  plot(NULL, xlab='', ylab='', xaxs='i', yaxs='i', axes=FALSE,
       xlim=c(mn, mx), ylim=c(0.2, length(fc.all)+0.2), type='n')
  grid(NULL, NA, lty=3, col='lightgrey')
  # plot bars
  barplot(fc.all, horiz=TRUE, width=0.6, space=2/3,
    col=binary.cols, axes=FALSE, add=TRUE, names.arg=FALSE)
  # gridlines and axes labels
  ticks <- seq(from=mn, to=mx, length.out=5)
  tick.labels <- formatC(ticks, digits=2)
  axis(side=1, at=ticks, labels=tick.labels, cex.axis=0.7)
  title(main='Fold change', xlab='Pseudo Fold Change')
  if(verbose>2) cat("+ finished associations.fcs.plot\n")
}

### Plot prevalence shifts in single panel
associations.pr.shift.plot <- function(pr.shifts, col, verbose=1){
  if(verbose>2) cat("+ starting associations.pr.shift.plot\n")
  # margins
  par(mar=c(5.1, 0, 4.1, 1.6))

  # plot background
  plot(NULL, xlab='', ylab='', xaxs='i', yaxs='i', axes=FALSE,
       xlim=c(0, 1), ylim=c(0.2, nrow(pr.shifts)+0.2), type='n')

  # plot bars
  row.names(pr.shifts) <- NULL
  barplot(t(pr.shifts[, c(2,3)]),
          horiz=TRUE, axes=FALSE, add=TRUE, space=c(0, 4/3),
          beside=TRUE, width=.3, col=c(col[1], col[2]))
  # gridlines and axes labels
  ticks <- seq(from=0, to=1, length.out=5)
  for (v in ticks) {
    abline(v=v, lty=3, col='lightgrey')
  }
  tick.labels <- formatC(ticks*100, digits=3)
  axis(side=1, at=ticks, labels=tick.labels, cex.axis=0.7)
  title(main='Prevalence shift', xlab='Prevalence [%]')
  if(verbose>2) cat("+ finished associations.pr.shift.plot\n")
}

# p-vals
associations.pvals.plot <- function(p.vals, alpha,  verbose=1){
  if(verbose>2) cat("+ starting associations.pvals.plot\n")
  # margins
  par(mar=c(5.1, .0, 4.1, 1.6))
  p.vals.log <- -log10(p.vals)
  # get minimum and maximum
  mx    <- max(ceiling(abs(range(p.vals.log, na.rm=TRUE, finite=TRUE))))
  mn    <- 0
  p.vals.log[is.infinite(p.vals.log)] <- mx
  # plot background
  plot(NULL, xlab='', ylab='', xaxs='i', yaxs='i', axes=FALSE,
       xlim=c(mn, mx), ylim=c(0.2, length(p.vals)+0.2), type='n')
  grid(NULL, NA, lty=3, col='lightgrey')
  # plot bars
  barplot(p.vals.log, horiz=TRUE, width=0.6, space=2/3,
          col='lightgrey', axes=FALSE, add=TRUE, names.arg=FALSE)
  # gridlines and axes labels
  ticks <- seq(from=mn, to=mx)
  abline(v=-log10(alpha), lty=1, col='red')
  tick.labels <- formatC(ticks, digits=2)
  axis(side=1, at=ticks, labels=tick.labels, cex.axis=0.7)
  title(main='Adj. P Value', xlab='-log10 Adj. P-Value')
  if(verbose>2) cat("+ finished associations.pvals.plot\n")
}


# check if a string is a valid r color reprensentation
# from stackoverflow: Check if character string is a valid color representation
# https://stackoverflow.com/questions/13289009
is.color <- function(x){
  sapply(x, function(z) {tryCatch(is.matrix(col2rgb(z)), error = function(e) FALSE)})
}
### check the user-supplied color scheme for validity
#     color scheme may either be a single RColorBrewer palette or a vector of
#     the same length as the number of classes containing interpretable colors
#     as strings
check.color.scheme <- function(color.scheme, label, meta.studies=NULL,  verbose=1){
  if(verbose>2) cat("+ starting check.color.scheme\n")
  n.classes = ifelse(label@info$type == 'BINARY', 2, length(unique(label@label)))

  if (length(color.scheme) == 1 && class(color.scheme) == 'character'){
    if (n.classes == 2){
      # if color scheme and binary label, make colors as before
      if (!color.scheme %in% row.names(brewer.pal.info)){
        warning("Not a valid RColorBrewer palette name, defaulting to RdBu.\nSee brewer.pal.info for more information about RColorBrewer palettes.")
        color.scheme <- 'RdYlBu'
      }
      colors <- rev(colorRampPalette(brewer.pal(brewer.pal.info[color.scheme,'maxcolors'], color.scheme))(2))
    } else {
      # if color scheme and multiclass label, make colors either directly out of the palette (if n.classes smaller than maxcolors) or like before
      if (!color.scheme %in% row.names(brewer.pal.info)){
        warning("Not a valid RColorBrewer palette name, defaulting to Set3.\n  See brewer.pal.info for more information about RColorBrewer palettes.")
        color.scheme <- 'Set3'
      }
      # if color scheme and multiclass label, check that the palette is not divergent or sequential, but qualitative. Only issue warning.
      if (brewer.pal.info[color.scheme,'category'] != 'qual'){warning("Using a divergent or sequential color palette for multiclass data.")}
      if (n.classes <= brewer.pal.info[color.scheme, 'maxcolors']){
        colors <- brewer.pal(n.classes, color.scheme)
      } else {
        warning("The data contains more classes than the color.palette provides.")
        colors <- rev(colorRampPalette(brewer.pal(brewer.pal.info[color.scheme,'maxcolors'], color.scheme))(n.classes))
      }
    }
  } else if (length(color.scheme == n.classes) && all(is.color(color.scheme))){
    # if colors, check that all strings are real colors and check that the same length as n classes
    # convert color names to hex representation
    colors <- sapply(color.scheme, FUN=function(x){rgb(t(col2rgb(x)), maxColorValue = 255)}, USE.NAMES = FALSE)
  } else {
    stop("Supplied colors do not match the number of classes or are no valid colors")
  }
  # add transparency
  colors <- sapply(colors, FUN=function(x){paste0(x, '85')}, USE.NAMES = FALSE)
  if(verbose>2) cat("+ finished check.color.scheme\n")
  return(colors)
}

associations.labels.plot <- function(labels, plot.type,  verbose=1){
  if(verbose>2) cat("+ starting associations.labels.plot\n")
  adj <- rep(0, length(labels))
  if (plot.type == 'quantile.rect') adj <- rep(-0.5, length(labels))
  if (plot.type == 'box') adj <- -0.5 + 1:length(labels)
  cex.org <- par()$cex
  par(cex=1)
  cex.labels <- min(.7,(((par()$pin[2]/length(labels))*.6)/max(strheight(labels, units = 'inches'))))
  for (i in 1:length(labels)){
    mtext(labels[i], 2, line=0, at=i+adj[i], las=1, cex=cex.labels)
  }
  par(cex=cex.org)
  if(verbose>2) cat("+ finished associations.labels.plot\n")
}


### maker analysis for two-class data
#   calculate p-value with Wilcoxon
#   fold change as normalized absolute difference between quantiles
#   prevalence shift
#   single marker AUC
analyse.binary.marker<- function(feat, label, detect.lim, colors,
                                   pr.cutoff, mult.corr, alpha,
                                   max.show, sort.by, probs.fc=seq(.1, .9, .05),  verbose=1){
  if(verbose>1) cat("+ starting analyse.binary.marker\n")
  s.time <- proc.time()[3]
  ##############################################################################
  ### Calculate wilcoxon, pseudo-FC, prevalence shift, and AUC for each feature
  if(verbose>1) cat('+++ calculating effect size for each feature.\n')
  if (is.null(detect.lim)){
    warning("Pseudo-count before log-transformation not supplied! Estimating it as 5% percentile.\n")
    detect.lim <- quantile(feat[feat!=0], 0.05)
  }
  if(verbose) pb = txtProgressBar(max=nrow(feat), style=3)
  effect.size <- t(apply(feat, 1, FUN=function(x){
    # pseudo-fold change as differential quantile area
    q.p <- quantile(log10(x[label@p.idx]+detect.lim), probs=probs.fc)
    q.n <- quantile(log10(x[label@n.idx]+detect.lim), probs=probs.fc)
    fc <- sum(q.p - q.n)/length(q.p)

    # wilcoxon
    p.val <- wilcox.test(x[label@n.idx], x[label@p.idx], exact = FALSE)$p.value

    # AU-ROC
    temp  <- roc(predictor=x, response=label@label, ci=TRUE, direction='<')
    aucs <- c(temp$ci)

    # prevalence shift
    temp.n <- sum(x[label@n.idx] >= pr.cutoff)/sum(label@n.idx)
    temp.p <- sum(x[label@p.idx] >= pr.cutoff)/sum(label@p.idx)
    pr.shift <- c(temp.p-temp.n, temp.n, temp.p)
    if(verbose) setTxtProgressBar(pb, (pb$getVal()+1))
    return(c('fc' = fc, 'p.val' = p.val, 'auc' = aucs[2], 'auc.ci.l' = aucs[1], 'auc.ci.h' = aucs[3],
             'pr.shift' = pr.shift[1], 'pr.n'=pr.shift[2], 'pr.p'=pr.shift[3]))
  }))
  cat('\n')
  bcol <- ifelse(effect.size[,'auc'] >= 0.5, colors[2], colors[1])

  ### Apply multi-hypothesis testing correction
  if(!tolower(mult.corr) %in% c('none','bonferroni','holm','fdr','bhy')) {
    stop("! Unknown multiple testing correction method:', mult.corr,' Stopping!\n  Must of one of c('none','bonferroni', 'holm','fdr','bhy')")
  }
  if (mult.corr == 'none') {
    warning('WARNING: No multiple hypothesis testing performed.')
    p.adj <- effect.size[,'p.val']
  } else {
    p.adj <- p.adjust(effect.size[,'p.val'], method=tolower(mult.corr))
  }

  if(verbose>1) cat('+++ found', sum(p.adj < alpha, na.rm=TRUE), 'significant associations at a significance level <', alpha, '\n')
  idx <- which(p.adj < alpha)

  if (length(idx) == 0){
    stop('No significant associations found. Stopping.\n')
  }

  idx <- idx[order(p.adj[idx], decreasing=TRUE)]

  # # truncated the list for the following plots
  truncated = FALSE
  if (length(idx) >= max.show) {
    truncated = TRUE
    idx <- idx[(length(idx)-max.show+1):length(idx)]
  if(verbose>1) cat('+++ truncating the list of significant associations to the top', max.show, '\n')
  }

  ### Sort features
  if(verbose>2) cat('+++ sorting features\n')
  if (sort.by == 'fc') {
    fc.sign <- ifelse(effect.size[idx,'fc'] == 0, 1, sign(effect.size[idx,'fc']))
    p.adj.log <- -log10(p.adj[idx])
    p.adj.log[fc.sign == -1] = -p.adj.log[fc.sign == -1]
    idx <- idx[order(p.adj.log, decreasing=FALSE)]
  } else if (sort.by == 'p.val') {
    idx <- idx[order(p.adj[idx], decreasing=TRUE)]
  } else if (sort.by == 'pr.shift') {
    pr.sign <- ifelse(effect.size[idx,'pr.shift'] == 0, 1, sign(effect.size[idx,'pr.shift']))
    p.adj.log <- -log10(p.adj[idx])
    p.adj.log[pr.sign == -1] = -p.adj.log[pr.sign == -1]
    idx <- idx[order(p.adj.log, decreasing=FALSE)]
  } else {
    if(verbose>1) cat('+++ Unknown sorting option:', sort.by, '. Instead order by fold change.\n')
    fc.sign <- ifelse(effect.size[idx,'fc'] == 0, 1, sign(effect.size[idx,'fc']))
    p.adj.log <- -log10(p.adj[idx])
    p.adj.log[fc.sign == -1] = -p.adj.log[fc.sign == -1]
    idx <- idx[order(p.adj.log, decreasing=FALSE)]
  }
  e.time <- proc.time()[3]
  if(verbose>1) cat("+ finished analyse.binary.markerin",e.time-s.time,"s\n")
  return(list("p.val" = effect.size[idx,'p.val'],
              "fc"=effect.size[idx,'fc'],
              "aucs"=effect.size[idx,c('auc', 'auc.ci.l', 'auc.ci.h')],
              "pr.shift"=effect.size[idx,c('pr.shift', 'pr.n', 'pr.p')],
              "bcol"=bcol[idx],
              "p.adj"=p.adj[idx],
              "feat.red"=feat[idx,],
              "truncated"=truncated,
              "detect.lim"=detect.lim))
}
