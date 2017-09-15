###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# R package flavor
#
# written by Georg Zeller
# with additions by Nicolai Karcher and Konrad Zych
# EMBL Heidelberg 2012-2017
#
# version 0.2.0
# file last updated: 12.06.2017
# GNU GPL-3.0
###

##### auxiliary function to trim whitespace from string
# returns string without leading or trailing whitespace
trim = function (x) {
  gsub("^\\s+|\\s+$", "", x)
}

### label.plot.horizontal() takes as input lists of (significantly) differentially abundant bacterial features and plots their names
### on the left side of a figur, parallel to each associated plot. inner.diff.x, inner.diff.y and outer.diff are numerical values that can be
### used to tweak the position of the text lines relatively to their plot. Specifically, inner.diff.y and inner.diff.x will shift the text
###  alongside the y-axis. outer.diff on the other hand is used as a multiplication factor which changes the distance between each different
### feature example combination globally.



### xv is vector containing values to draw a barplot, z and y determine upper and lower boundary of barplot, respectively.
#' @export
draw.error.bar <- function(plot.coords, z, y){
  g <- (max(plot.coords)-min(plot.coords))/(3*length(plot.coords))
  for (i in 1:length(plot.coords)) {
    lines(c(z[i],y[i]),c(plot.coords[i], plot.coords[i]))
    lines(c(z[i],z[i]),c(plot.coords[i]+g, plot.coords[i]-g))
    lines(c(y[i],y[i]),c(plot.coords[i]+g, plot.coords[i]-g))
  }
}

#' @export
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
#' @export
create.tints.rgb <- function(color.rgb, nr.tints, tint.steps = 1/nr.tints) {
  tints <- matrix(rep(0,(3*(nr.tints))),nrow=3, ncol=nr.tints)
  for (i in 1:nr.tints){
    tints[1,i] = color.rgb[1] + (1 - color.rgb[1]) * (tint.steps*i)
    tints[2,i] = color.rgb[2] + (1 - color.rgb[2]) * (tint.steps*i)
    tints[3,i] = color.rgb[3] + (1 - color.rgb[3]) * (tint.steps*i)
  }
  return (tints)
}

##### function to parse the header of a label file
### label.header - string in the format: #<TYPE>:<L1>=<class1>;<L2>=<class2>[;<L3>=<class3>]
###   where <TYPE> is a string specifying the type of label variable such as
###   BINARY (for binary classification), CATEGORICAL (for multi-class classification), or CONTINUOUS (for regression)
###   <L1> is a short numeric label for the first class with description <class1> (similarly for the other classes)
#' @export
parse.label.header <- function(label.header) {
  s = strsplit(label.header, ':')[[1]]
  type = trim(s[1])
  if (substr(type, 1, 1) == '#')
    type = trim(substr(type, 2, nchar(type)))
  class.descr = unlist(strsplit(strsplit(trim(s[2]), ';')[[1]], '='))
  l = class.descr[seq(2,length(class.descr),2)]
  class.descr = as.numeric(class.descr[seq(1,length(class.descr)-1,2)])
  names(class.descr) = l

  label.info = list()
  label.info$type = type
  label.info$class.descr = class.descr
  return(label.info)
}

##### function to parse the header of a model file
### TODO documentation
#' @export
parse.model.header <- function(model.header) {
  s = strsplit(model.header, ':')[[1]]
  type = trim(s[1])
  if (substr(type, 1, 1) == '#')
    type = trim(substr(type, 2, nchar(type)))
  label.header = trim(paste(s[2:length(s)], collapse=':'))
  if (substr(label.header, 1, 1) == '[') {
    stopifnot(substr(label.header, nchar(label.header), nchar(label.header)) == ']')
    label.header = substr(label.header, 2, nchar(label.header)-1)
  }
  p = grep('\\(.*\\)', type)
  properties = NULL
  if (length(p) > 0) {
    stopifnot(length(p) == 1)
    stopifnot(substr(type, nchar(type), nchar(type)) == ')')
    properties = substr(type, p+1, nchar(type)-1)
    type = trim(substr(type, 1, p-1))
  }

  model.info = list()
  model.info$type = type
  model.info$properties = properties
  model.info$label.header = label.header
  return(model.info)
}


### TODO docu!
#' @export
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

#' @title Read feature file
#' @description This file reads in the tsv file with features and converts it into a matrix
#' @param fn.in.feat name of the tsv file containing features
#' @export
#' @return matrix containing features from the file
read.features <- function(fn.in.feat){
  if (is.null(fn.in.feat)) stop("Filename for features file not provided!\n")
  if(!file.exists(fn.in.feat)) stop("Feature file ", fn.in.feat, " does not exist!\n")
  feat <- read.table(file = fn.in.feat, sep = '\t', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = '')
  feat <- as.matrix(feat)
  invisible(feat)
}

#' @title Read labels file
#' @description This file reads in the tsv file with labels and converts it into a matrix
#' @param fn.in.label name of the tsv file containing labels
#' @export
#' @return list with two values: - $header contains header of the file and - $label stores matrix containing features from the file
read.labels <- function(fn.in.label,feat=NULL){
  if (is.null(fn.in.label)) stop("Filename for labels file not provided!\n")
  if(!file.exists(fn.in.label)) stop("Label file ", fn.in.label, " does not exist!\n")
  label <- read.table(file=fn.in.label, sep='\t', header=TRUE, row.names=NULL, stringsAsFactors = FALSE,
                      check.names=FALSE, quote='', comment.char = "#")
  label <- as.matrix(label)
  if (dim(label)[1] > dim(label)[2]){
    temp            <- names(label)
    names(label)    <- NULL
    label           <- rbind(temp,label)
    rownames(label) <- label[,1]
    label[,1]       <- NULL
    label           <- t(label)
  }
  namesL          <- colnames(label)
  label           <- as.numeric(label)
  names(label)    <- namesL
  if(!all(is.null(feat))){
    if(any(names(label) != colnames(feat))) stop("read.labels: Names in label object and feat object do not match!\n")
  }

  # Check general suitablity of supplied dataset
  classes <- unique(label)
  for (i in classes){
    if(sum(label==i) <= 5) stop("Data set has only",sum(label==i), "training examples of class",i," This is not enough for SIAMCAT to proceed")
    if (sum(label==i) < 10){
      cat("Data set has only",sum(label==i), "training examples of class",i," . Note that a dataset this small/skewed is not necessarily suitable for analysis in this pipe line." )
    }
  }

  #Check label header!
  con               <- file(fn.in.label, 'rt')
  header            <- readLines(con, 1)
  if (substring(header,1,1) != "#"){
    stop("Label header seems to be missing or broken.")
  }
  close(con)
  label             <- list("label" = label, "header" = header)
  label$info <- parse.label.header(label$header)
  stopifnot(label$info$type == 'BINARY')
  label$positive.lab <- max(label$info$class.descr)
  label$negative.lab <- min(label$info$class.descr)
  label$n.idx <- label$label==label$negative.lab
  label$n.lab <- gsub('[_.-]', ' ', names(label$info$class.descr)[label$info$class.descr==label$negative.lab])
  label$p.idx <- label$label==label$positive.lab
  label$p.lab <- gsub('[_.-]', ' ', names(label$info$class.descr)[label$info$class.descr==label$positive.lab])
  invisible(label)
}

#' @title Read metadata file
#' @description This file reads in the tsv file with metadata and converts it into a matrix
#' @param fn.in.meta name of the tsv file containing metadata
#' @export
#' @return matrix containing features from the file
read.meta <- function(fn.in.meta){
  if (is.null(fn.in.meta) || toupper(fn.in.meta)=='NULL' || toupper(fn.in.meta)=='NONE' || toupper(fn.in.meta)=='UNKNOWN') {
    cat("Filename for metadata file not provided, continuing without it.\n")
  }else{
    if(!file.exists(fn.in.meta)) stop("Metadata file ", fn.in.meta, " does not exist!\n")
    meta <- read.table(file=fn.in.meta, sep='\t', header=TRUE, row.names=1, check.names=FALSE, quote='')
    meta <- as.matrix(meta)
  }
  invisible(meta)
}

#' @title Append source directory path
#' @description Append / at the end of the name of the directory path if it is not there
#' @param source.dir string with path to the source directory
#' @export
#' @return string with / at the end
appendDirName <- function(source.dir){
  if (substr(source.dir, nchar(source.dir), nchar(source.dir)) != '/') {
    source.dir <- paste(source.dir, '/', sep='')
  }
  invisible(source.dir)
}



#' This is a helper function to uniformly pass command line parameters forgeneric  filtering to SIAMCAT as well as the testing framework.
#' @keywords SIAMCAT make_filter_options
#' @export
#' @examples
#' make_filter_options()

make_filter_options <- function(){
  option_list = list(
  # make_option('--pkgdir', type='character', help='Source directory of dataprep'),
  make_option('--feat_in', type='character', help='Input file containing features'),
  make_option('--method', type='character', default='abundance', help='Filtering method (one of \"abundance\", \"cum.abundance\", or \"prevalence\")'),
  make_option('--cutoff', type='double', default=0.001, help='abundance / prevalence cutoff applied for filtering'),
  make_option('--recomp_prop', type='logical', default=FALSE, help='Should relative abundances be be recomputed?'),
  make_option('--rm_unmapped', type='logical', default=TRUE, help='Should the abundance of unmapped reads be removed?')
  )
  return(option_list)
}

#' This is a helper function to uniformly pass command line parameters for generic normalizing to SIAMCAT as well as the testing framework.
#' @keywords SIAMCAT make_normalizer_options
#' @export
#' @examples
#' make_normalizater_options()

make_normalizer_options <- function(){
  option_list = list(
  # make_option('--pkgdir', type='character', help='Source directory of dataprep'),
  make_option('--feat_in', type='character', help='Input file containing features'),
  make_option('--param_out', type='character', help='Output file to which normalization parameters will be written'),
  make_option('--method', type='character', default='log.std', help='Normalization method (either \"log.std\", \"rlog2.std\" or \"log.unit\")'),
  make_option('--log_n0', type='double', default=10^-8, help='Pseudocount that is added before log-transformation'),
  make_option('--sd_min_quantile', type='double', default=0.1, help='Quantile of the distribution of standard deviation of all feature that will be added to the denominator during standardization of each feature in order to avoid underestimation (only for metod==\"log.std\")'),
  make_option('--vector_norm', type='integer', default=2, help='Vector norm to use (either 1 or 2, only for method==\"log.unit\")'),
  make_option('--norm_feature', type='logical', default=FALSE, help='Normalize by feature (only for method==\"log.unit\")?'),
  make_option('--norm_sample', type='logical', default=TRUE, help='Normalize by sample (after feature normalization, only for method==\"log.unit\")?'),
  make_option('--norm_global', type='logical', default=FALSE, help='Normalize by global rescaling (only if both norm_feature and norm_sample are FALSE and only for method==\"log.unit\")?')
  )
  return(option_list)
}



#' This is a helper function to uniformly pass command line parameters for model evaluation to SIAMCAT (as well as potentially for the testing framework)
#' @keywords SIAMCAT make_evaler_options
#' @export
#' @examples
#' make_evaler_options()

make_evaler_options <- function(){
  option_list = list(
  #make_option('--pkgdir', type='character', help='Source directory of dataprep'),
  make_option(c('-s', '--srcdir'), type='character', help='Source directory of this and other utility scripts'),
  make_option('--label', type='character', help='Input file containing labels'),
  make_option('--pred', type='character', help='Input file containing the trained classification model(s)'),
  #make_option('--plot', type='character', help='Output file for plotting'),
  make_option('--write_eval_results', type='logical', default=FALSE, help='Should calculated parameters be written into tab-delimited file? (Necessary for generation of test files')
  #make_option('--output_results', type='character', help='Output file containing evaluation results (only necessary when write_eval_results is set to TRUE'),
  #make_option('--model_type', type='character', help='model type')
  )
  return(option_list)
}

#' This is a helper function to uniformly pass command line parameters for frozen normalization to SIAMCAT
#' @keywords SIAMCAT make_evaler_options
#' @export
#' @examples
#' make_frozen_normalizer_options()

make_frozen_normalizer_options <- function(){
  option_list = list(make_option('--feat_in', type='character', help='Input file containing features'),
                     make_option('--param_in', type='character', help='Input file from which normalization parameters will be loaded'))
  return(option_list)
}

#' This is a helper function to uniformly pass command line parameters for metadata addition to SIAMCAT
#' @keywords SIAMCAT make_evaler_options
#' @export
#' @examples
#' make_meta_adder_options()

make_meta_adder_options <- function(){
  option_list = list(
    make_option('--feat_in', type='character', help='Input file containing features'),
    make_option('--metadata_in', type='character', help='Input file containing metadata'),
    make_option('--pred_names', type='character', help='names (comma-separated list) of the metavariables to be added to the feature matrix as predictors'),
    make_option('--std_meta', type='logical', default=TRUE, help='Shall added (metadata) features be standardized?'))
  return(option_list)
}

#' This is a helper function to uniformly pass command line parameters for model interpretation in SIAMCAT
#' @keywords SIAMCAT make_evaler_options
#' @export
#' @examples
#' make_model_interpretor_options()

make_model_interpretor_options <- function(){
  option_list = list(
    make_option(c('-s', '--srcdir'), type='character', help='Source directory of this and other utility scripts'),
    make_option('--label', type='character', help='Input file containing labels'),
    make_option('--feat', type='character', help='Input file containing features'),
    make_option('--origin_feat', type='character', help='Input file containing unnormalized/filtered features'),
    make_option('--meta', type='character', default='NULL', help='Input file containing metadata'),
    make_option('--model', type='character', help='Input file containing the trained classification model(s)'),
    make_option('--pred', type='character', help='Input file containing the trained classification model(s)'),
    #make_option('--plot', type='character', help='Output file for plotting'),
    make_option('--col_scheme', type='character', default='RdYlBu', help='Color scheme'),
    make_option('--heatmap_type', type='character', default='zscore', help='which metric should be used to plot feature changes in heatmap? (zscore|fc)'),
    make_option('--consens_thres', type='double', default=0.5, help='specifies the minimal ratio of models incorporating a feature to include it into heatmap'))
  return(option_list)
}
