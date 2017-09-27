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

##### function to parse the header of a label file
### label.header - string in the format: #<TYPE>:<L1>=<class1>;<L2>=<class2>[;<L3>=<class3>]
###   where <TYPE> is a string specifying the type of label variable such as
###   BINARY (for binary classification), CATEGORICAL (for multi-class classification), or CONTINUOUS (for regression)
###   <L1> is a short numeric label for the first class with description <class1> (similarly for the other classes)
## #' @export
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

# ##### function to parse the header of a model file
# ### TODO documentation
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
# #' @keywords SIAMCAT make_filter_options
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
# #' @keywords SIAMCAT make_normalizer_options
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
# #' @keywords SIAMCAT make_evaler_options
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
# #' @keywords SIAMCAT make_evaler_options
#' @export
#' @examples
#' make_frozen_normalizer_options()

make_frozen_normalizer_options <- function(){
  option_list = list(make_option('--feat_in', type='character', help='Input file containing features'),
                     make_option('--param_in', type='character', help='Input file from which normalization parameters will be loaded'))
  return(option_list)
}

#' This is a helper function to uniformly pass command line parameters for metadata addition to SIAMCAT
# #' @keywords SIAMCAT make_evaler_options
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
# #' @keywords SIAMCAT make_evaler_options
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
