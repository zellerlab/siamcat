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
