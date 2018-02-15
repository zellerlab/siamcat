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

#' @title Split a dataset into training and a test sets.
#' @description This function prepares the cross-validation by splitting the data into \code{num.folds} training and test folds for \code{num.resample} times.
#' @param siamcat object of class \link{siamcat-class}
#' @param num.folds number of cross-validation folds (needs to be \code{>=2}), defaults to \code{2}
#' @param num.resample resampling rounds (values \code{<= 1} deactivate resampling), defaults to \code{1}
#' @param stratify boolean, should the splits be stratified s. t. an equal proportion of classes are present in each fold?, defaults to \code{TRUE}
#' @param inseparable column index or column name of metadata variable, defaults to \code{NULL}
#' @keywords SIAMCAT data.splitter
#' @export
#' @return object of class \link{siamcat-class}
# TODO add detail section for this function
data.splitter <- function(siamcat, num.folds=2, num.resample=1, stratify=TRUE, inseparable=NULL){
  ### read label and meta-data
  # (assuming the label file has 1 column)
  if (is.null(inseparable) || inseparable=='' || toupper(inseparable)=='NULL' || toupper(inseparable)=='NONE' || toupper(inseparable)=='UNKNOWN') {
    inseparable <- NULL
  #   cat('+++ Inseparable parameter not specified\n')
  }
  labelNum        <- as.numeric(siamcat@label@label)
  names(labelNum) <- names(siamcat@label@label)
  exm.ids         <- names(labelNum)

  # parse label description
  classes      <- sort(c(siamcat@label@negative.lab,siamcat@label@positive.lab))

  ### check arguments
  if (num.resample < 1) {
    cat('+++ Resetting num.resample = 1 (', num.resample, ' is an invalid number of resampling rounds)\n', sep='')
    num.resample  <- 1
  }
  if (num.folds < 2) {
    cat('+++ Resetting num.folds = 2 (', num.folds, ' is an invalid number of folds)\n', sep='')
    num.folds     <- 2
  }
  if (!is.null(inseparable) && stratify) {
    cat('+++ Resetting stratify to FALSE (Stratification is not supported when inseparable is given)\n')
    stratify      <- FALSE
  }
  if (num.folds >= length(labelNum)) {
    cat('+++ Performing un-stratified leave-one-out (LOO) cross-validation\n')
    stratify      <- FALSE
    num.folds     <- length(labelNum)-1
  }
  if (!is.null(inseparable) && is.null(siamcat@phyloseq@sam_data)){
    stop('Meta-data must be provided if the inseparable parameter is not NULL')
  }
  if (!is.null(inseparable)){
    if (is.numeric(inseparable) && length(inseparable) == 1){
      stopifnot(inseparable <= ncol(siamcat@phyloseq@sam_data))
    } else if (class(inseparable) == 'character' && length(inseparable == 1)){
      stopifnot(inseparable %in% colnames(siamcat@phyloseq@sam_data))
    } else {
      stop('Inseparable parameter must be either a single column index or a single column name of metadata matrix')
    }
  }

  train.list <- list(NULL)
  test.list  <- list(NULL)


  for (r in 1:num.resample) {
    labelNum      <- sample(labelNum)
    foldid        <- assign.fold(label = labelNum, num.folds=num.folds, stratified = stratify, 
                                 inseparable = inseparable, meta=siamcat@phyloseq@sam_data)
    names(foldid) <- names(labelNum)
    stopifnot(length(labelNum) == length(foldid))
    stopifnot(length(unique(foldid)) == num.folds)

    train.temp    <- list(NULL)
    test.temp     <- list(NULL)

    if (verbose) cat('\n+++ Splitting the dataset:\n')
    for (f in 1:num.folds) {
      # make sure each fold contains examples from all classes
      # for stratify==TRUE should be tested before assignment of test/training set
      if (stratify){
        stopifnot(all(sort(unique(labelNum[foldid==f])) == classes))
      }
      # select test examples
      test.idx        <- which(foldid == f)
      train.idx       <- which(foldid != f)
      train.temp[f] <- list(names(foldid)[train.idx])
      test.temp[f]  <- list(names(foldid)[test.idx])
      # for startify==FALSE, all classes must only be present in the training set
      # e.g. in leave-one-out CV, the test fold cannot contain all classes
      if (!stratify){
        stopifnot(all(sort(unique(labelNum[foldid != f])) == classes))
      }
      stopifnot(length(intersect(train.idx, test.idx)) == 0)
      if(verbose) cat('   + Fold ', f, ' contains ', sum(foldid==f), ' examples\n', sep='')
    }
    train.list[[r]] <- train.temp
    test.list[[r]]  <- test.temp
  }

  siamcat@dataSplit <- new("dataSplit",training.folds=train.list,
                           test.folds=test.list,
                           num.resample=num.resample,
                           num.folds=num.folds)
  return(siamcat)
}


##### function to partition training set into cross-validation folds for model selection
### Works analogous to the function used in data_splitter.r
#' @export
assign.fold <- function(label, num.folds, stratified, inseparable = NULL, meta=NULL) {
  foldid  <- rep(0, length(label))
  classes <- sort(unique(label))
  # Transform number of classes into vector of 1 to x for looping over.
  # stratify positive examples
  if (stratified) {
    # If stratify is TRUE, make sure that num.folds does not exceed the maximum number of examples for the class with the fewest training examples.
    if (any(as.data.frame(table(label))[,2] < num.folds)) {
      stop("+++ Number of CV folds is too large for this data set to maintain stratification. Reduce num.folds or turn stratification off. Exiting.\n")
      # q(status = 1)
    }
    for (c in 1:length(classes)) {
      idx         <- which(label==classes[c])
      foldid[idx] <- sample(rep(1:num.folds, length.out=length(idx)))
    }
  } else {
    # If stratify is not TRUE, make sure that num.sample is not bigger than number.folds
    if (length(label) <= num.folds){
      cat("+++ num.samples is exceeding number of folds, setting CV to (k-1) unstratified CV\n")
      num.folds   <- length(label)-1
    }
    if (!is.null(inseparable)) {
      strata      <- unique(meta[,inseparable])
      sid         <- sample(rep(1:num.folds, length.out=length(strata)))
      for (s in 1:length(strata)) {
        idx          <- which(meta[,inseparable] == strata[s])
        foldid[idx]  <-sid[s]
      }
      stopifnot(all(!is.na(foldid)))
    } else {
      foldid      <- sample(rep(1:num.folds, length.out=length(label)))
    }

  }

  # make sure that for each test fold the training fold
  # (i.e. all other folds together) contain examples from all classes
  # except for stratified CV
  if (!stratified){
    for (f in 1:num.folds) {
      stopifnot(all(sort(unique(label[foldid!=f])) == classes))
    }
  } else {
    for (f in 1:num.folds) {
      stopifnot(all(sort(unique(label[foldid==f])) == classes))
    }
  }

  stopifnot(length(label) == length(foldid))

  return(foldid)
}
