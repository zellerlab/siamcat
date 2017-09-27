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
#' @param label label object
#' @param num.folds number of cross-validation folds (needs to be \code{>=2})
#' @param num.resample resampling rounds (values \code{<= 1} deactivate resampling)
#' @param stratify boolean, should the splits be stratified s. t. an equal proportion of classes are present in each fold?
#' @param inseparable defaults to NULL (obsolete)
#' @param meta metadata object (obsolete)
#' @keywords SIAMCAT data.splitter
#' @export
#' @return list containing the indices of the training and test folds and the parameters of the splits: \itemize{
#'  \item \code{$training.folds} nested list, containing for \code{length(num.folds)} the sample names of the \code{length(num.resample)} training folds;
#'  \item \code{$test.folds} nested list, containing for \code{length(num.folds)} the sample names of the \code{length(num.resample)} test folds;
#'  \item \code{$num.resample} = number of repeated samplings;
#'  \item \code{$num.folds} = number of folds
#'}
data.splitter <- function(label, num.folds, num.resample, stratify, inseparable, meta){
  ### read label and meta-data
  # (assuming the label file has 1 column)
  # TODO delete inseparable parameter (or at least add an default)
  # TODO delete meta parameter

  if (is.null(inseparable) || inseparable=='' || toupper(inseparable)=='NULL' || toupper(inseparable)=='NONE' || toupper(inseparable)=='UNKNOWN') {
    inseparable <- NULL
    cat('+++ Inseparable parameter not specified\n')
  }
  labelNum        <- as.numeric(label$label)
  names(labelNum) <- names(label$label)
  exm.ids         <- names(labelNum)

  # parse label description
  classes      <- sort(label$info$class.descr)


  ### check arguments
  if (num.resample < 1) {
    cat('+++ Resetting num.resample = 1 (', num.resample, ' is an invalid number of resampling rounds)\n', sep='')
    num.resample   <- 1
  }
  if (num.folds < 2) {
    cat('+++ Resetting num.folds = 2 (', num.folds, ' is an invalid number of folds)\n', sep='')
    num.folds      <- 2
  }
  if (!is.null(inseparable) && stratify) {
    cat('+++ Stratification is not supported when inseparable is given\n')
    stratify       <- FALSE
  }
  if (num.folds >= length(labelNum)) {
    cat('+++ Performing leave-one-out (LOO) cross-validation\n')
    if (stratify) {
      cat('   WARNING: Stratification is not possible with LOO cross-validation\n')
      stratify     <- FALSE
    }
    num.folds   <- length(labelNum)
  }


  train.list <- list(NULL)
  test.list  <- list(NULL)


  for (r in 1:num.resample) {
    labelNum      <- sample(labelNum)
    foldid        <- rep(0, length(labelNum))
    foldid        <- assign.fold(label = labelNum, num.folds, stratified = stratify, inseparable = NULL, foldid = foldid)
    names(foldid) <- names(labelNum)
    stopifnot(length(labelNum) == length(foldid))
    stopifnot(length(unique(foldid)) == num.folds)

    train.temp    <- list(NULL)
    test.temp     <- list(NULL)

    cat('\n+++ Splitting the dataset:\n')
    for (f in 1:num.folds) {
      # make sure each fold contains examples from all classes
      stopifnot(all(sort(unique(labelNum[foldid==f])) == classes))

      # select test examples
      test.idx        <- which(foldid == f)
      train.idx       <- which(foldid != f)
      train.temp[f] <- list(names(foldid)[train.idx])
      test.temp[f]  <- list(names(foldid)[test.idx])
      stopifnot(length(intersect(train.idx, test.idx)) == 0)
      cat('   + Fold ', f, ' contains ', sum(foldid==f), ' examples\n', sep='')
    }
    train.list[[r]] <- train.temp
    test.list[[r]]  <- test.temp
  }


  return(list("training.folds" = train.list, "test.folds" = test.list, "num.resample" = num.resample, "num.folds" = num.folds))
}
