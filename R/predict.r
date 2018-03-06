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

#' @title Prediction on the test set
#' @description This function takes the test set instances and the model trained by \link{plm.trainer} in order to predict the classes.
#' @param feat features object
#' @param label label object
#' @param data.split filename containing the test samples or list of test instances produced by \link{data.splitter()}, defaults to \code{NULL} leading to testing on the complete dataset
#' @param models.list list of models trained by \link{plm.trainer}
#' @export
#' @keywords SIAMCAT plm.predictor
#' @return list containing the precitions \itemize{
#'  \item \code{$pred};
#'  \item \code{$mat}
#'}
make.predictions <- function(feat, label, data.split=NULL, models.list){

  feat         <- t(feat)

  ### subselect training examples as specified in fn.train.sample (if given)
  foldList     <- get.foldList(data.split, label, mode="test", model=models.list)
  fold.name    <- foldList$fold.name
  fold.exm.idx <- foldList$fold.exm.idx
  num.runs     <- foldList$num.runs
  num.folds    <- foldList$num.folds

  ### apply one LASSO model per test sample (i.e. CV fold)
  # predictions are made on a concatenation of test examples from all test samples
  pred = NULL
  predList = list()
  fold.pred.idx = list()

  for (r in 1:num.runs) {

    data <- as.data.frame(feat[fold.exm.idx[[r]],])
    data$label <- as.factor(1)
    model <- models.list[[r]]
    cat('Applying ', colnames(model$W)[r], ' on ', fold.name[r], ' (', r, ' of ', num.runs, ')...\n', sep='')

    # subselect test examples
    task      <- makeClassifTask(data = data, target = "label")
    pdata    <- predict(model,  task = task)

    p        <- label$negative.lab+abs(label$positive.lab-label$negative.lab)*pdata$data[,4]
    names(p) <- rownames(pdata$data)

    pred     <- c(pred, p)
    fold.pred.idx[[r]] = (length(pred)-length(p)+1):length(pred)
  }

  cat('\nTotal number of predictions made:', length(pred), '\n')


  ### reformat predictions in case models were trained in repeated cross-validation
  if (length(unique(names(pred))) < length(pred)) {
    ref.names = NULL
    if (any(substr(fold.name,1,14) == 'whole data set')) {
      r.idx = c(1:num.runs) #as.numeric(sapply(strsplit(fold.name, 'predicted by model '), '[[', 2))
      runs = sort(unique(r.idx))
      stopifnot(all(runs == 1:length(runs)))
      if (!is.null(data.split)) {
        ref.names = names(label$label)
      } else {
        ref.names = unique(names(pred))
      }
    } else {
      r.idx = as.numeric(sapply(strsplit(fold.name, 'rep'), '[[', 2))
      runs = sort(unique(r.idx))
      stopifnot(all(runs == 1:length(runs)))
      if (!is.null(data.split)) {
        ref.names = names(label$label)
      } else {
        ref.names = names(pred)[unlist(fold.pred.idx[r.idx==1])]
      }
    }
    #  cat(ref.names, '\n\n')
    #  cat(names(label), '\n\n')
    #  cat(names(pred), '\n\n')

    pred.mat = matrix(data=NA, nrow=length(ref.names), ncol=length(runs))
    rownames(pred.mat) = ref.names
    if (any(substr(fold.name,1,14) == 'whole data set')) {
      colnames(pred.mat) = paste('Model', runs, sep='')
    } else {
      colnames(pred.mat) = paste('CV_rep', runs, sep='')
    }

    for (r in runs) {
      idx = which(r.idx == r)
      p = unlist(fold.pred.idx[idx])
      m = match(names(pred)[p], ref.names)
      #    cat(sort(m), '\n\n')
      #    cat(length(m), '\n\n')
      #    cat(length(label), '\n\n')
      #if (!is.null(data.split)) {
      #   stopifnot(all(sort(m) == 1:length(label$label)))
      #}
      pred.mat[m,r] = pred[p]
      stopifnot(all(names(pred)[p] == rownames(pred.mat)[m]))
    }
    correlation <- cor(pred.mat, method='spearman')
    cat('\nCorrelation between predictions from repeated CV:\n')
    cat('Min: ', min(correlation), ', Median: ', median(correlation), ', Mean: ', mean(correlation), '\n', sep='')
  } else {
    pred.mat = as.matrix(pred,byrow=TRUE)
  }

  invisible(list(pred = pred, mat = pred.mat))
}
