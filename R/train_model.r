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

#' @title Model training
#' @description This function trains the a machine learning model on the training data
#' @param feat feature object
#' @param label label object
#' @param method string, specifies the type of model to be trained, may be one of these: \code{c("lasso", "enet", "ridge", "lasso_ll", "ridge_ll", "randomForest")}
#' @param data.split filename containing the training samples or list of training instances produced by \link{data.splitter()}, defaults to \code{NULL} leading to training on the complete dataset
#' @param stratify boolean, should the folds in the internal cross-validation be stratified?
#' @param modsel.crit list, specifies the model selection criterion during internal cross-validation, may contain these: \code{c("auc", "f1", "acc", "pr")}
#' @param min.nonzero.coeff integer number of minimum nonzero coefficients that should be present in the model (only for \code{"lasso"}, \code{"ridge"}, and \code{"enet"}
#' @param param.set a list of extra parameters for mlr run, may contain: \code{cost} - for lasso_ll and ridge_ll; \code{alpha} for enet and \code{ntree, mtry} for RandomForrest
#' @export
#' @keywords SIAMCAT plm.trainer
#' @return list containing \itemize{
#'  \item \code{models.list} the list of trained models for each cross-validation fold
#'  \item \code{out.matrix} ?
#'  \item \code{W.mat} ?
#'}
# TODO add details section for this function
train.model <- function(feat, label,  method = c("lasso", "enet", "ridge", "lasso_ll", "ridge_ll", "randomForest"),
                        data.split=NULL, stratify = TRUE,
                        modsel.crit=list("auc"),  min.nonzero.coeff = 1, param.set=NULL){
  # TODO 1: modsel.criterion should be implemented
  # check modsel.crit
  if (!all(modsel.crit %in% c("auc", "f1", "acc", "pr", "auprc"))){
    cat("Unkown model selection criterion... Defaulting to AU-ROC!\n")
    measure <- list(mlr::auc)
  } else {
    measure <- list()
  }
  for (m in modsel.crit){
    if (m == 'auc'){
      measure[[length(measure)+1]] <- mlr::auc
    } else if (m == 'acc'){
      measure[[length(measure)+1]] <- mlr::acc
    } else if (m == 'f1'){
      measure[[length(measure)+1]] <- mlr::f1
    } else if (m == 'pr' || m == 'auprc'){
      auprc <- makeMeasure(id = "auprc", minimize = FALSE, best = 1, worst = 0,
                           properties = c("classif", "req.pred", "req.truth", "req.prob"),
                           name = "Area under the Precision Recall Curve",
                           fun = function(task, model, pred, feats, extra.args) {
                           #if (anyMissing(pred$data$response) || length(unique(pred$data$truth)) == 1L)
                           #  return(NA_real_)
                           measureAUPRC(getPredictionProbabilities(pred), pred$data$truth, pred$task.desc$negative, pred$task.desc$positive)
                           })
      measure[[length(measure)+1]] <- auprc
    }
  }
  # TODO 2: instead of filename containing the traning sample indices, provide the list from data.splitter

  # transpose feature matrix as a convenience preprocessing
  feat         <- t(feat)
  ### subselect training examples as specified in fn.train.sample (if given)
  foldList     <- get.foldList(data.split, label, mode="train")
  fold.name    <- foldList$fold.name
  fold.exm.idx <- foldList$fold.exm.idx
  num.runs     <- foldList$num.runs
  num.folds    <- foldList$num.folds

  cat('\nPreparing to train', method,  'models on', num.runs, 'training set samples...\n\n')


  # Create matrix with hyper parameters.
  hyperpar.list   <- list()

  # Create List to save models.
  models.list     <- list()
  power           <- NULL

  for (r in 1:num.runs) {
    cat('Training on ', fold.name[r], ' (', r, ' of ', num.runs, ')\n', sep='')
    ### subselect examples for training
    label.fac         <- factor(label$label, levels=c(label$negative.lab, label$positive.lab))
    train.label       <- label.fac
    train.label       <- label.fac[fold.exm.idx[[r]]]
    data              <- as.data.frame(feat[fold.exm.idx[[r]],])
    stopifnot(nrow(data)         == length(train.label))
    stopifnot(all(rownames(data) == names(train.label)))
    data$label                     <- train.label

    ### internal cross-validation for model selection
    model             <- train.plm(data=data, method = method, measure=measure, min.nonzero.coeff=min.nonzero.coeff,param.set=param.set)
    if(!all(model$feat.weights == 0)){
       models.list[[r]]  <- model
    }else{
      warning("Model without any features selected!\n")
    }
    cat('\n')
  }
 
  models.list$model.type <- method
  invisible(models.list)
}

measureAUPRC <- function(probs, truth, negative, positive){
  pr <- pr.curve(scores.class0 = probs[which(truth == positive)],
                 scores.class1 = probs[which(truth == negative)])
  return(pr$auc.integral)
}
