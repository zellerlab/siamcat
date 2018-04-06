#!/usr/bin/Rscript
###
# SIAMCAT -  Statistical Inference of Associations between
#   Microbial Communities And host phenoTypes
# R flavor
# EMBL Heidelberg 2012-2018
# GNU GPL 3.0
###

#' @title Model training
#' @description This function trains the a machine learning model on the
#'        training data
#' @param siamcat object of class \link{siamcat-class}
#' @param method string, specifies the type of model to be trained, may be one
#'        of these: \code{c("lasso", "enet", "ridge", "lasso_ll", "ridge_ll",
#'        "randomForest")}
#' @param stratify boolean, should the folds in the internal cross-validation be
#'        stratified?, defaults to \code{TRUE}
#' @param modsel.crit list, specifies the model selection criterion during
#'        internal cross-validation, may contain these: \code{c("auc", "f1",
#'        "acc", "pr")}, defaults to \code{list("auc")}
#' @param min.nonzero.coeff integer number of minimum nonzero coefficients that
#'        should be present in the model (only for \code{"lasso"},
#'        \code{"ridge"}, and \code{"enet"}, defaults to \code{1}
#' @param param.set a list of extra parameters for mlr run, may contain:
#'        \itemize{
#'          \item \code{cost} - for lasso_ll and ridge_ll
#'          \item \code{alpha} for enet
#'          \item \code{ntree} and \code{mtry} for RandomForrest.
#'        } Defaults to \code{NULL}
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'        for only information about progress and success, \code{2} for normal
#'        level of information and \code{3} for full debug information,
#'        defaults to \code{1}
#' @export
#' @keywords SIAMCAT plm.trainer
#' @return object of class \link{siamcat-class} with added \code{model_list}
#' @details This functions performs the training of the machine learning model
#'        and functions as an interface to the \code{mlr}-package.
#'
#'        The function expects a \link{siamcat-class}-object with a prepared
#'        cross-validation (see \link{create.data.split}) in the
#'        \code{data_split}-slot of the object. It then trains a model for
#'        each fold of the datasplit.
#'
#'        For the machine learning methods that require additional
#'        hyperparameters (e.g. \code{lasso_ll}), the optimal hyperparameters
#'        are tuned with the function \link[mlr]{tuneParams} within the
#'        \code{mlr}-package.
#'
#'        The methods \code{"lasso"}, \code{"enet"}, and \code{"ridge"} are
#'        implemented as mlr-taks using the \code{"classif.cvglmnet"} Learner,
#'        \code{"lasso_ll"} and \code{"ridge_ll"} use the
#'        \code{"classif.LiblineaRL1LogReg"} and the
#'        \code{"classif.LiblineaRL2LogReg"} Learners respectively. The
#'        \code{"randomForest"} method is implemented via the
#'        \code{"classif.randomForest"} Learner.
#' @examples
#'
#'  data(siamcat_example)
#'  # simple working example
#'  siamcat_validated <- train.model(siamcat_example, method='lasso')
#'
train.model <- function(siamcat,  method = c("lasso", "enet", "ridge", "lasso_ll", "ridge_ll", "randomForest"),
                        stratify = TRUE, modsel.crit=list("auc"),  min.nonzero.coeff = 1, param.set=NULL, verbose=1){

  if(verbose>1) cat("+ starting train.model\n")
  s.time <- proc.time()[3]
  # check modsel.crit
  if (!all(modsel.crit %in% c("auc", "f1", "acc", "pr", "auprc"))){
    warning("Unkown model selection criterion... Defaulting to AU-ROC!\n")
    measure <- list(mlr::auc)
  } else {
    measure <- list()
  }
  if(verbose>2) cat("+++ preparing selection measures\n")
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
                           measureAUPRC(getPredictionProbabilities(pred), pred$data$truth, pred$task.desc$negative, pred$task.desc$positive)
                           })
      measure[[length(measure)+1]] <- auprc
    }
  }

  # Create matrix with hyper parameters.
  hyperpar.list   <- list()

  # Create List to save models.
  models.list     <- list()
  power           <- NULL
  num.runs        <- siamcat@data_split@num.folds * siamcat@data_split@num.resample
  bar             <- 0
  if(verbose>1) cat('+ training', method,  'models on', num.runs, 'training sets\n')

  if(verbose==1 || verbose==2) pb <- txtProgressBar(max=num.runs, style=3)

  for (fold in 1:siamcat@data_split@num.folds) {

    if(verbose>2) cat('+++ training on cv fold:', fold, '\n')

    for (resampling in 1:siamcat@data_split@num.resample) {

      if(verbose>2) cat('++++ repetition:', resampling, '\n')

      fold.name    <- paste0('cv_fold', as.character(fold), '_rep', as.character(resampling))
      fold.exm.idx <- match(siamcat@data_split@training.folds[[resampling]][[fold]], names(siamcat@label@label))

      ### subselect examples for training
      label.fac         <- factor(siamcat@label@label, levels=c(siamcat@label@negative.lab, siamcat@label@positive.lab))
      train.label       <- label.fac[fold.exm.idx]
      data              <- as.data.frame(t(siamcat@phyloseq@otu_table)[fold.exm.idx,])
      stopifnot(nrow(data)         == length(train.label))
      stopifnot(all(rownames(data) == names(train.label)))
      data$label        <- train.label

      ### internal cross-validation for model selection
      model             <- train.plm(data=data,  method = method, measure=measure, min.nonzero.coeff=min.nonzero.coeff,
                                     param.set=param.set, neg.lab=siamcat@label@negative.lab)
      bar <- bar+1

      if(!all(model$feat.weights == 0)){
        models.list[[bar]]  <- model
      }else{
        warning("Model without any features selected!\n")
      }

      if(verbose==1 || verbose==2) setTxtProgressBar(pb, bar)
    }
  }

  siamcat@model_list <- new("model_list",models=models.list,model.type=method)
  e.time            <- proc.time()[3]

  if(verbose>1)  cat("\n+ finished train.model in",e.time-s.time,"s\n")
  if(verbose==1) cat("\nTrained",method,"models successfully.\n")

  return(siamcat)
}

measureAUPRC <- function(probs, truth, negative, positive){
  pr <- pr.curve(scores.class0 = probs[which(truth == positive)],
                 scores.class1 = probs[which(truth == negative)])
  return(pr$auc.integral)
}
