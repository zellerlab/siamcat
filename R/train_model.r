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
#' @description This function trains the a machine learning model on the training data, using a \code{num.folds}-fold internal cross-validation scheme to find the optimal hyper-parameters of the model.
#' @param feat features object
#' @param label label object
#' @param cl class of learner, directly passed to \link[mlr]{makeLearner}
#' @param data.split filename containing the training samples or list of training instances produced by \link{data.splitter()}, defaults to \code{NULL} leading to training on the complete dataset
#' @param stratify boolean, should the folds in the internal cross-validation be stratified?
#' @param min.nonzero.coeff integer number of minimum nonzero coefficients that should be present in the model
#' @export
#' @keywords SIAMCAT plm.trainer
#' @return an object of class \link[mlr]{makeWrappedModel}
# TODO add details section for this function
train.model2 <- function(feat, label,  method = c("lasso", "enet", "ridge", "lasso_ll", "ridge_ll", "randomForest"),
                        data.split=NULL, stratify = TRUE,
                        modsel.crit=list("auc"),  min.nonzero.coeff = 1){
  # TODO 1: modsel.criterion should be implemented
  # check modsel.crit
  if (!all(modsel.crit %in% c("auc", "f1", "acc", "pr"))){
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
    } else if (m == 'pr'){
      measure[[length(measure)+1]] <- auprc
    }
  }
  # TODO 2: instead of filename containing the traning sample indices, provide the list from data.splitter

  # transpose feature matrix as a convenience preprocessing
  feat         <- t(feat)
  ### subselect training examples as specified in fn.train.sample (if given)
  foldList     <- get.foldList(data.split)
  fold.name    <- foldList$fold.name
  fold.exm.idx <- foldList$fold.exm.idx
  num.runs     <- foldList$num.runs
  num.folds    <- foldList$num.folds

  cat('\nPreparing to train', method,  'models on', num.runs, 'training set samples...\n\n')

  ### train one model per training sample (i.e. CV fold)
  # feat has structure: examples in rows; features in columns!
  W.mat           <- matrix(data=NA, nrow=ncol(feat), ncol=num.runs)
  rownames(W.mat) <- c(colnames(feat))
  colnames(W.mat) <- paste('M', fold.name, sep='_')

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
    model             <- train.plm(data=data, method = method, measure=measure, min.nonzero.coeff=min.nonzero.coeff)
    if(!all(model$feat.weights == 0)){
       models.list[[r]]  <- model
    }else{
      stop("OOOPS!!\n")
    }
    stopifnot(all(names(model$W) == rownames(W.mat)))
    W.mat[,r]          <- as.numeric(c(model$feat.weights))
    cat('\n')
  }
  # Preprocess hyper parameters
  ### Write models into matrix to reload in plm_predictor.r
  for (i in 1:length(models.list)){
    if(method %in% c("lasso", "enet", "ridge")){
      beta <- models.list[[i]]$learner.model$glmnet.fit$beta
      nRowVec <- 1
      if(!all(is.null(dim(beta))))  nRowVec <- nrow(beta)
      vec <- rep(NA, nRowVec + 2)
      vec[1] <- 0
      vec[2] <- beta
      vec[3:length(vec)] <- as.numeric(beta)
      if (i == 1) {
        out.matrix <- matrix(vec)
      } else {
        out.matrix <- cbind(out.matrix, vec)
      }
      # This overwrites rownames everytime, but doesnt need an additional conditional statement.
      # paste0 pastes two equal-length string vectors element-wise
      rownames(out.matrix) <- c("lambda", "a0", rownames(beta))
      # In the case of glmnet, make the coefficient matrix from a sparse matrix into a regular one.
      #models.list[[i]]$original.model$beta <- as.matrix(models.list[[i]]$learner.model$glmnet.fit$beta)
    } else if(method == "lasso_ll" || method == "ridge_ll"){
            # Liblinear needs C, W (intercept term is included in W).
      # Furthermore, it needs an element called "ClassNames" which is important in determining which class label is positive or negative.
      vec <- rep(NA, length(models.list[[i]]$learner.model$W) + 3)
      vec[1] <- 0
      vec[2:3] <- as.numeric(models.list[[i]]$learner.model$ClassNames)
      vec[4:length(vec)] <- as.numeric(models.list[[i]]$learner.model$W)
      if (i == 1) {
        out.matrix <- matrix(vec)
      } else {
        out.matrix <- cbind(out.matrix, vec)
      }
      # This overwrites rownames everytime, but doesnt need an additional conditional statement.
      # paste0 pastes two equal-length string vectors element-wise
      # Note that the weight-vector is a row-vector here.
      rownames(out.matrix) <- c("C", "negative label", "positive label", colnames(models.list[[i]]$learner.model$W))

    }else if(method == "randomForest"){
      vec <- rep(NA, length(models.list[[i]]$learner.model$importance) + 3)
      vec[1] <- 0
      vec[2:3] <- as.numeric(models.list[[i]]$learner.model$classes)
      vec[4:length(vec)] <- as.numeric(models.list[[i]]$learner.model$importance)
      if (i == 1) {
        out.matrix <- matrix(vec)
      } else {
        out.matrix <- cbind(out.matrix, vec)
      }
      # This overwrites rownames everytime, but doesnt need an additional conditional statement.
      # paste0 pastes two equal-length string vectors element-wise
      # Note that the weight-vector is a row-vector here.
      rownames(out.matrix) <- c("C", "negative label", "positive label", rownames(models.list[[i]]$learner.model$importance))
    }
  }
  colnames(out.matrix) = paste('M', fold.name, sep='_')
  #save(power,file="power.RData")
  invisible(list(out.matrix=out.matrix, W.mat=W.mat, models.list=models.list))
}

auprc <- makeMeasure(id = "auprc", minimize = FALSE, best = 1, worst = 0,
  properties = c("classif", "req.pred", "req.truth", "req.prob"),
  name = "Area under the Precision Recall Curve",
  fun = function(task, model, pred, feats, extra.args) {
    #if (anyMissing(pred$data$response) || length(unique(pred$data$truth)) == 1L)
    #  return(NA_real_)
    measureAUPRC(getPredictionProbabilities(pred), pred$data$truth, pred$task.desc$negative, pred$task.desc$positive)
  }
)

measureAUPRC <- function(probs, truth, negative, positive){
  pr <- pr.curve(scores.class0 = probs[which(truth == positive)],
                 scores.class1 = probs[which(truth == negative)])
  return(pr$auc.integral)
}
