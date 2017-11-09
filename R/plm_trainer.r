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
plm.trainer <- function(feat, label,  method = c("lasso", "enet", "ridge", "libLineaR", "randomForest"), 
                        data.split=NULL, stratify = TRUE, 
                        modsel.crit  = "auc",  min.nonzero.coeff = 1){
  # TODO 1: modsel.criterion should be implemented
  # TODO 2: instead of filename containing the traning sample indices, provide the list from data.splitter
  # TODO 3: add model.type as parameter
  # transpose feature matrix as a convenience preprocessing

  feat         <- t(feat)

  label.fac                  <- factor(label$label, levels=c(label$negative.lab, label$positive.lab)) 

  data                       <- cbind(feat,label.fac[rownames(feat)])
  data                       <- as.data.frame(data)
  data[,ncol(data)]          <- as.factor(data[,ncol(data)]) ### ???!!!
  colnames(data)             <- paste0("Sample_",1:ncol(data))
  colnames(data)[ncol(data)] <- "label"
  ### subselect training examples as specified in fn.train.sample (if given)
  num.runs     <- 1
  fold.name    <- list()
  fold.exm.idx <- list() 
  num.folds    <- 2

  # print(fn.train.sample)
  if (is.null(data.split)){
    # train on whole data set
    fold.name[[1]]    <- 'whole data set'
    fold.exm.idx[[1]] <- names(label$label)
  } else {
    if (class(data.split) == 'character') {
      # read in file containing the training instances
      num.runs      <- 0
      con           <- file(data.split, 'r')
      input         <- readLines(con)
      close(con)
      #print(length(input))
      
      num.folds     <- as.numeric(strsplit(input[[3]],"#")[[1]][2])
      for (i in 1:length(input)) {
        l               <- input[[i]]
        if (substr(l, 1, 1) != '#') {
          num.runs                 <- num.runs + 1
          s                        <- unlist(strsplit(l, '\t'))
          fold.name[[num.runs]]    <- substr(s[1], 2, nchar(s[1]))
          ### Note that the %in%-operation is order-dependend.
          fold.exm.idx[[num.runs]] <- which(names(label$label) %in% as.vector(s[2:length(s)]))
          cat(fold.name[[num.runs]], 'contains', length(fold.exm.idx[[num.runs]]), 'training examples\n')
          #      cat(fold.exm.idx[[num.runs]], '\n\n')
          #    } else {
          #      cat('Ignoring commented line:', l, '\n\n')
        }
      }
    } else if (class(data.split) == 'list') {
      # use training samples as specified in training.folds in the list
      num.folds <- data.split$num.folds
      num.runs <- 0
      for (cv in 1:data.split$num.folds){
        for (res in 1:data.split$num.resample){
          num.runs <- num.runs + 1

          fold.name[[num.runs]] = paste0('cv_fold', as.character(cv), '_rep', as.character(res))
          fold.exm.idx[[num.runs]] <- match(data.split$training.folds[[res]][[cv]], names(label$label))
          cat(fold.name[[num.runs]], 'contains', length(fold.exm.idx[[num.runs]]), 'training examples\n')
        }
      }
    } else {
      stop('Wrong input for training samples!...')
    }
  }
  #stop()
  fold.name     <- unlist(fold.name)
  stopifnot(length(fold.name) == num.runs)
  stopifnot(length(fold.exm.idx) == num.runs)
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
    cat('Training on ', fold.name[r], ' (', r, ' of ', num.runs, ')', sep='')
    ### subselect examples for training
    train.label       <- label
    train.label$label <- train.label$label[fold.exm.idx[[r]]]
    train.feat        <- feat[fold.exm.idx[[r]],]
    stopifnot(nrow(train.feat)         == length(train.label$label))
    stopifnot(all(rownames(train.feat) == names(train.label$label)))

    # reorder training examples so that class order is the same for all models
    exm.order         <- sort(train.label$label, index.return=TRUE)$ix
    train.label$label <- train.label$label[exm.order]
    train.feat        <- train.feat[exm.order,]

    train.label.exp   <- sample(train.label$label)
    foldid            <- assign.fold(label = train.label.exp, num.folds, stratified = stratify, inseparable=inseparable, meta=meta)

    ### internal cross-validation for model selection
    model             <- train.plm(data=data, method = method, subset=fold.exm.idx[[r]])
    if(!all(model$feat.weights == 0)){
       models.list[[r]]  <- model
    }

    ### TODO Important: the 'mh' variable gets written into the coefficient matrix.
    ### This needs to be changed ASAP, as the check in plm_predictor.r is obsolete with a hard-coded string like this.
    mh = paste('#LASSO (L1-regularized logistic regression (L1R_LR)',': [BINARY:',
               label$negative.lab, '=negative',';',
               label$positive.lab, '=positive', ']', sep='')
    ### collect model parameters (feature weights)
    if (r==1) {
      model.header = mh
    } else {
      stopifnot(model.header == mh)
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
      vec <- rep(NA, nrow(beta) + 2)
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
    } else if(method == "libLineaR"){
            # Liblinear needs C, W (intercept term is included in W).
      # Furthermore, it needs an element called "ClassNames" which is important in determining which class label is positive or negative.
      vec <- rep(NA, length(models.list[[i]]$original.model$W) + 3)
      vec[1] <- 0
      vec[2:3] <- as.numeric(models.list[[i]]$original.model$ClassNames)
      vec[4:length(vec)] <- as.numeric(models.list[[i]]$original.model$W)
      if (i == 1) {
        out.matrix <- matrix(vec)
      } else {
        out.matrix <- cbind(out.matrix, vec)
      }
      # This overwrites rownames everytime, but doesnt need an additional conditional statement.
      # paste0 pastes two equal-length string vectors element-wise
      # Note that the weight-vector is a row-vector here.
      rownames(out.matrix) <- c("C", "negative label", "positive label", colnames(models.list[[i]]$original.model$W))

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
  invisible(list(out.matrix=out.matrix, model.header=model.header, W.mat=W.mat, models.list=models.list))
}
