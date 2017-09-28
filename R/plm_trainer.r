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
#' @param training.samples filename containing the training samples or list of training instances produced by \link{data.splitter()}, defaults to \code{NULL} leading to training on the complete dataset
#' @param num.folds integer, number of folds for the internal model cross-validation
#' @param stratify boolean, should the folds in the internal cross-validation be stratified?
#' @param modsel.crit model selection criterion (not used at the moment)
#' @param min.nonzero.coeff integer number of minimum nonzero coefficients that should be present in the model
#' @export
#' @keywords SIAMCAT plm.trainer
#' @return list containing \itemize{
#' \item \code{$out.matrix};
#' \item \code{$model.header};
#' \item \code{$W.mat};
#' \item \code{$hyperpar.mat};
#' \item \code{$model}
#'}
# TODO add details section for this function
plm.trainer <- function(feat, label, training.samples=NULL, num.folds=5, stratify, modsel.crit, min.nonzero.coeff, model.type='lasso', inseparable=NULL, meta=NULL){
  # TODO 1: modsel.criterion should be implemented
  # TODO 2: instead of filename containing the traning sample indices, provide the list from data.splitter
  # TODO 3: add model.type as parameter
  # transpose feature matrix as a convenience preprocessing
  feat         <- t(feat)
  label.fac                  <- factor(label$label, levels=c(label$negative.lab, label$positive.lab))

  data                       <- cbind(feat,label.fac[rownames(feat)])
  data                       <- as.data.frame(data)
  data[,ncol(data)]          <- as.factor(data[,ncol(data)])
  colnames(data)             <- paste0("Sample_",1:ncol(data))
  colnames(data)[ncol(data)] <- "cancer"
  ### subselect training examples as specified in fn.train.sample (if given)
  num.runs     <- 1
  fold.name    <- list()
  fold.exm.idx <- list()

  # print(fn.train.sample)
  if (is.null(training.samples)){
    # train on whole data set
    fold.name[[1]]    <- 'whole data set'
    fold.exm.idx[[1]] <- names(label$label)
  } else {
    if (class(training.samples) == 'character') {
      # read in file containing the training instances
      num.runs      <- 0
      con           <- file(training.samples, 'r')
      input         <- readLines(con)
      close(con)
      print(length(input))
      for (i in 1:length(input)) {
        l               <- input[[i]]
        if (substr(l, 1, 1) != '#') {
          num.runs                 <- num.runs + 1
          print(num.runs)
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
    } else if (class(training.samples) == 'list') {
      # use training samples as specified in training.folds in the list
      num.runs <- 0
      for (cv in 1:training.samples$num.folds){
        for (res in 1:training.samples$num.resample){
          num.runs <- num.runs + 1
          print(num.runs)

          fold.name[[num.runs]] = paste0('cv_fold', as.character(cv), '_rep', as.character(res))
          fold.exm.idx[[num.runs]] <- match(training.samples$training.folds[[res]][[cv]], names(label$label))
          cat(fold.name[[num.runs]], 'contains', length(fold.exm.idx[[num.runs]]), 'training examples\n')
        }
      }
    } else {
      stop('Wrong input for training samples!...')
    }
  }
  print(num.runs)
  #stop()
  fold.name     <- unlist(fold.name)
  stopifnot(length(fold.name) == num.runs)
  stopifnot(length(fold.exm.idx) == num.runs)
  cat('\nPreparing to train', model.type,  'models on', num.runs, 'training set samples...\n\n')

  ### train one model per training sample (i.e. CV fold)
  # feat has structure: examples in rows; features in columns!
  W.mat           <- matrix(data=NA, nrow=ncol(feat), ncol=num.runs)
  rownames(W.mat) <- c(colnames(feat))
  colnames(W.mat) <- paste('M', fold.name, sep='_')

  # Create matrix with hyper parameters.
  hyperpar.list   <- list()

  # Create List to save models.
  models.list     <- list()

  for (r in 1:num.runs) {
    cat('Training on ', fold.name[r], ' (', r, ' of ', num.runs, ')...\n', sep='')
    ### subselect examples for training
    train.label       <- label
    train.label$label <- train.label$label[fold.exm.idx[[r]]]
    train.feat        <- feat[fold.exm.idx[[r]],]
    stopifnot(nrow(train.feat)         == length(train.label$label))
    stopifnot(all(rownames(train.feat) == names(train.label$label)))

    # reorder training examples so that class order is the same for all models
    #print(train.label$label)
    exm.order         <- sort(train.label$label, index.return=TRUE)$ix
    train.label$label <- train.label$label[exm.order]
    train.feat        <- train.feat[exm.order,]


    ### assign training data to internal folds for model selection
    ### For structure of foldid, see data_splitter.r
    # foldid            <- rep(0, length(train.label$label))
    # perm              <- sample(1:length(train.label$label), length(train.label$label)) / length(train.label$label)
    # for (f in num.folds:1) {
    #   foldid[perm <= f/num.folds] = f
    # }
    # ####
    # Uncommented the lines above, because i wasn't sure what exactly they do. Foldid will be overwritten by assing.fold anyway, won't it?
    # ####
    train.label.exp   <- sample(train.label$label)
    foldid            <- assign.fold(label = train.label.exp, num.folds, stratified = stratify, inseparable=inseparable, meta=meta)

    ### internal cross-validation for model selection
    sqrt.mdim         <- sqrt(nrow(feat))
    hyper.par <- list(C      = c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10),  # all LiblineaR methods
                      lambda = c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3),            # LASSO & ENet
                      alpha  = c(0.7,0.8,0.9)                                    # ENet
                      #                 ntree  = c(250, 500, 1000),                                # RF (not functional for now)
                      #                 mtry   = c(sqrt.mdim/2, sqrt.mdim, sqrt.mdim*2)            # RF (not functional for now)
    )
    opt.hyper.par    <- select.model(train.feat, train.label, model.type, hyper.par, min.nonzero = min.nonzero.coeff,
                                     num.folds=num.folds, stratified=FALSE, foldid=foldid, data=data)
    # cat('  optimal C=', opt.hyper.par$lambda, ' (', which(opt.C$lambda==C.vec), ' of ', length(hyper.par$lambda), ')\n', sep='')
    ### retrain whole training set with best parameter setting (after every internal CV run!)
    model            <- train.plm(train.feat, train.label, model.type, opt.hyper.par, data=data, subset=fold.exm.idx[[r]])
    models.list[[r]] <- model


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
    save(model,file="model.RData")
    stopifnot(all(names(model$W) == rownames(W.mat)))
    W.mat[,r]          <- as.numeric(c(model$feat.weights))
    hyperpar.list[[r]] <- unlist(opt.hyper.par)
    stopifnot(!all(model$feat.weights == 0))
    cat('\n')
  }
  # Preprocess hyper parameters
  hyperpar.mat           <- matrix(unlist(hyperpar.list), ncol = num.runs, nrow = length(hyperpar.list[[1]]), byrow = FALSE)
  print(dim(hyperpar.mat))
  rownames(hyperpar.mat) <- names(hyperpar.list[[1]])
  ### Write models into matrix to reload in plm_predictor.r
  for (i in 1:length(models.list)){
    if (model.type == 'lasso') {
      # glmnet needs lambda, a0 and beta.
      vec <- rep(NA, nrow(models.list[[i]]$learner.model$glmnet.fit$beta) + 2)
      print(i)
      vec[1] <- hyperpar.mat[1, i]
      vec[2] <- models.list[[i]]$learner.model$glmnet.fit$a0
      vec[3:length(vec)] <- as.numeric(models.list[[i]]$learner.model$glmnet.fit$beta)
      if (i == 1) {
        out.matrix <- matrix(vec)
      } else {
        out.matrix <- cbind(out.matrix, vec)
      }
      # This overwrites rownames everytime, but doesnt need an additional conditional statement.
      # paste0 pastes two equal-length string vectors element-wise
      rownames(out.matrix) <- c("lambda", "a0", rownames(models.list[[i]]$learner.model$glmnet.fit$beta))
      # In the case of glmnet, make the coefficient matrix from a sparse matrix into a regular one.
      #models.list[[i]]$original.model$beta <- as.matrix(models.list[[i]]$learner.model$glmnet.fit$beta)
    } else if (model.type == 'enet'){
      # glmnet needs lambda, a0 and beta.
      vec <- rep(NA, length(models.list[[i]]$original.model$beta) + 3)
      vec[1:2] <- hyperpar.mat[, i]
      vec[3] <- models.list[[i]]$original.model$a0
      vec[4:length(vec)] <- as.numeric(models.list[[i]]$original.model$beta)
      if (i == 1) {
        out.matrix <- matrix(vec)
      } else {
        out.matrix <- cbind(out.matrix, vec)
      }
      # This overwrites rownames everytime, but doesnt need an additional conditional statement.
      # paste0 pastes two equal-length string vectors element-wise
      rownames(out.matrix) <- c("lambda", "alpha", "a0", rownames(models.list[[i]]$original.model$beta))
      # In the case of glmnet, make the coefficient matrix from a sparse matrix into a regular one.
      models.list[[i]]$original.model$beta <- as.matrix(models.list[[i]]$original.model$beta)
    } else if (model.type == 'lasso_ll' || model.type == 'ridge_ll') {
      # Liblinear needs C, W (intercept term is included in W).
      # Furthermore, it needs an element called "ClassNames" which is important in determining which class label is positive or negative.
      vec <- rep(NA, length(models.list[[i]]$original.model$W) + 3)
      vec[1] <- hyperpar.mat[1, i]
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

    } else if (model.type == 'gelnet') {
      # gelnet needs b (intercept), w, as well as alpha and lambda.
      vec <- rep(NA, length(models.list[[i]]$original.model$w) + 3)
      vec[1:2] <- hyperpar.mat[1, i]
      vec[3] <- models.list[[i]]$original.model$b
      vec[4:length(vec)] <- as.numeric(models.list[[i]]$original.model$w)
      if (i == 1) {
        out.matrix <- matrix(vec)
      } else {
        out.matrix <- cbind(out.matrix, vec)
      }
      # This overwrites rownames everytime, but doesnt need an additional conditional statement.
      # paste0 pastes two equal-length string vectors element-wise
      # Note that the feature weight containing object is a vector
      rownames(out.matrix) <- c("alpha", "lambda", "b", names(models.list[[i]]$original.model$w))
    }
  }
  colnames(out.matrix) = paste('M', fold.name, sep='_')
  invisible(list(out.matrix=out.matrix, model.header=model.header, W.mat=W.mat, hyperpar.mat=hyperpar.mat, model=model))
}
