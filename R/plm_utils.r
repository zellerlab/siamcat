###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# RScript flavor
#
# written by Georg Zeller
# with additions by Nicolai Karcher and Konrad Zych
# EMBL Heidelberg 2012-2017
#
# version 0.2.0
# file last updated: 25.06.2017
# GNU GPL 3.0
###

#liblinear.type = 6
#class.weights  = c(1, 1) # TODO reset!!!
#model.tag      = 'LASSO'
#eps            = 1e-8

# TODO: to be removed, since it is not used by any other function
##### function to draw a stratified sample from the label vector
# #' @export
# sample.strat = function(label, frac) {
#   classes = unique(label)
#   s = NULL
#   for (c in classes) {
#     idx = which(label==c)
#     s = c(s, sample(which(label==c), round(frac * sum(label==c))))
#   }
#   return(s)
# }

##### function to train a LASSO model for a single given C
#' @export
train.plm <- function(data, method = c("lasso", "enet", "ridge", "libLineaR", "randomForest"), subset) {
  #model <- list(original.model=NULL, feat.weights=NULL)

  ## 1) Define the task
  ## Specify the type of analysis (e.g. classification) and provide data and response variable
  task      <- makeClassifTask(data = data, target = "label")

  ## 2) Define the learner
  ## Choose a specific algorithm (e.g. linear discriminant analysis)
  cl        <- "classif.cvglmnet" ### the most ocommon learner defined her so taht this does not have done multiple times
  paramSet  <- NULL

  if(method == "lasso"){
    lrn       <- makeLearner(cl, predict.type="prob", 'nlambda'=10, 'alpha'=1)
  }else if(method == "ridge"){
    lrn       <- makeLearner(cl, predict.type="prob", 'nlambda'=10, 'alpha'=0)
  }else if(method == "enet"){
    lrn       <- makeLearner(cl, predict.type="prob", 'nlambda'=10)
    paramSet  <- makeParamSet(makeNumericParam('alpha', lower=0, upper=1))
  }else if(method == "libLineaR"){
    cl        <- "classif.LiblineaRL1LogReg"
    lrn       <- makeLearner(cl, predict.type="prob")
  } else if(method == "randomForest"){
    sqrt.mdim <- sqrt(nrow(data))
    cl        <- "classif.randomForest"
    lrn       <- makeLearner(cl, predict.type = "prob", fix.factors.prediction = TRUE)
    paramSet  <- makeParamSet(makeNumericParam('ntree', lower=100, upper=1000),
                        makeDiscreteParam('mtry', values=c(round(sqrt.mdim/2), round(sqrt.mdim), round(sqrt.mdim*2))))
  } else {
    stop(method, " is not a valid method, currently supported: lasso, enet, ridge, libLineaR, randomForest.\n")
  }
  

  ## 3) Fit the model
  ## Train the learner on the task using a random subset of the data as training set
  if(!all(is.null(paramSet))){
    hyperPars <- tuneParams(learner = lrn, 
                         task = task,
                         resampling =  makeResampleDesc('CV', iters=5L, stratify=TRUE),
                         par.set = paramSet, 
                         control = makeTuneControlGrid(resolution = 10L), 
                         measures=list(acc))
    lrn       <- setHyperPars(lrn, par.vals=hyperPars$x)
  }

  model     <- train(lrn, task, subset=subset)

  if(cl == "classif.cvglmnet"){
    coef                <- coefficients(model$learner.model)
    bias.idx            <- which(rownames(coef) == '(Intercept)')
    coef                <- coef[-bias.idx,]
    model$feat.weights  <- (-1) *  as.numeric(coef) ### check!!!
  } else if(cl == "classif.LiblineaRL1LogReg"){
    model$feat.weights  <-model$learner.model$W[-which(colnames(model$learner.model$W)=="Bias")]
  } else if(cl == "classif.randomForest"){
    model$feat.weights  <-model$learner.model$importance
  }

  model$lrn          <- lrn ### ???
  model$task         <- task

return(model)
}


#' @export
select.model <- function(feat, label, method, hyper.par, min.nonzero=1,
                         num.folds=5, stratified=FALSE, foldid=foldid, data) {
  print(method)
  method        <- tolower(method)
  opt.hyper.par <- NULL
  nonzero.coeff <- matrix(0, num.folds, length(hyper.par$lambda))
  # here we use the area under the ROC curve as a model selection criterion,
  # but there are alternatives (e.g. accuracy, partial AU-ROC or the area under the precision-recall curve)
  fold.id       <- foldid
  p             <- rep(NA, length(label$label))
  if (method == 'lasso') {
    aucs <- rep(0, length(hyper.par$lambda))
    for (i in 1){#:length(hyper.par$lambda)) {
      hp <- NULL
      hp$lambda <- hyper.par$lambda[i]
      for (f in 1:num.folds) {
        test.idx <- which(fold.id == f)
        # Replace this with train.idx <- which(fold.id != f)? For better readability.
        train.idx         <- setdiff(1:length(label$label), test.idx)
        label.train       <- label
        label.train$label <- label.train$label[train.idx]
        #print(feat[train.idx,])
        model             <- train.plm(feat[train.idx,], label.train, method, hp, data=data, subset=train.idx)
        pred              <- predict.plm(feat[test.idx,], model, method, hp, data=data, subset=test.idx)
        #print(model)
        nonzero.coeff[f,i] = sum(model$feat.weights[1:(ncol(feat)-1)] != 0)
      }
      aucs[i]   <- performance(pred, measures = auc)#roc(response=label$label, predictor=p)$auc
      #cat('    ', method, ' model selection: (lambda=', hyper.par$lambda[i],
      #    ') AU-ROC=', format(aucs[i], digits=3), '\n', sep='')
      print(performance(pred, measures = auc))
      cat("    AU-ROC = ", format(aucs[i], digits=3), '\n', sep='')
    }
    suff.nonzero         <- apply(nonzero.coeff, 2, min) > min.nonzero
    nonzero.coeff        <- nonzero.coeff[,suff.nonzero]
    aucs                 <- aucs[suff.nonzero]
    opt.idx              <- which.max(aucs)
    opt.hyper.par$lambda <- hyper.par$lambda[opt.idx]
    # cat('    optimal lambda =', hyper.par$lambda[opt.idx], '\n')

  } else if (method == 'lasso_ll' || method == 'ridge_ll' || method == 'l1_svm' || method == 'l2_svm') {
    aucs <- rep(0, length(hyper.par$C))
    for (i in 1:length(hyper.par$C)) {
      hp <- NULL
      hp$C <- hyper.par$C[i]
      for (f in 1:num.folds) {
        test.idx    <- which(fold.id == f)
        train.idx   <- setdiff(1:length(label$label), test.idx)
        label.train <- label
        label.train$label <- label.train$label[train.idx]
        model       <- train.plm(feat[train.idx,], label.train, method, hp)
        p[test.idx] <- predict.plm(feat[test.idx,], model, method, hp)
      }
      aucs[i] <- roc(response=label, predictor=p)$auc
      cat('    ', method, ' model selection: (C=', hyper.par$C[i],
          ') AU-ROC=', format(aucs[i], digits=3), '\n', sep='')
    }
    opt.idx <- which.max(aucs)
    opt.hyper.par$C <- hyper.par$C[opt.idx]
    cat('    optimal C =', hyper.par$C[opt.idx], '\n')

  } else if (method == 'enet' || method == 'gelnet') {
    aucs <- matrix(0, nrow=length(hyper.par$alpha),
                   ncol=length(hyper.par$lambda))
    for (i in 1:length(hyper.par$alpha)) {
      for (j in 1:length(hyper.par$lambda)) {
        hp <- NULL
        hp$alpha <- hyper.par$alpha[i]
        hp$lambda <- hyper.par$lambda[j]
        for (f in 1:num.folds) {
          test.idx <- which(fold.id == f)
          train.idx <- setdiff(1:length(label$label), test.idx)
          label.train <- label
          label.train$label <- label.train$label[train.idx]
          model       <- train.plm(feat[train.idx,], label.train, method, hp, data, train.idx)
          p[test.idx] <- predict.plm(feat[test.idx,], model, method, hp, data, test.idx)
        }
        aucs[i,j] <- roc(response=label, predictor=p)$auc
        cat('    ', method, ' model selection: (alpha=', hyper.par$alpha[i],
            ', lambda=', hyper.par$lambda[j], ') AU-ROC=',
            format(aucs[i,j], digits=3), '\n', sep='')
      }
    }
    opt.idx              <- arrayInd(which.max(aucs), dim(aucs))
    opt.hyper.par$alpha  <- hyper.par$lambda[opt.idx[1]]
    opt.hyper.par$lambda <- hyper.par$lambda[opt.idx[2]]
    cat('    optimal alpha =' , hyper.par$alpha[opt.idx[1]],
        ', lambda =', hyper.par$lambda[opt.idx[2]], '\n')

  } else {
    stop('unknown method')
  }
  return(opt.hyper.par)
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
    print(num.folds)
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
      print("+++ num.samples is exceeding number of folds, setting CV to (k-1) unstratified CV\n")
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
