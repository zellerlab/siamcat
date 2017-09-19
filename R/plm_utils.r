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


##### function to draw a stratified sample from the label vector
#' @export
sample.strat = function(label, frac) {
  classes = unique(label)
  s = NULL
  for (c in classes) {
    idx = which(label==c)
    s = c(s, sample(which(label==c), round(frac * sum(label==c))))
  }
  return(s)
}

##### function to train a LASSO model for a single given C
#' @export
train.plm <- function(feat, label, method, hyper.par, data, subset) {
  method <- tolower(method)
  model <- list(original.model=NULL, feat.weights=NULL)
  # note that some of the logit models are set up inversely to each other,
  # requiring to swap coefficient signs
  if (method == 'lasso') {
    # here we will ignore any hyper.par$alpha to ensure a LASSO model (and not an Elastic Net) is trained
    lambda              <- hyper.par$lambda
    # Note that ordering of label vector is important!
    #print(label.fac)
    #saveList <- list(feat=feat, label.fac=label.fac)
    #save(saveList, file="featwhat.RData")
    #model$original.model <- glmnet(x=feat, y=label.fac, family='binomial', standardize=FALSE, alpha=1, lambda=lambda)

    ## 1) Define the task
    ## Specify the type of analysis (e.g. classification) and provide data and response variable
    task                <- makeClassifTask(data = data, target = "cancer")

    ## 2) Define the learner
    ## Choose a specific algorithm (e.g. linear discriminant analysis)
    lrn                 <- makeLearner("classif.cvglmnet", predict.type="prob")

    ## 3) Fit the model
    ## Train the learner on the task using a random subset of the data as training set
    model               <- train(lrn, task)

    coef                <- coefficients(model$learner.model)
    bias.idx            <- which(rownames(coef) == '(Intercept)')
    coef                <- coef[-bias.idx,]
    
    # Remove bias/intercept term when returning model feature weights
    model$feat.weights <- (-1) *  as.numeric(coef) ### check!!!
    model$lrn          <- lrn
    model$task         <- task
    
  } else if (method == 'lasso_ll') {
    liblin.type <- 6
    model$original.model <- LiblineaR(feat, label$label, 
                                      type=liblin.type, bias=TRUE, epsilon=1e-6,
                                      cost=hyper.par$C)
    stopifnot((is.numeric(label$positive.lab) &&  is.numeric(label$negative.lab)) || label$positive.lab < label$negative.lab)
    # Apparently, LiblineaR orders class labels differently, depending on "type". I need to check this more thoroughly. Works for now.
    if(model$original.model$ClassNames[2] != label$positive.lab){
      # Since prediction function is mirrored on the y-axis, we also need to generate "swapped" models.
      # This means model$ClassNames has to have structure c(NL, PL) and the model features have to be mirrored as well (neg. numbers become pos. and vice verca)
      # If the upper is true, ensure that swapping takes place.
      temp = model$original.model$ClassNames[1]
      model$original.model$ClassNames[1] = model$original.model$ClassNames[2]
      model$original.model$ClassNames[2] = temp
      model$original.model$W = model$original.model$W * -1
    }
    bias.idx <- which(colnames(model$original.model$W) == 'Bias')
    model$feat.weights <- as.numeric(model$original.model$W[,-bias.idx])
    
  } else if (method == 'ridge_ll') {
    liblin.type <- 0
    model$original.model <- LiblineaR(feat, label$label, 
                                      type=liblin.type, bias=TRUE, epsilon=1e-6,
                                      cost=hyper.par$C)
    stopifnot((is.numeric(label$postive.lab) &&  is.numeric(label$negative.lab)) || label$postive.lab < label$negative.lab)
    # Apparently, LiblineaR orders class labels differently, depending on "type". I need to check this more thoroughly. Works for now.
    if(model$original.model$ClassNames[2] != label$postive.lab){
      # Since prediction function is mirrored on the y-axis, we also need to generate "swapped" models.
      # This means model$ClassNames has to have structure c(NL, PL) and the model features have to be mirrored as well (neg. numbers become pos. and vice verca)
      # If the upper is true, ensure that swapping takes place.
      temp = model$original.model$ClassNames[1]
      model$original.model$ClassNames[1] = model$original.model$ClassNames[2]
      model$original.model$ClassNames[2] = temp
      model$original.model$W = model$original.model$W * -1
    }
    bias.idx <- which(colnames(model$original.model$W) == 'Bias')
    model$feat.weights <- as.numeric(model$original.model$W[,-bias.idx])
    
  } else if (method == 'enet') {
    stopifnot(hyper.par$alpha < 1) # otherwise we're going to train a LASSO model
    model$original.model <- glmnet(feat, factor(label$label, levels=c(label$negative.lab, label$postive.lab)),
                                   family='binomial', standardize=FALSE,
                                   alpha=hyper.par$alpha, lambda=hyper.par$lambda)
    c <- coefficients(model$original.model)
    c <- c[,ncol(c)]
    bias.idx <- which(names(c) == '(Intercept)')
    model$feat.weights <- -1 * as.numeric(c[-bias.idx])    
    
  } else if (method == 'gelnet') {
    stopifnot(hyper.par$alpha < 1) # otherwise we're going to train a LASSO model
    model$original.model <- gelnet(feat, factor(label$label, levels=c(label$negative.lab, label$postive.lab)),
                                   l2=hyper.par$alpha, l1=hyper.par$lambda,
                                   silent=TRUE)
    # for consistency with the other models the coefficient sign needs to be flipped
    model$feat.weights <- model$original.model$w
    
  } else {
    stop('unknown method')
  }
  
  return(model)
}


#' @export
predict.plm <- function(feat, model, method, opt.hyper.par, data, subset) {
  method <- tolower(method)
  
  # note that some of the logit models are set up inversely to each other,
  # requiring to select coefficient/prediction columns accordingly using col.idx 
  if (method == 'lasso') {
    col.idx <- 1
    # glmnet's predict function needs to be given a lambda value
    pred <- predict(model, task = model$task, subset = subset)
    
  } else if (method == 'lasso_ll' || method == 'ridge_ll') {
    col.idx <- 2 # this is a bit counter-intuitive given the column names of the predict data frame
    pred <- predict(model$original.model, feat, 
                    proba=TRUE)$probabilities[,col.idx]
    # Adding rownames here is important.
    names(pred) <- rownames(feat)
    
  } else if (method == 'enet') {
    col.idx <- 1
    # glmnet's predict function needs to be given a lambda value 
    pred <- predict(model$original.model, feat,
                    alpha=opt.hyper.par$alpha, s=opt.hyper.par$lambda,
                    type="response")[,col.idx]
    
  } else if (method == 'gelnet') {
    # gelnet's predictions need to be calculated "by hand"
    m = model$original.model
    pred <- 1.0 / (1.0 + exp(feat %*% m$w + m$b))
    names(pred) <- rownames(feat)
    
  } else {
    stop('unknown method')
  }
  return(pred)
}

#' @export
select.model <- function(feat, label, method, hyper.par, min.nonzero=1,
                         num.folds=5, stratified=FALSE, foldid=foldid, data) {
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
      cat('    ', method, ' model selection: (lambda=', hyper.par$lambda[i],
          ') AU-ROC=', format(aucs[i], digits=3), '\n', sep='')
    }
    suff.nonzero         <- apply(nonzero.coeff, 2, min) > min.nonzero
    nonzero.coeff        <- nonzero.coeff[,suff.nonzero]
    aucs                 <- aucs[suff.nonzero]
    opt.idx              <- which.max(aucs)
    opt.hyper.par$lambda <- hyper.par$lambda[opt.idx]
    cat('    optimal lambda =', hyper.par$lambda[opt.idx], '\n')
    
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
assign.fold <- function(label, num.folds, stratified, inseparable = NULL, foldid) {
  classes <- sort(unique(label))
  # Transform number of classes into vector of 1 to x for looping over.
  # stratify positive examples
  if (stratified) {
    # If stratify is TRUE, make sure that num.folds does not exceed the maximum number of examples for the class with the fewest training examples.
    if (any(as.data.frame(table(label))[,2] < num.folds)) {
      print("+++ Number of CV folds is too large for this data set to maintain stratification. Reduce num.folds or turn stratification off. Exiting.\n")
      q(status = 1)
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
      strata      <- unique(meta.data[,inseparable])
      sid         <- sample(rep(1:num.folds, length.out=length(strata)))
      foldid      <- rep(NA, length(label))
      for (s in 1:length(strata)) {
        idx          <- which(meta.data[,inseparable] == strata[s])
        foldid[idx]  <-sid[s]
      }
      stopifnot(all(!is.na(foldid)))
    } else {
      foldid      <- sample(rep(1:num.folds, length.out=length(label)))
    }
    
  }
  # make sure each fold contains examples from all classes
  for (f in 1:num.folds) {
    stopifnot(all(sort(unique(label[foldid==f])) == classes))
  }
  stopifnot(length(label) == length(foldid))
  return(foldid)
}

##### end function
