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
#' @param test.samples filename containing the test samples or list of test instances produced by \link{data.splitter()}, defaults to \code{NULL} leading to testing on the complete dataset
#' @param model model object trained by \link{plm.trainer}
#' @param model.mat model matrix needed to rebuild the model
#' @param hyperpars (not used)
#' @param model.type string, type of the model that was trained
#' @export
#' @keywords SIAMCAT plm.predictor
#' @return list containing the precitions \itemize{
#'  \item \code{$pred};
#'  \item \code{$mat}
#'}
plm.predictor <- function(feat, label, test.samples=NULL, models.list, model.mat, hyperpars, model.type){
  # TODO hyperpars is not used at the moment, as far as i see
  # TODO 2: instead of feat and label containing the test sample indices, provide all features and the list from data.splitter
  feat         <- t(feat)
  label.fac                  <- factor(label$label, levels=c(label$negative.lab, label$positive.lab))

  data                       <- cbind(feat,label.fac[rownames(feat)])
  data                       <- as.data.frame(data)
  data[,ncol(data)]          <- as.factor(data[,ncol(data)])
  colnames(data)             <- paste0("Sample_",1:ncol(data))
  colnames(data)[ncol(data)] <- "cancer"

  num.runs = length(models.list)
  ### subselect test examples as specified in test.samples (if given)
  fold.name = list()
  fold.exm.idx = list()
  if (is.null(test.samples)){
    # apply each LASSO model on whole data set when only single test set is given
    for (r in 1:num.runs) {
      fold.name[[r]] = paste('whole data set predicted by model', r)
      fold.exm.idx[[r]] = rownames(feat)
    }
  } else {
    if (class(test.samples) == 'character'){
      con = file(test.samples, 'r')
      input = readLines(con)
      m.idx = 0
      for (i in 1:length(input)) {
        l = input[[i]]
        if (substr(l, 1, 1) != '#') {
          m.idx = m.idx + 1
          s = unlist(strsplit(l, '\t'))
          fold.name[[m.idx]] = substr(s[1], 2, nchar(s[1]))
          fold.exm.idx[[m.idx]] = which(rownames(feat) %in% as.vector(s[2:length(s)]))
          #      cat(fold.name[[m.idx]], 'contains', length(fold.exm.idx[[m.idx]]), 'test examples\n')
          #      cat(fold.exm.idx[[m.idx]], '\n\n')
        }
      }
      close(con)
      stopifnot(length(fold.name) == num.runs)
      stopifnot(length(fold.exm.idx) == num.runs)
      stopifnot(all(paste('M', unlist(fold.name), sep='_') == colnames(models.list[[1]]$W)))
  } else if (class(test.samples) == 'list') {
    i = 1
    for (cv in 1:test.samples$num.folds){
      for (res in 1:test.samples$num.resample){
        fold.name[[i]] = paste0('cv_fold', as.character(cv), '_rep', as.character(res))
        fold.exm.idx[[i]] <- match(test.samples$test.folds[[res]][[cv]], names(label$label))
        i = i + 1
        # cat(fold.name[[num.runs]], 'contains', length(fold.exm.idx[[num.runs]]), 'training examples\n')
      }
    }
    }
  }
  fold.name = unlist(fold.name)
  cat('\nPreparing to make predictions with', num.runs,   model.type, ' model(s)...\n')

  ### apply one LASSO model per test sample (i.e. CV fold)
  # predictions are made on a concatenation of test examples from all test samples
  pred = NULL
  predList = list()
  fold.pred.idx = list()

  # Init hyperpar list
  opt.hp <- list(lambda = NULL, C = NULL, alpha = NULL, ntree = NULL)

  for (r in 1:num.runs) {
    model <- models.list[[r]]
    cat('Applying ', colnames(model$W)[r], ' on ', fold.name[r], ' (', r, ' of ', num.runs, ')...\n', sep='')
    curr.model.new <- list()
    curr.model.new$original.model <- list()

    # Rebuild model
    if (model.type == 'lasso'){
      curr.model.new$original.model$lambda <- model.mat[1, r]
      curr.model.new$original.model$a0 <- model.mat[2, r]
      curr.model.new$original.model$beta <- model.mat[3:dim(model.mat)[1], r]
      # glmnet needs "offset" value
      curr.model.new$original.model$offset <- FALSE
      # glmnet needs sparse matrix
      curr.model.new$original.model$beta <- Matrix(curr.model.new$original.model$beta, sparse = TRUE)
      class(curr.model.new$original.model) <- c("lognet", "glmnet")

      # Create hp list
      opt.hp$lambda <- curr.model.new$original.model$lambda

    } else if (model.type == 'enet') {
      curr.model.new$original.model$lambda <- model.mat[1, r]
      curr.model.new$original.model$lambda <- model.mat[2, r]
      curr.model.new$original.model$a0 <- model.mat[3, r]
      curr.model.new$original.model$beta <- model.mat[4:dim(model.mat)[1], r]
      # glmnet needs "offset" value
      curr.model.new$original.model$offset <- FALSE
      # glmnet needs sparse matrix
      curr.model.new$original.model$beta <- Matrix(curr.model.new$original.model$beta, sparse = TRUE)
      class(curr.model.new$original.model) <- c("lognet", "glmnet")

      # Create hp list
      opt.hp$lambda <- curr.model.new$original.model$lambda
      opt.hp$alpha <- curr.model.new$original.model$alpha

    } else if (model.type == 'lasso_ll'){
      curr.model.new$original.model$W <- t(as.matrix(model.mat[4:dim(model.mat)[1], r]))
      curr.model.new$original.model$Bias <- "TRUE"
      curr.model.new$original.model$Type <- 6
      curr.model.new$original.model$TypeDetail <- "L1-regularized logistic regression (L1R_LR)"
      curr.model.new$original.model$ClassNames <- model.mat[2:3, r]
      curr.model.new$original.model$NbClass <- 2
      class(curr.model.new$original.model) <- c("LiblineaR")

      # Create hp list
      opt.hp$C <- curr.model.new$original.model$C

    } else if (model.type == 'ridge_ll'){
      curr.model.new$original.model$W <- t(as.matrix(model.mat[4:dim(model.mat)[1], r]))
      curr.model.new$original.model$Bias <- "TRUE"
      curr.model.new$original.model$Type <- 0
      curr.model.new$original.model$TypeDetail <- "L2-regularized logistic regression primal (L2R_LR)"
      curr.model.new$original.model$ClassNames <- model.mat[2:3, r]
      curr.model.new$original.model$NbClass <- 2
      class(curr.model.new$original.model) <- c("LiblineaR")

      # Create hp list
      opt.hp$C <- curr.model.new$original.model$C
    } else if (model.type == 'gelnet') {
      curr.model.new$original.model$alpha <- as.vector(na.omit(model.mat[,1]))
      curr.model.new$original.model$lambda <- as.vector(na.omit(model.mat[,2]))
      curr.model.new$original.model$b <- as.vector(na.omit(model.mat[3,r]))
      curr.model.new$original.model$w <- as.vector(na.omit(model.mat[4:dim(model.mat)[1], r]))

      # Create hp list
      opt.hp$alpha <- curr.model.new$original.model$alpha
      opt.hp$lambda <- curr.model.new$original.model$lambda
    }
    # subselect appropriate model
    m = models.list[[r]]
    m$W = m$W[,r]

    # subselect test examples
    test.feat = feat[fold.exm.idx[[r]],,drop=FALSE]

    pdata    <- predict.plm(test.feat, m, model.type, opt.hp, data = data, subset=fold.exm.idx[[r]])
    p        <- label$negative.lab+abs(label$positive.lab-label$negative.lab)*pdata$data[,4]
    names(p) <- rownames(pdata$data)

    pred     <- c(pred, p)
    fold.pred.idx[[r]] = (length(pred)-length(p)+1):length(pred)
  }

  cat('\nTotal number of predictions made:', length(pred), '\n')

  if (!is.null(test.samples)) {
    ### if test labels are given do some evaluation as well
    # get the appropriate labels for all test sets
    test.label = NULL
    aucs = vector('numeric', num.runs)
    for (r in 1:num.runs) {
      lab        <- label$label[fold.exm.idx[[r]]]
      test.label <- c(test.label, lab)
      lab.p.idx  <- (length(test.label)-length(lab)+1):length(test.label)
      # accuracy of individual test sets
      if (length(unique(test.label[lab.p.idx])) == 2) {
        ev = eval.classifier(pred[lab.p.idx], as.vector(test.label[lab.p.idx]), label)
        aucs[r] = calc.auroc(ev)
      }
    }
    stopifnot(length(test.label) == length(pred))
    stopifnot(names(test.label) == names(pred))

    # in case of cross-validation there should be exactly one prediction per labeled example,
    # so we reorder them according to the order of label
    if (length(label$label) == length(pred) && all(names(label$label) %in% names(pred)) && all(names(pred) %in% names(label$label))) {
      m = match(names(label$label), names(pred))
      pred = pred[m]
      test.label = test.label[m]
      stopifnot(all(names(label$label) == names(pred)))
    }

    # test accuracy of combined test set
    c.auc = NA
    if (length(unique(test.label)) == 2) {
      ev = eval.classifier(pred, as.vector(test.label), label)
      c.auc = calc.auroc(ev)
    }
    cat('Combined test AUC = ', format(c.auc, digits=3),
        ' (m=', format(mean(aucs, na.rm=TRUE), digits=3),
        ', s.d.=', format(sd(aucs, na.rm=TRUE), digits=3), ')\n', sep='')
  }


  ### reformat predictions in case models were trained in repeated cross-validation
  if (length(unique(names(pred))) < length(pred)) {
    ref.names = NULL
    if (any(substr(fold.name,1,14) == 'whole data set')) {
      r.idx = as.numeric(sapply(strsplit(fold.name, 'predicted by model '), '[[', 2))
      runs = sort(unique(r.idx))
      stopifnot(all(runs == 1:length(runs)))
      if (!is.null(test.samples)) {
        ref.names = names(label$label)
      } else {
        ref.names = unique(names(pred))
      }
    } else {
      r.idx = as.numeric(sapply(strsplit(fold.name, 'rep'), '[[', 2))
      runs = sort(unique(r.idx))
      stopifnot(all(runs == 1:length(runs)))
      if (!is.null(test.samples)) {
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
      if (!is.null(test.samples)) {
        stopifnot(all(sort(m) == 1:length(label$label)))
      }
      pred.mat[m,r] = pred[p]
      stopifnot(all(names(pred)[p] == rownames(pred.mat)[m]))
    }
    correlation <- cor(pred.mat, method='spearman')
    cat('\nCorrelation between predictions from repeated CV:\n')
    cat('Min: ', min(correlation), ', Median: ', median(correlation), ', Mean: ', mean(correlation), '\n', sep='')
  }else{
    pred.mat = as.matrix(pred,byrow=TRUE)
  }
  #print(pred.mat[1:3,1:3])
  invisible(list(pred = pred, mat = pred.mat))
}
