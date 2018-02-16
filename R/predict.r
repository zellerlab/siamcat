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
#' @param siamcat object of class \link{siamcat-class}
#' @export
#' @keywords SIAMCAT plm.predictor
#' @return object of class \link{siamcat-class}
#'
make.predictions <- function(siamcat, verbose=1){


  feat         <- t(siamcat@phyloseq@otu_table)

  ### subselect training examples as specified in fn.train.sample (if given)
  foldList     <- get.foldList(siamcat@dataSplit, siamcat@label, mode="test", model=siamcat@modelList)
  fold.name    <- foldList$fold.name
  fold.exm.idx <- foldList$fold.exm.idx
  num.runs     <- foldList$num.runs
  num.folds    <- foldList$num.folds

  ### apply one LASSO model per test sample (i.e. CV fold)
  # predictions are made on a concatenation of test examples from all test samples
  pred = NULL
  predList = list()
  fold.pred.idx = list()

  # Init hyperpar list
  opt.hp <- list(lambda = NULL, C = NULL, alpha = NULL, ntree = NULL)
  if(verbose==1 || verbose==2) pb <- txtProgressBar(max=num.runs, style=3)
  for (r in 1:num.runs) {
    label.fac         <- factor(siamcat@label@label, levels=c(siamcat@label@negative.lab, siamcat@label@positive.lab))
    test.label        <- label.fac
    test.label        <- label.fac[fold.exm.idx[[r]]]
    data              <- as.data.frame(feat[fold.exm.idx[[r]],])
    stopifnot(nrow(data)         == length(test.label))
    stopifnot(all(rownames(data) == names(test.label)))
    data$label                     <- test.label
    model <- siamcat@modelList@models[[r]]
    if(verbose>2) cat('Applying ', siamcat@modelList@model.type, ' on ', fold.name[r], ' (', r, ' of ', num.runs, ')...\n', sep='')
    # subselect appropriate model

    # subselect test examples
    test.feat = feat[fold.exm.idx[[r]],,drop=FALSE]
    task      <- makeClassifTask(data = data, target = "label")
    pdata    <- predict(model,  task = task)
    # save(pdata,file="pdata.RData") # not very good, saves the data somewhere, depending on working directory
    p        <- siamcat@label@negative.lab+abs(siamcat@label@positive.lab-siamcat@label@negative.lab)*pdata$data[,4]
    names(p) <- rownames(pdata$data)

    pred     <- c(pred, p)
    fold.pred.idx[[r]] = (length(pred)-length(p)+1):length(pred)
    if(verbose==1 || verbose==2) setTxtProgressBar(pb, r)
  }

  if(verbose) cat('\nTotal number of predictions made:', length(pred), '\n')

  if (!is.null(siamcat@dataSplit)) {
    ### if test labels are given do some evaluation as well
    # get the appropriate labels for all test sets
    test.label = NULL
    aucs = vector('numeric', num.runs)
    for (r in 1:num.runs) {
      lab        <- siamcat@label@label[fold.exm.idx[[r]]]
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
    if (length(siamcat@label@label) == length(pred) && all(names(siamcat@label@label) %in% names(pred)) && all(names(pred) %in% names(siamcat@label@label))) {
      m = match(names(siamcat@label@label), names(pred))
      pred = pred[m]
      test.label = test.label[m]
      stopifnot(all(names(siamcat@label@label) == names(pred)))
    }

    # test accuracy of combined test set
    c.auc = NA
    if (length(unique(test.label)) == 2) {
      ev = eval.classifier(pred, as.vector(test.label), label)
      c.auc = calc.auroc(ev)
    }
    if(verbose)cat('Combined test AUC = ', format(c.auc, digits=3),
        ' (m=', format(mean(aucs, na.rm=TRUE), digits=3),
        ', s.d.=', format(sd(aucs, na.rm=TRUE), digits=3), ')\n', sep='')
  }


  ### reformat predictions in case models were trained in repeated cross-validation
  if (length(unique(names(pred))) < length(pred)) {
    ref.names = NULL
    if (any(substr(fold.name,1,14) == 'whole data set')) {
      r.idx = c(1:num.runs) #as.numeric(sapply(strsplit(fold.name, 'predicted by model '), '[[', 2))
      runs = sort(unique(r.idx))
      stopifnot(all(runs == 1:length(runs)))
      if (!is.null(siamcat@dataSplit)) {
        ref.names = names(siamcat@label@label)
      } else {
        ref.names = unique(names(pred))
      }
    } else {
      r.idx = as.numeric(sapply(strsplit(fold.name, 'rep'), '[[', 2))
      runs = sort(unique(r.idx))
      stopifnot(all(runs == 1:length(runs)))
      if (!is.null(siamcat@dataSplit)) {
        ref.names = names(siamcat@label@label)
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
      #if (!is.null(siamcat@dataSplit)) {
      #   stopifnot(all(sort(m) == 1:length(label$label)))
      #}
      pred.mat[m,r] = pred[p]
      stopifnot(all(names(pred)[p] == rownames(pred.mat)[m]))
    }
    correlation <- cor(pred.mat, method='spearman')
    if(verbose)cat('\nCorrelation between predictions from repeated CV:\n')
    if(verbose)cat('Min: ', min(correlation), ', Median: ', median(correlation), ', Mean: ', mean(correlation), '\n', sep='')
  }else{
    pred.mat = as.matrix(pred,byrow=TRUE)
  }
  #print(pred.mat[1:3,1:3])
  siamcat@predMatrix <- pred.mat
  return(siamcat)
}
