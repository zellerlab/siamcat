#!/usr/bin/Rscript
###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# R flavor
# EMBL Heidelberg 2012-2018
# GNU GPL 3.0
###

#' @title Prediction on the test set
#' @description This function takes the test set instances and the model trained by \link{plm.trainer} in order to predict the classes.
#' @param siamcat object of class \link{siamcat-class}
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'        for only information about progress and success, \code{2} for normal
#'        level of information and \code{3} for full debug information, defaults to \code{1}
#' @export
#' @keywords SIAMCAT plm.predictor
#' @return object of class \link{siamcat-class}
#'
make.predictions <- function(siamcat, siamcat.holdout=NULL, normalize.holdout=TRUE, verbose=1){

  s.time <- proc.time()[3]

  # if holdout is NULL, make predictions on data in siamcat
  if (is.null(siamcat.holdout)){

    if(verbose>1) cat("+ starting make.predictions on siamcat object\n")

    feat <- t(matrix(siamcat@phyloseq@otu_table)
    label.fac <- factor(siamcat@label@label, levels=c(siamcat@label@negative.lab, siamcat@label@positive.lab))

    # assert that there is a split
    stopifnot(!is.null(siamcat@data_split))

    num.folds <- siamcat@data_split@num.folds
    num.resample <- siamcat@data_split@num.resample

    pred <- matrix(NA, ncol=num.resample, nrow=length(label.fac),
                   dimnames=list(names(label.fac), paste0('CV_rep', 1:num.resample)))
    i = 1
    if(verbose==1 || verbose==2) pb <- txtProgressBar(max=num.folds*num.resample, style=3)
    for (f in 1:num.folds){
      for (r in 1:num.resample){

        test.label <- label.fac[siamcat@data_split@test.folds[[r]][[f]]]
        data       <- as.data.frame(feat[siamcat@data_split@test.folds[[r]][[f]],])

        # assert stuff
        stopifnot(nrow(data) == length(test.label))
        stopifnot(all(rownames(data) == names(test.label)))

        data$label <- test.label
        model <- siamcat@model_list@models[[i]]

        stopifnot(!any(rownames(model$task$env$data) %in% rownames(data)))
        if(verbose>2) cat('Applying ', siamcat@model_list@model.type, ' on cv_fold', f, '_rep', r,
            ' (', i, ' of ', num.resample*num.folds, ')...\n', sep='')

        task <- makeClassifTask(data = data, target = "label")
        pdata <- predict(model,  task = task)

        # rescale posterior probabilities between -1 and 1 (this works only for binary data!!!!)
        # TODO: Will need adjustment and generalization in the future)
        p <- siamcat@label@negative.lab+abs(siamcat@label@positive.lab-siamcat@label@negative.lab)*pdata$data[,4]
        names(p) <- rownames(pdata$data)
        pred[names(p), r] <- p
        i <- i+1
        if(verbose==1 || verbose==2) setTxtProgressBar(pb, i)
      }
    }
    stopifnot(!any(is.na(pred)))
    siamcat@pred_matrix <- pred
    return.object <- siamcat
  } else {

    if(verbose>1) cat("+ starting make.predictions on external dataset\n")

    if (normalize.holdout){
      if(verbose>1) cat("+ Performing frozen normalization on holdout set\n")
      siamcat.holdout <- normalize.features(siamcat.holdout, norm_param=siamcat@norm_param, verbose=verbose)
    } else {
      cat("WARNING: holdout set is not being normalized!\n")
    }
    feat.test <- t(siamcat.holdout@phyloseq@otu_table)
    feat.ref <- t(siamcat@phyloseq@otu_table)

    # data sanity checks
    stopifnot(all(colnames(feat.ref) %in% colnames(feat.test)))

    # prediction
    num.models <- siamcat@data_split@num.folds * siamcat@data_split@num.resample

    pred <- matrix(NA, ncol=num.models, nrow=nrow(feat.test),
                   dimnames=list(rownames(feat.test), paste0('Model_', 1:num.models)))
    if(verbose==1 || verbose==2) pb <- txtProgressBar(max=num.folds*num.resample, style=3)
    for (i in 1:num.models){

     data <- as.data.frame(feat.test)
     model <- siamcat@model_list@models[[i]]

     data <- data[,model$features]
     data$label <- as.factor(siamcat.holdout@label@label)

     if(verbose>2) cat('Applying ', siamcat@model_list@model.type, ' on complete external dataset',
         ' (', i, ' of ', num.models, ')...\n', sep='')

     task <- makeClassifTask(data = data, target = "label")
     pdata <- predict(model,  task = task)

     p <- siamcat@label@negative.lab+abs(siamcat@label@positive.lab-siamcat@label@negative.lab)*pdata$data[,4]
     names(p) <- rownames(pdata$data)
     pred[names(p), i] <- p
     if(verbose==1 || verbose==2) setTxtProgressBar(pb, i)
    }
    return.object <- pred
  }

  # print correlation matrix
  if(verbose>1) cat('Total number of predictions made:', length(pred), '\n')
  correlation <- cor(pred, method='spearman')
  if(verbose>1) cat('Correlation between predictions from repeated CV:\n')
  if(verbose>1) cat('Min: ', min(correlation), ', Median: ', median(correlation), ', Mean: ', mean(correlation), '\n', sep='')

  # print out time
  e.time <- proc.time()[3]
  if(verbose>1) cat("+ finished make.predictions in",e.time-s.time,"s\n")
  if(verbose==1) cat("\nMade predictions successfully.\n")

  return(return.object)
}
