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

#' @title Validate samples in metadata, labels, and features
#' @description This function checks if metadata is available for all samples in labels and vice versa. If yes, order metadata according to the order found in labels.
#' @param feat features object
#' @param label labels object
#' @param meta metadata object
#' @keywords SIAMCAT validate.data
#' @export
#' @return a list containing values after filtering: \itemize{
#'  \item \code{$feat} = features;
#'  \item \code{$label} = labels;
#'  \item \code{$meta} = metadata
#' }
validate.data <- function(feat, label, meta = NULL){
  # TODO attempt multiple reading attempts!
  # TODO 2: should this function also validate the features? (at the moment, nothing happens to the features...)
  # TODO 3: this functions takes label$label at the moment, not an label object. also does not return an label obejct. Should be changed, i guess?

  if (!is.null(meta)) {
  stopifnot(all(names(label) %in% rownames(meta)) && all(rownames(meta) %in% names(label)))
  m    <- match(names(label), rownames(meta))
  meta <- meta[m,]
  stopifnot(all(names(label) == rownames(meta)))
  return(list("meta"=meta,
              "label"=label,
              "feat"=feat))
  }
  return(list("label"=label,
              "feat"=feat))
}
  ### run data checks / validation / format conversion
  # TODO !!!

#' @title Check for potential confounders in the metadata
#' @description This function checks for associations between class labels and potential confounders (e.g. age, sex, or BMI) that are present in the metadata.
#' Statistical testing is performed with Fisher's exact test or Wilcoxon test, while associations are visualized either as barplot or Q-Q plot, depending on the type of metadata.
#' @param meta metadata object
#' @param label labels object
#' @keywords SIAMCAT confounder.check
#' @export
#' @examples
#' pdf(filename.pdf)
#' confounder.check(meta, label)
#' dev.off()
#' @return Does not return anything, but produces a single plot for each metadata category
confounder.check <- function(meta, label){

  for (m in 1:ncol(meta)) {
    mname <- gsub('[_.-]', ' ', colnames(meta)[m])
    mname <- paste(toupper(substring(mname, 1, 1)), substring(mname, 2), sep="")
    cat('checking', mname, 'as a potential confounder...\n')

    mvar  <- as.numeric(meta[,m])
    u.val <- unique(mvar)
    u.val <- u.val[!is.na(u.val)]
    if (length(u.val) == 1) {
      cat('  skipped because all subjects have the same value\n')
    } else if (length(u.val) <= 5) {
      cat('  using a bar plot\n')
      par(mar=c(6.1,4.1,4.1,4.1))
      ct     <- matrix(NA, nrow=2, ncol=length(u.val))
      colors <- rainbow(length(u.val), s=0.4, v=0.6, alpha=1, start=0.0)
      for (i in 1:length(u.val)) {
        ct[1,i] = sum(mvar[label$n.idx] == u.val[i], na.rm=TRUE)
        ct[2,i] = sum(mvar[label$p.idx] == u.val[i], na.rm=TRUE)
      }
      freq   <- t(ct)
      for (i in 1:dim(freq)[2]) {
        freq[,i] <- freq[,i] / sum(freq[,i])
      }
      barplot(freq, ylim=c(0,1), main=mname, names.arg=c(label$n.lab, label$p.lab), col=colors)
      p.val  <- fisher.test(ct)$p.value
      mtext(paste('Fisher test p-value:', format(p.val, digits=4)), side=1, line=3, at=1, adj=0)
    } else {
      cat('  using a Q-Q plot\n')
      par(mar=c(5.1,4.1,4.1,4.1))
      ax.int <- c(min(mvar, na.rm=TRUE), max(mvar, na.rm=TRUE))
      qqplot(mvar[label$n.idx], mvar[label$p.idx], xlim=ax.int, ylim=ax.int, pch=16, cex=0.6,
             xlab=label$n.lab, ylab=label$p.lab, main=paste('Q-Q plot for', mname))
      abline(0, 1, lty=3)
      p.val  <- wilcox.test(mvar[label$n.idx], mvar[label$p.idx], exact=FALSE)$p.value
      text(ax.int[1]+0.9*(ax.int[2]-ax.int[1]), ax.int[1]+0.1*(ax.int[2]-ax.int[1]),
           paste('MWW test p-value:', format(p.val, digits=4)), pos=2)
    }
  }
}

#' @title Perform unsupervised feature filtering.
#' @description This function may convert absolute abundances into relative abundances and then performs unsupervised feature filtering. Features can be filtered based on abundance or prevalence. Additionally, unmapped reads may be removed.
#' @param feat feature object
#' @param filter.method method used for filtering the features, can be one of these: \code{c("abundance", "cum.abundance", "prevalence")}
#' @param cutoff float, abundace or prevalence cutoff
#' @param recomp.prop boolean, should absolute abundances be converted into relative abundances?
#' @param rm.unmapped boolean, should unmapped reads be discarded?
#' @details Currently, there are three filtering methods implemented:
#' \itemize{
#'  \item \code{"abundance"} remove features for which the abundance is never above the threshold value (e.g. 0.5\%) in any of the samples
#'  \item \code{"cum.abundance"} remove features with very low abundance in all samples, i.e. features that are never among the most abundant entities that collectively make up \code{(1-cutoff)} of the reads in any samples
#'  \item \code{"prevalence"} remove features with low prevalence across samples, i.e. features that are undetected in more than \code{(1-cutoff)} proportion of samples
#' }
#' @keywords SIAMCAT filter.feat
#' @export
#' @return Returns the filtered feature matrix
filter.feat <- function(feat, filter.method, cutoff, recomp.prop, rm.unmapped){
  ### this statement does not have the purpose to calculate relative abundances on the fly and return them.
  ### Instead, it's purpose is to be able to calculate f.idx (specifying the indices of features which are to be kept)
  ### when feature list has already been transformed to relative abundances, but e.g. certain features have been removed manually.
  ## TODO check filter.method, add default value for cutoff, recomp.prop, and rm.unmapped?
  if (recomp.prop) {
    # recompute relative abundance values (proportions)
    ra.feat = prop.table(feat, 2)
  } else {
    ra.feat = feat
  }

  ### apply filters
  if (filter.method == 'abundance') {
    # remove features whose abundance is never above the threshold value (e.g. 0.5%) in any of the samples
    f.max = apply(ra.feat, 1, max)
    f.idx = which(f.max >= cutoff)
  } else if (filter.method == 'cum.abundance') {
    # remove features with very low abundance in all samples i.e. ones that are never among the most abundant
    # entities that collectively make up (1-cutoff) of the reads in any sample
    f.idx = vector('numeric', 0)
    # sort features per sample and apply cumsum to identify how many collectively have weight K
    for (s in 1:ncol(ra.feat)) {
      srt = sort(ra.feat[,s], index.return=TRUE)
      cs = cumsum(srt$x)
      m = max(which(cs < cutoff))
      f.idx = union(f.idx, srt$ix[-(1:m)])
    }
    # an index of those features that collectively make up more than 1-K of the read mass in any sample
    f.idx = sort(f.idx)
  } else if (filter.method == 'prevalence') {
    # remove features with low prevalence across samples
    # i.e. ones that are 0 (undetected) in more than (1-cutoff) proportion of samples
    f.idx = which(rowSums(ra.feat > 0) / ncol(ra.feat) > cutoff)
  } else {
    cat('\nunrecognized filter.method, exiting!\n')
    quit(save='no', status = 1)
  }

  cat('Removed ', nrow(feat)-length(f.idx), ' features whose values did not exceed ', cutoff,
      ' in any sample (retaining ', length(f.idx), ')\n', sep='')
  feat = feat[f.idx,]

  ### postprocessing and output generation
  if (rm.unmapped) {
    # remove 'unmapped' feature
    unm.idx = rownames(feat) == 'UNMAPPED' | rownames(feat) == 'unmapped' | rownames(feat) == '-1' | rownames(feat) == 'UNCLASSIFIED' | rownames(feat) == 'unclassified' | rownames(feat) == 'UNASSIGNED' | rownames(feat) == 'unassigned'
    if (any(unm.idx)) {
      feat = feat[!unm.idx,]
      cat('Removed ', sum(unm.idx), ' features corresponding to UNMAPPED reads',
          ' (retaining ', nrow(feat), ')\n', sep='')
    } else {
      cat('tried to remove unmapped reads, but could not find them. Continue anyway.')
    }
  }
  return(feat)
}

#' @title Perform feature normalization according to specified parameters.
#' @description This function performs feature normalization according to user-specified parameters.
#' @param feat feature object
#' @param norm.method normalization method, can be one of these: \code{c("rank.unit", "rank.std", "log.std", "log.unit")}
#' @param log.n0 pseudocount to be added before log-transformation
#' @param sd.min.q quantile of the distribution of standard deviation of all feature that will be added to the denominator during standardization of each feature in order to avoid underestimation (only for methods=\code{c("rank.std", "log.std")})
#' @param n.p vector norm to use, can be either \code{1} or \code{2}, (only for method=\code{"log.unit"})
#' @param n.sample boolean, normalize by feature? (only for method=\code{"log.unit"})
#' @param n.feature boolean, normalize by sample? (only for method=\code{"log.unit"})
#' @param n.global boolean, Normalize by global rescaling? (only for method=\code{"log.unit"})
#' @details There are four different normalization methods implemented:
#' \itemize{
#'  \item \code{"rank.unit"} converts features to ranks and normalizes each column by the square root of the sum of ranks
#'  \item \code{"rank.std"} converts features to ranks and applies z-score standardization
#'  \item \code{"log.unit"} log-transforms features (after addition of pseudocounts), converts features to ranks and normalizes by features or samples
#'  \item \code{"log.std"} log-transforms features (after addition of pseudocounts) and applies z-score standardization
#' }
#' The other parameters are depending on the normalization method of choice, e.g. \code{log.n0} as pseudocount is only needed for the methods making use of the log-transformation.
#' @keywords SIAMCAT normalize.feat
#' @export
#' @return list containing the matrix of normalized features and a list of normalization parameters: \itemize{
#'  \item \code{$par} = parameters utilized in the normalization;
#'  \item \code{$feat} = normalized features
#'}
normalize.feat <- function(feat, norm.method, log.n0, sd.min.q, n.p, n.sample, n.feature, n.global) {
  ### remove features with missing values
  # TODO there may be better ways of dealing with NA features
  # TODO 2 add defaults for the parameters!!! Not all parameters are needed for all normalization methods
  num.orig.feat = nrow(feat)
  keep.idx = rowSums(is.na(feat) == 0)
  if (any(!keep.idx)) {
    feat = feat[keep.idx,]
    cat('Removed ', nrow(feat)-num.orig.feat, ' features with missing values (retaining ', nrow(feat),' )\n', sep='')
  }

  ### keep track of normalization parameters
  par = list()
  par$norm.method = norm.method
  par$log.n0 = log.n0
  par$n.p = n.p
  par$n.sample = n.sample
  par$n.feature = n.feature
  par$n.global = n.global
  par$retained.feat = rownames(feat)

  ### apply normalization
  if (norm.method == 'rank.unit') {
    for (c in 1:ncol(feat)) {
      feat[,c] = rank(feat[,c], ties.method='average')
    }
    stopifnot(!any(is.na(feat)))
    for (c in 1:ncol(feat)) {
      feat[,c] = feat[,c] / sqrt(sum(feat[,c]^2))
    }
  } else if (norm.method == 'rank.std') {
    for (c in 1:ncol(feat)) {
      feat[,c] = rank(feat[,c], ties.method='average')
    }
    m = apply(feat, 1, mean)
    s = apply(feat, 1, sd)
    q = quantile(s, sd.min.q, names=FALSE)
    stopifnot(q > 0)
    # TODO needs an apply-style rewrite!
    for (r in 1:nrow(feat)) {
      feat[r,] = (feat[r,] - m[r]) / (s[r] + q)
    }
    par$feat.mean = m
    par$feat.adj.sd = s + q
    stopifnot(!any(is.na(feat)))
  } else if (norm.method == 'log.std') {
    feat = log10(feat + log.n0)
    m = apply(feat, 1, mean)
    s = apply(feat, 1, sd)
    q = quantile(s, sd.min.q, names=FALSE)
    #cat(sort(s, decreasing=TRUE), '\n')
    stopifnot(q > 0)
    # TODO needs an apply-style rewrite!
    for (r in 1:nrow(feat)) {
      feat[r,] = (feat[r,] - m[r]) / (s[r] + q)
    }
    par$feat.mean = m
    par$feat.adj.sd = s + q
    stopifnot(!any(is.na(feat)))
  } else if (norm.method == 'log.unit') {
    cat('Feature sparsity before normalization: ', 100*mean(feat==0), '%\n', sep='')
    feat = log10(feat + log.n0)
    if (n.p == 1) {
      if (n.feature) {
        feat.norm.denom = vector('numeric', nrow(feat))
        # TODO needs an apply-style rewrite!
        for (r in 1:nrow(feat)) {
          feat.norm.denom[r] = sum(feat[r,])
          feat[r,] = feat[r,] / feat.norm.denom[r]
        }
        par$feat.norm.denom = feat.norm.denom
      }
      if (n.sample) {
        for (c in 1:ncol(feat)) {
          feat[,c] = feat[,c] / sum(feat[,c])
        }
      }
    } else if (n.p == 2) {
      if (n.feature) {
        feat.norm.denom = vector('numeric', nrow(feat))
        for (r in 1:nrow(feat)) {
          feat.norm.denom[r] = sqrt(sum(feat[r,]^2))
          feat[r,] = feat[r,] / feat.norm.denom[r]
        }
        par$feat.norm.denom = feat.norm.denom
      }
      if (n.sample) {
        for (c in 1:ncol(feat)) {
          feat[,c] = feat[,c] / sqrt(sum(feat[,c]^2))
        }
      }
    } else {
      stop('unknown norm!')
    }
    if (!n.feature && !n.sample && n.global) {
      global.norm.denom = max(feat)
      feat = feat / global.norm.denom
      par$global.norm.denom = global.norm.denom
    }
    cat('Feature sparsity after normalization: ', 100*mean(feat==0), '%\n', sep='')
    stopifnot(!any(is.na(feat)))
  } else {
    cat('\nunrecognized norm.method, exiting!\n')
    quit(save='no', status = 1)
  }
  return(list("par" = par, "feat" = feat))
}


#' @title Split a dataset into training and a test sets.
#' @description This function prepares the cross-validation by splitting the data into \code{num.folds} training and test folds for \code{num.resample} times.
#' @param label label object
#' @param num.folds number of cross-validation folds (needs to be \code{>=2})
#' @param num.resample resampling rounds (values \code{<= 1} deactivate resampling)
#' @param stratify boolean, should the splits be stratified s. t. an equal proportion of classes are present in each fold?
#' @param inseparable defaults to NULL (obsolete)
#' @param meta metadata object (obsolete)
#' @keywords SIAMCAT data.splitter
#' @export
#' @return list containing the indices of the training and test folds and the parameters of the splits: \itemize{
#'  \item \code{$training.folds} nested list, containing for \code{length(num.folds)} the sample names of the \code{length(num.resample)} training folds;
#'  \item \code{$test.folds} nested list, containing for \code{length(num.folds)} the sample names of the \code{length(num.resample)} test folds;
#'  \item \code{$num.resample} = number of repeated samplings;
#'  \item \code{$num.folds} = number of folds
#'}
data.splitter <- function(label, num.folds, num.resample, stratify, inseparable, meta){
  ### read label and meta-data
  # (assuming the label file has 1 column)
  # TODO delete inseparable parameter (or at least add an default)
  # TODO delete meta parameter

  if (is.null(inseparable) || inseparable=='' || toupper(inseparable)=='NULL' || toupper(inseparable)=='NONE' || toupper(inseparable)=='UNKNOWN') {
    inseparable <- NULL
    cat('+++ Inseparable parameter not specified\n')
  }
  labelNum        <- as.numeric(label$label)
  names(labelNum) <- names(label$label)
  exm.ids         <- names(labelNum)

  # parse label description
  classes      <- sort(label$info$class.descr)


  ### check arguments
  if (num.resample < 1) {
    cat('+++ Resetting num.resample = 1 (', num.resample, ' is an invalid number of resampling rounds)\n', sep='')
    num.resample   <- 1
  }
  if (num.folds < 2) {
    cat('+++ Resetting num.folds = 2 (', num.folds, ' is an invalid number of folds)\n', sep='')
    num.folds      <- 2
  }
  if (!is.null(inseparable) && stratify) {
    cat('+++ Stratification is not supported when inseparable is given\n')
    stratify       <- FALSE
  }
  if (num.folds >= length(labelNum)) {
    cat('+++ Performing leave-one-out (LOO) cross-validation\n')
    if (stratify) {
      cat('   WARNING: Stratification is not possible with LOO cross-validation\n')
      stratify     <- FALSE
    }
    num.folds   <- length(labelNum)
  }


  train.list <- list(NULL)
  test.list  <- list(NULL)


  for (r in 1:num.resample) {
    labelNum      <- sample(labelNum)
    foldid        <- rep(0, length(labelNum))
    foldid        <- assign.fold(label = labelNum, num.folds, stratified = stratify, inseparable = NULL, foldid = foldid)
    names(foldid) <- names(labelNum)
    stopifnot(length(labelNum) == length(foldid))
    stopifnot(length(unique(foldid)) == num.folds)

    train.temp    <- list(NULL)
    test.temp     <- list(NULL)

    cat('\n+++ Splitting the dataset:\n')
    for (f in 1:num.folds) {
      # make sure each fold contains examples from all classes
      stopifnot(all(sort(unique(labelNum[foldid==f])) == classes))

      # select test examples
      test.idx        <- which(foldid == f)
      train.idx       <- which(foldid != f)
      train.temp[f] <- list(names(foldid)[train.idx])
      test.temp[f]  <- list(names(foldid)[test.idx])
      stopifnot(length(intersect(train.idx, test.idx)) == 0)
      cat('   + Fold ', f, ' contains ', sum(foldid==f), ' examples\n', sep='')
    }
    train.list[[r]] <- train.temp
    test.list[[r]]  <- test.temp
  }


  return(list("training.folds" = train.list, "test.folds" = test.list, "num.resample" = num.resample, "num.folds" = num.folds))
}

#' @title Evaluate prediction results
#' @description This function takes the correct labels and predictions for all samples and evaluates the results using the \itemize{
#'  \item Area Under the Receiver Operating Characteristic (ROC) Curve (AU-ROC)
#'  \item and the Precision-recall Curve (PR)
#' }
#' as metric. Predictions can be supplied either for a single case or as matrix after resampling of the dataset.
#' @param label label object
#' @param pred prediction for each sample by the model, should be a matrix with dimensions \code{length(label) x 1} or \code{length(label) x num.resample}
#' @keywords SIAMCAT eval.result
#' @export
#' @return list containing \itemize{
#'  \item \code{$rov.average};
#'  \item \code{$auc.average};
#'  \item \code{$ev.list};
#'  \item \code{$pr.list};
#' }. If \code{prediction} had more than one column, the function will additonally return \itemize{
#'  \item \code{$auc.all};
#'  \item \code{$aucspr}
#'}
eval.result <- function(label, pred){

  # TODO compare header to label
  ### make sure that label and prediction are in the same order
  #stopifnot(all(names(label) %in% rownames(pred)) && all(rownames(pred) %in% names(label)))
  m    <- match(names(label$label), rownames(pred))
  #cat(m, '\n')
  pred <- pred[m,,drop=FALSE]
  stopifnot(all(names(label$label) == rownames(pred)))

  # ROC curve
  auroc = 0
  if (ncol(pred) > 1) {
    rocc = list(NULL)
    aucs = vector('numeric', ncol(pred))
    for (c in 1:ncol(pred)) {
      rocc[c] = list(roc(response=label$label, predictor=pred[,c], ci=FALSE))
      aucs[c] = rocc[[c]]$auc
    }
    l.vec = rep(label$label, ncol(pred))
  } else {
    l.vec = label$label
  }
  # average data for plotting one mean prediction curve
  summ.stat = 'mean'
  rocsumm = list(roc(response=label$label, predictor=apply(pred, 1, summ.stat),
                 ci=TRUE, of="se", sp=seq(0, 1, 0.05)))
  auroc = list(rocsumm[[1]]$auc)
  # precision recall curve
  pr = list(NULL)
  ev = list(NULL)
  if (ncol(pred) > 1) {
    aucspr = vector('numeric', dim(pred)[2])
    for (c in 1:ncol(pred)) {
      print(c)
      ev[c] = list(eval.classifier(pred[,c], label$label, label))
      pr[c] = list(get.pr(ev[[c]]))
      aucspr[c] = calc.aupr(ev[[c]])
    }
    ev = append(ev,list(eval.classifier(apply(pred, 1, summ.stat), label$label, label)))
  } else {
    ev[1] = list(eval.classifier(as.vector(pred), label$label, label))
    pr[1] = list(get.pr(ev[[1]]))
  }
  if (ncol(pred) > 1) {
    return(list("roc.all" = rocc,
                "auc.all"=aucs,
                "roc.average"=rocsumm,
                "auc.average"=auroc,
                "ev.list"=ev,
                "pr.list"=pr,
                "aucspr"=aucspr))
  } else {
    return(list("roc.average"=rocsumm,
                "auc.average"=auroc,
                "ev.list"=ev,
                "pr.list"=pr))
  }
}

#' @title Select samples
#' @description This functions filters samples based on metadata.
#' Provided with a column in the metadata and a range or a set of allowed values, the function will return the labels, metadata, and features for the samples matching the allowed range or set.
#' @param meta metadata object
#' @param feat features object
#' @param label labels object
#' @param filter name of the column from metadata on which the selection should be done
#' @param allowed.set a vector of allowed values
#' @param allowed.range a range of allowed values
#' @keywords SIAMCAT select.samples
#' @export
#' @return list containing values after selection: \itemize{
#'  \item \code{$feat} = features;
#'  \item \code{$label} = labels;
#'  \item \code{$meta} = metadata
#' }
select.samples  <- function(meta, feat, label, filter, allowed.set = NULL, allowed.range = NULL){
  # TODO at the moment, this functions takes label$label, not an label-object. Should be corrected
  # parse interval
  if(!is.null(allowed.range )) {
    allowed.range  <- gsub('\\[|\\]','', allowed.range)
    allowed.range  <- as.numeric(unlist(strsplit(allowed.range,',')))
    stopifnot(length(allowed.range) == 2)
    stopifnot(allowed.range [1] <= allowed.range[2])
  }
  #allowed.set  <- opt$allowed.set # TODO should be removed, i guess?
  if (!is.null(allowed.set)) {
    # parse set
    allowed.set  <- gsub('\\{|\\}','', allowed.set)
    allowed.set  <- as.numeric(unlist(strsplit(allowed.set,',')))
    allowed.set  <- sort(unique(allowed.set))
  }
  if (!xor(is.null(allowed.range ), is.null(allowed.set))) {
    cat('Neither allowed.range  nor allowed.set (or both at the same time) have been provided, exiting!\n')
    quit(save='no', status = 1) # TODO this really sucks in Rstudio, maybe it can be less harsh?
  } else {
    if (!is.null(allowed.range )) {
      cat('allowed.range  = [', paste(allowed.range , collapse=','), ']\n', sep='')
    } else {
      cat('allowed.set = {', paste(allowed.set, collapse=','), '}\n', sep='')
    }
  }
  cat('\n')
  # TODO throws an error when called in Rstudio, should be removed i guess
  # if (substr(source.dir, nchar(source.dir), nchar(source.dir)) != '/') {
  #   source.dir   <- paste(source.dir, '/', sep='')
  # }
  if(!filter %in% colnames(meta)) stop("The filter name is not present in colnames of the metadata. Stopping.\n")
  filter.var    <- meta[,filter]

  if (!is.null(allowed.range )) {
    s.idx <- !is.na(filter.var) & filter.var >= allowed.range [1] & filter.var <= allowed.range [2]
    cat('Removed ', sum(!s.idx), ' samples with ', filter,
        ' not in [', paste(allowed.range , collapse=', '), '] (retaining ', sum(s.idx), ')\n', sep='')
  } else {
    s.idx <- !is.na(filter.var) & filter.var %in% allowed.set
    cat('Removed ', sum(!s.idx), ' samples with ', filter,
        ' not in {', paste(allowed.set, collapse=', '), '} (retaining ', sum(s.idx), ')\n', sep='')
  }
  results         <- list(NULL) # why?
  results$feat    <- feat[,s.idx]
  results$label   <- label[s.idx]
  results$meta    <- meta[s.idx,]
  # why invisible return?
  invisible(results)
}

#' @title Check and visualize associations between features and classes
#' @description This function calculates for each feature the median log10 fold change between the different classes found in labels.
#'
#' Significance of the differences is computed for each feature using a wilcoxon test followed by multiple hypothesis testing correction.
#'
#' Additionally, the Area under the Receiver Operating Characteristic Curve (AU-ROC) is computed for the features found to be associated with the two different classes at a user-specified significance level \code{alpha}.
#'
#' Finally, the function produces a plot of the top \code{max.show} associated features, showing the distribution of the log10-transformed abundances for both classes, the median fold change, and the AU-ROC.
#' @param feat feature object
#' @param label label object
#' @param fn.plot filename for the pdf-plot
#' @param color.scheme valid R color scheme, defaults to \code{'RdYlBu'}
#' @param alpha float, significance level, defaults to \code{0.05}
#' @param min.fc float, minimum log10 fold change for significantly associated features, defaults to \code{0}
#' @param mult.corr multiple hypothesis correction method, see \code{\link[stats]{p.adjust}}, defaults to \code{"fdr"}
#' @param detect.lim float, pseudocount to be added before log-transormation of the data, defaults to \code{1e-08}
#' @param max.show integer, how many associated features should be shown, defaults to \code{50}
#' @param plot.type string, specify how the abundance should be plotted, must be one of these: \code{c("bean", "box", "quantile.box", "quantile.rect")}, defaults to \code{"bean"}
#' @return Does not return anything, but produces an association plot
#' @keywords SIAMCAT check.associations
#' @export
check.associations <- function(feat, label, fn.plot, color.scheme="RdYlBu", alpha=0.05, min.fc=0, mult.corr="fdr",
                               detect.lim=10^-8, max.show=50, plot.type="bean"){

  sort.by <- 'pv'
  ### some color pre-processing
  if (color.scheme == 'matlab') {
    color.scheme <- matlab.like(100)
  } else {
    # TODO check for valid param!
    color.scheme <- rev(colorRampPalette(brewer.pal(11, color.scheme))(100))
  }
  col.p <- color.scheme[length(color.scheme)-4]
  col.n <- color.scheme[1+4]

  ### Define set of vectors that have the indeces and "description" of all positively and negatively labeled training examples.

  p.val <- vector('numeric', nrow(feat))
  fc    <- vector('numeric', nrow(feat))

  for (i in 1:nrow(feat)) {
    fc[i]    <- median(log10(feat[i,label$p.idx] + detect.lim)) - median(log10(feat[i,label$n.idx] + detect.lim))
    p.val[i] <- wilcox.test(feat[i,label$n.idx], feat[i,label$p.idx], exact = FALSE)$p.value
  }

  ### Apply multi-hypothesis testing correction
  if(!tolower(mult.corr) %in% c('none','bonferroni','holm','fdr','bhy')) stop('Unknown multiple testing correction method:', mult.corr,' Stopping!\n')
  if (mult.corr == 'none') {
    p.adj <- p.val
  } else {
    p.adj <- p.adjust(p.val, method=tolower(mult.corr))
  }

  cat('Found', sum(p.adj < alpha, na.rm=TRUE), 'significant associations at a significance level <', alpha, '\n')

  idx <- which(p.adj < alpha)
  if (min.fc > 0) {
    idx <- which(p.adj < alpha & abs(fc) > min.fc)
    cat('Found', length(idx), 'significant associations with absolute log10 fold change >', min.fc, '\n')
  }


  # TODO sort.by not an option of the function
  if (length(idx) > 0) {
    if (sort.by == 'fc') {
      idx <- idx[order(fc[idx], decreasing=FALSE)]
    } else if (sort.by == 'pv') {
      idx <- idx[order(p.adj[idx], decreasing=TRUE)]
    } else {
      cat('Unknown sorting option:', sort.by, 'order by p-value...\n')
      idx <- idx[order(p.adj[idx], decreasing=TRUE)]
    }
    for (i in idx) {
      cat(sprintf('%-40s', rownames(feat)[i]), 'p-value:', format(p.adj[i], digits=4), '\n')
    }
    # truncated the list for the following plots
    if (length(idx) > max.show) {
      idx <- idx[(length(idx)-max.show+1):length(idx)]
      cat('Truncating the list of significant associations to the top', max.show, '\n')
    }


    # compute single-feature AUCs
    cat('\nCalculating the area under the ROC curve for each significantly associated feature\n')
    aucs <- vector('numeric', nrow(feat))
    for (i in idx) {
      f       <- feat[i,]
      ev      <- eval.classifier(f, label$label, label)
      aucs[i] <- calc.auroc(ev)
      if (aucs[i] < 0.5) {
        aucs[i] <- 1-aucs[i]
      }
    }
    for (i in idx) {
      cat(sprintf('%-40s', rownames(feat)[i]), aucs[i], '\n')
    }

    ### generate plots with significant associations between features and labels
    pdf(fn.plot, paper='special', height=8.27, width=11.69) # format: A4 landscape

    lmat  <- cbind(1,2,3,4)
    layout(lmat, widths=c(0.6,0.075,0.2,0.2))

    x <- log10(as.matrix(feat[idx, label$p.idx, drop=FALSE]) + detect.lim)
    y <- log10(as.matrix(feat[idx, label$n.idx, drop=FALSE]) + detect.lim)

    col <- c(paste(col.n, '77', sep=''), paste(col.p, '77', sep=''), 'gray')
    if (plot.type == 'box') {
      par(mar=c(5.1, 25.1, 4.1, 0))
      box.colors <- rep(c(col[1],col[2]),nrow(x))
      plot.data <- data.frame()
      for (i in 1:nrow(x)){
        temp <- as.data.frame(rbind(cbind(x[i,],rep(paste(label$n.lab, rownames(x)[i]), length(x[i,]))), cbind(y[i,], rep(paste(label$p.lab, rownames(x)[i]), length((y[i,]))))))
        temp[,1] <- as.numeric(as.character(temp[,1]))
        plot.data <- rbind(plot.data, temp)
        if (i == nrow(x)) {
          plot(NULL, xlab='', ylab='',xaxs='i', yaxs='i', axes=FALSE,
               xlim=c(min(plot.data[,1]-0.2), max(plot.data[,1]) + 1), ylim=c(+0.5, length(idx)*2+0.5), type='n')
          boxplot(plot.data[,1] ~ plot.data[,ncol(plot.data)],horizontal=TRUE,
                  names = c(""), show.names = FALSE, col = box.colors, axes = FALSE, outcol = c(col[1], col[2]), add = TRUE)
          mn          <- as.integer(c(min(plot.data[,1])))
          mx          <- as.integer(c(max(plot.data[,1])))
          ticks       <- mn:mx
          for (v in ticks) {
            abline(v=v, lty=3, col='lightgrey')
          }
          tick.labels <- formatC(10^ticks, format='E', digits=0)
          axis(side=1, at=ticks, labels=tick.labels, cex.axis=0.7)
          ### function label.plot.horizontal has been written in utils.r.
          label.plot.horizontal(x, y, rownames(feat)[idx], x.suff=paste(' (', label$p.lab, ')', sep=''),
                                y.suff=paste(' (', label$n.lab, ')', sep=''), outer.diff = 2, inner.diff.x = 0, inner.diff.y = -1)
        }
      }
    }
    else if (plot.type == "quantile.box"){
      par(mar=c(5.1, 25.1, 4.1, 0))
      plot.data.range(x, y, rownames(feat)[idx], x.col=col[2], y.col=col[1])
      label.plot.horizontal(x, y, rownames(feat)[idx], x.suff=paste(' (', label$p.lab, ')', sep=''),
                            y.suff=paste(' (', label$n.lab, ')', sep=''), outer.diff = 1, inner.diff.x = 0.15, inner.diff.y = -0.15)
    }

    else if (plot.type == "quantile.rect"){
      par(mar=c(5.1, 25.1, 4.1, 0))
      quantiles.vector <- c(0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9)
      x.q = apply(x, 1, function (x) quantile(x, quantiles.vector, na.rm=TRUE, names=FALSE))
      x.medians = apply(x,1,function (x) median(x))
      y.q = apply(y, 1, function (y) quantile(y, quantiles.vector, na.rm=TRUE, names=FALSE))
      y.medians = apply(y,1,function (y) median(y))

      p.m = min(c(min(x, na.rm=TRUE), min(y, na.rm=TRUE)))
      plot(rep(p.m, dim(x)[1]), 1:dim(x)[1],
           xlab='', ylab='', yaxs='i', axes=FALSE,
           xlim=c(min(x,y), max(x,y+2)), ylim=c(0, dim(x)[1]), frame.plot=FALSE, type='n')
      for (v in seq(p.m,0,1)) {
        abline(v=v, lty=3, col='lightgrey')
      }

      tck = floor(p.m):0
      axis(1, tck, formatC(10^tck, format='E', digits=0), las=1, cex.axis=0.7)
      for (i in 1:(nrow(x.q)/2)){
        if (i == 1) {
          rect(x.q[i,], 0.5:dim(x)[1], x.q[nrow(x.q)+1-i,], (0.5:dim(x)[1])+0.3, col = c("white"), border = c("black"), lwd = 0.9)
          rect(y.q[i,], 0.5:dim(y)[1], y.q[nrow(y.q)+1-i,], (0.5:dim(y)[1])-0.3, col = c("white"), border = c("black"), lwd = 0.9)

        }
        else {
          rect(x.q[i,], 0.5:dim(x)[1], x.q[nrow(x.q)+1-i,], (0.5:dim(x)[1])+0.3, col = col[2], border = c("black"), lwd = 0.9)
          rect(y.q[i,], 0.5:dim(y)[1], y.q[nrow(y.q)+1-i,], (0.5:dim(y)[1])-0.3, col = col[1], border = c("black"), lwd = 0.9)
        }
      }
      points(x.medians, y=(0.5:dim(x)[1])+0.15, pch=18, cex = min(35/nrow(x),4))
      points(y.medians, y=(0.5:dim(y)[1])-0.15, pch=18, cex = min(35/nrow(x),4))
      mtext('Quantiles', 3, line=0, at=1, adj = 1.675, padj = 0.45, las=1, cex=0.7)
      ### create.tints.rgb is in utils.r
      red.tints  <- create.tints.rgb(col2rgb(col[2])/255, nr.tints=5, tint.steps = 1/5)
      blue.tints <- create.tints.rgb(col2rgb(col[1])/255, nr.tints=5, tint.steps = 1/5)
      legend(-1.75, nrow(x), legend = c("80%","60%","40%","20%","median","","","","",""),
             bty='n', cex=1, fill=c(
               rgb(matrix(red.tints[,5], ncol = 3)),
               rgb(matrix(red.tints[,3], ncol = 3)),
               rgb(matrix(red.tints[,2], ncol = 3)),
               rgb(matrix(red.tints[,1], ncol = 3)),
               0,
               rgb(matrix(blue.tints[,5], ncol = 3)),
               rgb(matrix(blue.tints[,3], ncol = 3)),
               rgb(matrix(blue.tints[,2], ncol = 3)),
               rgb(matrix(blue.tints[,1], ncol = 3)),
               0),
             lty    <- c(0,0,0,0,0,0,0,0,0,0),
             lwd    <- c(1.3,1.3,1.3,1.3,2,1.3,1.3,1.3,1.3,1.3), ncol = 2,
             border <- c("black", "black","black","black","white","black","black","black","black","white"))
      legend(-1.675, nrow(x), legend = c("","","","",""),
             bty='n', lty = c(0,0,0,0,0),
             # cap legend size for diamond (should look symmetric to other symbols)
             pch = 18, cex = 1, pt.cex = c(0,0,0,0, min(35/nrow(x), 2.25)))

      label.plot.horizontal(x, y, rownames(feat)[idx], x.suff=paste(' (', label$p.lab, ')', sep=''),
                            y.suff=paste(' (', label$n.lab, ')', sep=''), outer.diff = 1, inner.diff.x = -0.3, inner.diff.y = -0.6)
    }
    else if (plot.type == "bean"){
      par(mar=c(5.1, 25.1, 4.1, 0))
      bean.data <- data.frame()
      for (i in 1:nrow(x)){
        temp      <- as.data.frame(rbind(cbind(x[i, ], rep(paste(label$n.lab, rownames(x)[i]), length(x[i, ]))),
                                    cbind(y[i, ], rep(paste(label$p.lab, rownames(x)[i]), length((y[i, ]))))))
        temp[,1]  <- as.numeric(as.character(temp[,1]))
        bean.data <- rbind(bean.data, temp)
        if (i == nrow(x)){
          plot(NULL, xlab='', ylab='',xaxs='i', yaxs='i', axes=FALSE,
               xlim = c(as.integer(min(x))-1.5,as.integer(max(x))+1), ylim=c(0.45, length(idx)+0.6), type='n')
          beanplot(bean.data[, 1] ~ bean.data[, ncol(bean.data)], side = "both", bw="nrd0", col = list(col[1],
                                                                                                       col[2]), horizontal = TRUE, names = c(""), show.names = FALSE, beanlines = "median", maxstripline = 0.2, what = c(FALSE,TRUE,TRUE,FALSE),
                   axes = FALSE, add = TRUE )
          mn    <- as.integer(c(min(bean.data[,1])-1.5))
          mx    <- as.integer(c(max(bean.data[,1])+1))
          ticks <- mn:mx
          for (v in ticks) {
            abline(v=v, lty=3, col='lightgrey')
          }
          tick.labels <- formatC(10^ticks, format='E', digits=0)
          axis(side=1, at=ticks, labels=tick.labels, cex.axis=0.7)
          label.plot.horizontal(x, y, rownames(feat)[idx], x.suff=paste(' (', label$p.lab, ')', sep=''),
                                y.suff=paste(' (', label$n.lab, ')', sep=''), outer.diff = 1, inner.diff.x = 0.15, inner.diff.y = -0.15)
        }
      }
    }
    else {
      print("plot type has not been specified properly; continue with quantileplot")
      plot.type <- "quantile.box"
      par(mar=c(5.1, 25.1, 4.1, 0))
      plot.data.range(x, y, rownames(feat)[idx], x.col=col[2], y.col=col[1])
      label.plot.horizontal(x, y, rownames(feat)[idx], x.suff=paste(' (', label$p.lab, ')', sep=''),
                            y.suff=paste(' (', label$n.lab, ')', sep=''), outer.diff = 1, inner.diff.x = 0.15, inner.diff.y = -0.15)
    }

    p.val.annot <- formatC(p.adj[idx], format='E', digits=2)
    if (sum(p.adj < alpha, na.rm=TRUE) <= max.show) {
      title(main='Differentially abundant features', xlab='Abundance (log10-scale)')
    } else {
      title(main=paste('Differentially abundant features\ntruncated to the top', max.show),
            xlab='Abundance (log10-scale)')
    }
    par(mar=c(5.1,0,4.1, 0))
    for (i in 1:length(p.val.annot)) {
      if (plot.type == 'box'){
        mtext(p.val.annot[i], 4, line=2, at=(2*i)-0.5, las=1, cex=min(0.7, 1-(length(idx)/100)))
      }
      else if (plot.type == "quantile.rect"){
        mtext(p.val.annot[i], 4, line=2, at=i-0.5, las=1, cex=min(0.7, 1-(length(idx)/100)))
      }
      else {
        mtext(p.val.annot[i], 4, line=2, at=i, las=1, cex=min(0.7, 1-(length(idx)/100)))
      }
    }
    plot(NULL, xlab='', ylab='',xaxs='i', yaxs='i', axes=FALSE,
         type='n', xlim=c(0,10), ylim=c(0,length(p.val.annot)+0.5))
    title(main='Adj. p-value')

    # plot fold changes
    par(mar=c(5.1, 2.1, 4.1, 2.1))
    bcol  <- ifelse(fc[idx] > 0, col[2], col[1])
    mn    <- floor(min(fc[idx]))
    mx    <- ceiling(max(fc[idx]))
    mx    <- max(abs(mn), abs(mx))
    if (!is.finite(mx)) {
      mx    <- 10
    }
    mn    <- -mx
    plot(NULL, xlab='', ylab='', xaxs='i', yaxs='i', axes=FALSE,
         xlim=c(mn, mx), ylim=c(0.2, length(idx)+0.2), type='n')
    barplot(fc[idx], horiz=TRUE, width=0.6, space=2/3, col=bcol, axes=FALSE, add=TRUE)

    ticks <- mn:mx
    for (v in ticks) {
      abline(v=v, lty=3, col='lightgrey')
    }
    tick.labels <- formatC(10^ticks, format='E', digits=0)
    axis(side=1, at=ticks, labels=tick.labels, cex.axis=0.7)
    title(main='Fold change', xlab='FC (log10-scale)')

    # plot single-feature AUCs
    par(mar=c(5.1, 1.1, 4.1, 3.1))
    plot(NULL, xlab='', ylab='',xaxs='i', yaxs='i', axes=FALSE,
         xlim=c(0.5,1), ylim=c(0.5, length(idx)+0.5), type='n')
    ticks       <- seq(0.5, 1.0, length.out=6)
    for (v in ticks) {
      abline(v=v, lty=3, col='lightgrey')
    }
    for (b in 1:length(idx)) {
      i <- idx[b]
      points(aucs[i], b, pch=18, col=bcol[b])
      points(aucs[i], b, pch=5, col='black', cex=0.9)
    }
    axis(side=1, at=ticks, cex.axis=0.7)
    title(main='Feature AUCs', xlab='AU-ROC')
    # close pdf device
    tmp <- dev.off()
  }

}

#' @title Add metadata as predictors
#' @description This function adds metadata to the feature matrix to be later used as predictors
#' @param feat features object
#' @param meta metadata object
#' @param pred.names vector of names of the metavariables to be added to the feature matrix as predictors
#' @param std.meta boolean, should added metadata features be standardized?
#' @keywords SIAMCAT add.meta.pred
#' @export
#' @return features object with added metadata
add.meta.pred <- function(feat, meta, pred.names, std.meta){
  ### add metadata as predictors to the feature matrix
  cnt <- 0
  if (pred.names != '' && pred.names != 'NULL') {
    for (p in pred.names) {
      if(!p%in%colnames(meta)) stop("There is no meta variable called ",p,"\n")
      idx <- which(colnames(meta) == p)
      stopifnot(length(idx) == 1)
      cat('adding ', p, '\n', sep='')
      m = meta[,idx]
      if (!all(is.finite(m))) {
        na.cnt = sum(!is.finite(m))
        cat('filling in', na.cnt, 'missing values by mean imputation\n')
        mn = mean(m, na.rm=TRUE)
        m[!is.finite(m)] = mn
      }
      if (std.meta) {
        cat('standardize metadata feature', p, '\n')
        m.mean = mean(m, na.rm = TRUE)
        m.sd = sd(m, na.rm = TRUE)
        stopifnot(!m.sd == 0)
        m = (m - m.mean)/m.sd
      }
      feat = rbind(feat, m)
      rownames(feat)[nrow(feat)] = paste('META-', toupper(p), sep='')
      cnt = cnt + 1
    }
  }
  stopifnot(all(!is.na(feat)))
  cat('added', cnt, 'meta-variables as predictors to the feature matrix\n')
  invisible(feat)
}

#' @title Model training
#' @description This function trains the a machine learning model on the training data, using a \code{num.folds}-fold internal cross-validation scheme to find the optimal hyper-parameters of the model.
#' @param feat features object
#' @param label label object
#' @param fn.train.sample file containing the training samples
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
plm.trainer <- function(feat, label, fn.train.sample, num.folds=5, stratify, modsel.crit, min.nonzero.coeff){
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
  print(fn.train.sample)
  if (!is.null(fn.train.sample)) {
    num.runs      <- 0
    con           <- file(fn.train.sample, 'r')
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

  } else {
    # train on whole data set
    fold.name[[1]]    <- 'whole data set'
    fold.exm.idx[[1]] <- names(label$label)
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
    foldid            <- rep(0, length(train.label$label))
    perm              <- sample(1:length(train.label$label), length(train.label$label)) / length(train.label$label)
    for (f in num.folds:1) {
      foldid[perm <= f/num.folds] = f
    }
    train.label.exp   <- sample(train.label$label)
    foldid            <- assign.fold(label = train.label.exp, num.folds, stratified = stratify, foldid = foldid)
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

#' @title Prediction on the test set
#' @description This function takes the test set instances and the model trained by \link{plm.trainer} in order to predict the classes.
#' @param feat features object
#' @param label label object
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
plm.predictor <- function(feat, label, model, model.mat, hyperpars, model.type){
  # TODO hyperpars is not used at the moment, as far as i see
  # TODO 2: instead of feat and label containing the test sample indices, provide all features and the list from data.splitter
  feat         <- t(feat)
  label.fac                  <- factor(label$label, levels=c(label$negative.lab, label$positive.lab))

  data                       <- cbind(feat,label.fac[rownames(feat)])
  data                       <- as.data.frame(data)
  data[,ncol(data)]          <- as.factor(data[,ncol(data)])
  colnames(data)             <- paste0("Sample_",1:ncol(data))
  colnames(data)[ncol(data)] <- "cancer"


  ### subselect test examples as specified in fn.test.sample (if given)
  fold.name = list()
  fold.exm.idx = list()
  if (!is.null(fn.test.sample)) {
    con = file(fn.test.sample, 'r')
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
    stopifnot(all(paste('M', unlist(fold.name), sep='_') == colnames(model$W)))
  } else {
    # apply each LASSO model on whole data set when only single test set is given
    for (r in 1:num.runs) {
      fold.name[[r]] = paste('whole data set predicted by model', r)
      fold.exm.idx[[r]] = rownames(feat)
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
    m = model
    m$W = m$W[,r]

    # subselect test examples
    test.feat = feat[fold.exm.idx[[r]],,drop=FALSE]

    pdata    <- predict.plm(test.feat, model, model.type, opt.hp, data = data, subset=fold.exm.idx[[r]])
    p        <- label$negative.lab+abs(label$positive.lab-label$negative.lab)*pdata$data[,4]
    names(p) <- rownames(pdata$data)

    pred     <- c(pred, p)
    fold.pred.idx[[r]] = (length(pred)-length(p)+1):length(pred)
  }

  cat('\nTotal number of predictions made:', length(pred), '\n')

  if (!is.null(fn.test.label)) {
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
      if (!is.null(fn.test.label)) {
        ref.names = names(label$label)
      } else {
        ref.names = unique(names(pred))
      }
    } else {
      r.idx = as.numeric(sapply(strsplit(fold.name, 'rep'), '[[', 2))
      runs = sort(unique(r.idx))
      stopifnot(all(runs == 1:length(runs)))
      if (!is.null(fn.test.label)) {
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
      if (!is.null(fn.test.label)) {
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
