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
