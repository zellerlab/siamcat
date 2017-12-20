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
#'  \item \code{$par} <- parameters utilized in the normalization;
#'  \item \code{$feat} <- normalized features
#'}
normalize.feat <- function(feat, norm.method = c("rank.unit", "rank.std", "log.std", "log.unit"), log.n0 = 10^-8, sd.min.q = 0.1, n.p = 2,
                           n.sample = FALSE, n.feature = TRUE, n.global = FALSE) {
  ### remove features with missing values
  # TODO there may be better ways of dealing with NA features
  num.orig.feat <- nrow(feat)
  keep.idx <- rowSums(is.na(feat) == 0)
  if (any(!keep.idx)) {
    feat <- feat[keep.idx,]
    cat('Removed ', nrow(feat)-num.orig.feat, ' features with missing values (retaining ', nrow(feat),' )\n', sep='')
  }

  ### keep track of normalization parameters
  par <- list()
  par$norm.method <- norm.method
  par$log.n0 <- log.n0
  par$n.p <- n.p
  par$n.sample <- n.sample
  par$n.feature <- n.feature
  par$n.global <- n.global
  par$retained.feat <- rownames(feat)

  ### apply normalization
  if (norm.method == 'rank.unit') {

    feat <- apply(feat, 2, rank, ties.method='average')

    stopifnot(!any(is.na(feat)))

    feat <- apply(feat, 2, FUN=function(x){x/sqrt(sum(x^2))})

  } else if (norm.method == 'rank.std') {

    feat <- apply(feat, 2, rank, ties.method='average')

    m <- apply(feat, 1, mean)
    s <- apply(feat, 1, sd)
    q <- quantile(s, sd.min.q, names=FALSE)
    stopifnot(q > 0)

    feat <- t(apply(feat, 1, FUN=function(x, q){(x - mean(x))/(sd(x) + q)}, q=q))

    par$feat.mean <- m
    par$feat.adj.sd <- s + q
    stopifnot(!any(is.na(feat)))

  } else if (norm.method == 'log.std') {
    feat <- log10(feat + log.n0)
    m <- apply(feat, 1, mean)
    s <- apply(feat, 1, sd)
    q <- quantile(s, sd.min.q, names=FALSE)
    stopifnot(q > 0)

    feat <- t(apply(feat, 1, FUN=function(x, q){(x - mean(x))/(sd(x) + q)}, q=q))

    par$feat.mean <- m
    par$feat.adj.sd <- s + q
    stopifnot(!any(is.na(feat)))

  } else if (norm.method == 'log.unit') {
    cat('Feature sparsity before normalization: ', 100*mean(feat==0), '%\n', sep='')
    feat <- log10(feat + log.n0)

    if (n.p == 1) {
      if (n.feature) {
        feat.norm.denom <- rowSums(feat)
        feat <- apply(feat, 1, FUN=function(x){x/sum(x)})
        par$feat.norm.denom <- feat.norm.denom
      }
      if (n.sample) {
        feat <- apply(feat, 2, FUN=function(x){x/sum(x)})
      }

    } else if (n.p == 2) {
      if (n.feature) {
        feat.norm.denom <- sqrt(rowSums(feat^2))
        feat <- apply(feat, 1, FUN=function(x){x/sqrt(sum(x^2))})
        par$feat.norm.denom <- feat.norm.denom
      }
      if (n.sample) {
        feat <- apply(feat, 2, FUN=function(x){x/sqrt(sum(x^2))})
      }

    } else {
      stop('unknown norm!')
    }

    if (!n.feature && !n.sample && n.global) {
      global.norm.denom <- max(feat)
      feat <- feat / global.norm.denom
      par$global.norm.denom <- global.norm.denom
    }

    cat('Feature sparsity after normalization: ', 100*mean(feat==0), '%\n', sep='')
    stopifnot(!any(is.na(feat)))
  } else {
    stop('unrecognized norm.method, exiting!\n')
  }
  return(list("par"= par, "feat" =feat))
}
