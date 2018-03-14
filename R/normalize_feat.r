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

#' @title Perform feature normalization
#' @description This function performs feature normalization according to user-
#'  specified parameters.
#' @param siamcat an object of class \link{siamcat}
#' @param norm.method string, normalization method, can be one of these:
#'  '\code{c("rank.unit", "rank.std", "log.std", "log.unit", "clr")}
#' @param norm.param list, specifying the parameters of the different
#'  normalization methods, see details for more information
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'        for only information about progress and success, \code{2} for normal
#'        level of information and \code{3} for full debug information, defaults to \code{1}
#' @details There are five different normalization methods available:
#' \itemize{
#'  \item \code{"rank.unit"} converts features to ranks and normalizes each
#'        column (=sample) by the square root of the sum of ranks
#'  \item \code{"rank.std"} converts features to ranks and applies z-score
#'        standardization
#'  \item \code{"clr"} centered log-ratio transformation (with the addition of
#'        pseudocounts)
#'  \item \code{"log.std"} log-transforms features (after addition of
#'        pseudocounts) and applies z-score standardization
#'  \item \code{"log.unit"} log-transforms features (after addition of
#'        pseudocounts) and normalizes by features or samples with different norms
#' }
#'
#' The list entries in \code{"norm.param"} specify the normalzation parameters,
#' which are dependant on the normalization method of choice:
#' \itemize{
#'  \item \code{"rank.unit"} does not require any other parameters
#'  \item \code{"rank.std"} requires \code{sd.min.q}, quantile of the
#'        distribution of standard deviations of all features that will be added
#'        to the denominator during standardization in order to avoid
#'        underestimation of the standard deviation, defaults to 0.1
#'  \item \code{"clr"} requires \code{log.n0}, which is the pseudocount to be
#'        added before log-transformation, defaults to \code{NULL} leading to the
#'        estimation of \code{log.n0} from the data
#'  \item \code{"log.std"} requires both \code{log.n0} and \code{sd.min.q}, using
#'        the same default values
#'  \item \code{"log.unit"} requires next to \code{log.n0} also the parameters
#'        \code{n.p} and \code{norm.margin}. \code{n.p} specifies the vector norm
#'        to be used, can be either \code{1} for \code{x/sum(x)} or \code{2} for
#'        \code{x/sqrt(sum(x^2))}. The parameter \code{norm.margin} specifies the
#'        margin over which to normalize, similarly to the \code{apply}-syntax:
#'        Allowed values are \code{1} for normalization over features, \code{2}
#'        over samples, and \code{3} for normalization by the global maximum.
#'}
#'
#' The function allows to perform a frozen normalization on a different dataset.
#' After normalizing the first dataset, the output list \code{$par} contains all
#' parameters of the normalization. Supplying this list together with a new dataset
#' will normalize the second dataset in a comparable way to the first dataset (e.g.
#' by using the same mean for the features for z-score standardization)
#'
#' @keywords SIAMCAT normalize.features
#' @export
#' @return an object of class \link{siamcat}

normalize.features   <- function(siamcat, norm.method=c("rank.unit", "rank.std", "log.std", "log.unit", "clr"),
                           norm.param=list(log.n0=1e-08, sd.min.q=0.1, n.p=2, norm.margin=1), verbose=1) {

  if(verbose>1) cat("+ starting normalize.features\n")
  s.time <- proc.time()[3]
  feat <- siamcat@phyloseq@otu_table

  if (is.null(norm.param$norm.method)){
    # de novo normalization
    if(verbose>1) cat('+++ performing de novo normalization using the ', norm.method, ' method\n')
    ### remove features with missing values
    # TODO there may be better ways of dealing with NA features
    keep.idx <- rowSums(is.na(feat) == 0)
    if (any(!keep.idx)) {
      feat.red <- feat[keep.idx,]
      if(verbose>1) cat('+++ removed ', nrow(feat.red)-nrow(feat), ' features with missing values (retaining ', nrow(feat.red),' )\n', sep='')
    } else {
      feat.red <- feat
    }

    ## check if one of the allowed norm.methods have been supplied
    if (!norm.method %in% c("rank.unit", "rank.std", "log.std", "log.unit", "clr")){
      stop("Unknown normalization method! Exiting...")
    }

    ### keep track of normalization parameters
    par <- list()
    par$norm.method <- norm.method
    par$retained.feat <- rownames(feat.red)
    if(verbose>2) cat("+++ checking is parameters are compatible with each other\n")
    ## check if the right set of normalization parameters have been supplied for the chosen norm.method
    if (norm.method == 'rank.std' && is.null(norm.param$sd.min.q)){
      stop("The rank.std method requires the parameter sd.min.q, which is not supplied. Exiting ...")
    }
    if (norm.method != 'rank.std' && is.null(norm.param$log.n0)){
      warning("Pseudo-count before log-transformation not supplied! Estimating it as 5% percentile...\n")
      norm.param$log.n0 <- quantile(feat.red[feat.red!=0], 0.05)
    }
    if (norm.method == 'log.std' && (is.null(norm.param$sd.min.q))){
      stop("The log.std method requires the parameters sd.min.q, which is not supplied. Exiting ...")
    }
    if (norm.method == 'log.unit' && (is.null(norm.param$n.p) || is.null(norm.param$norm.margin))){
      stop("The log.std method requires the parameters n.p and norm.margin, which are not supplied. Exiting ...")
    }

    ### keep track of normalization parameters
    par$log.n0 <- norm.param$log.n0
    par$n.p <- norm.param$n.p
    par$norm.margin <- norm.param$norm.margin

    if(verbose>1) cat('+ feature sparsity before normalization: ', 100*mean(feat.red == 0), '%\n', sep='')
    # normalization
    if(verbose>2) cat("+++ performing normalization\n")
    if (norm.method == 'rank.unit'){
      feat.rank <- apply(feat.red, 2, rank, ties.method='average')
      stopifnot(!any(is.na(feat)))
      feat.norm <- apply(feat.rank, 2, FUN=function(x){x/sqrt(sum(x^2))})
    } else if (norm.method == 'clr'){
      feat.log <- feat.red + norm.param$log.n0
      gm <- apply(feat.log, 1, FUN=function(x){exp(mean(log(x)))})
      feat.norm <- t(apply(feat.log, 1, FUN=function(x){log(x/exp(mean(log(x))))}))
      par$geometric.mean <- gm
    } else if (norm.method == 'rank.std'){
      feat.rank <- apply(feat.red, 2, rank, ties.method='average')
      m <- apply(feat.rank, 1, mean)
      s <- apply(feat.rank, 1, sd)
      q <- quantile(s, norm.param$sd.min.q, names=FALSE)
      stopifnot(q > 0)
      feat.norm <- t(apply(feat.rank, 1, FUN=function(x, q){(x - mean(x))/(sd(x) + q)}, q=q))
      par$feat.mean <- m
      par$feat.adj.sd <- s + q
    } else if (norm.method == 'log.std'){
      feat.log <- log10(feat.red + norm.param$log.n0)
      m <- apply(feat.log, 1, mean)
      s <- apply(feat.log, 1, sd)
      q <- quantile(s, norm.param$sd.min.q, names=FALSE)
      stopifnot(q > 0)
      feat.norm <- t(apply(feat.log, 1, FUN=function(x, q){(x - mean(x))/(sd(x) + q)}, q=q))
      par$feat.mean <- m
      par$feat.adj.sd <- s + q
    } else if (norm.method == 'log.unit'){
      feat.log <- log10(feat.red + norm.param$log.n0)
      if (norm.param$n.p == 1){
        norm.fun <- function(x){x/sum(x)}
        par$norm.fun <- norm.fun
        par$feat.norm.denom <- rowSums(feat.log)
      } else if (norm.param$n.p == 2){
        norm.fun <- function(x){x/sqrt(sum(x^2))}
        par$norm.fun <- norm.fun
        par$feat.norm.denom <- sqrt(rowSums(feat.log^2))
      } else {
        stop('Unknown vector norm, must be either 1 or 2. Exiting...')
      }
      if (norm.param$norm.margin == 1){
        feat.norm <- t(apply(feat.log, 1, FUN=norm.fun))
      } else if (norm.param$norm.margin == 2){
        feat.norm <- apply(feat.log, 2, FUN=norm.fun)
      } else if (norm.param$norm.margin == 3){
        par$global.norm <- max(feat.log)
        feat.norm <- feat.log/max(feat.log)
      } else {
        stop('Unknown margin for normalization, must be either 1 (for features), 2 (for samples), or 3 (global). Exiting...')
      }
    }
    if(verbose>1) cat('+++ feature sparsity after normalization: ', 100*mean(feat.norm == 0), '%\n', sep='')
    stopifnot(!any(is.na(feat.norm)))
    siamcat@norm.param         <- par
    siamcat@norm.param$norm.method <- norm.method
  } else {
    # frozen normalization
    if(verbose>1) cat('+ performing frozen ', norm.param$norm.method, ' normalization using the supplied parameters\n')
    # check if all retained.feat in norm.params are also found in features
    stopifnot(all(norm.param$retained.feat %in% row.names(feat)))
    feat.red <- feat[norm.param$retained.feat,]

    if(verbose>1) cat('+ feature sparsity before normalization: ', 100*mean(feat.red == 0), '%\n', sep='')

    # normalization
    if(verbose>2) cat("+++ performing normalization\n")
    if (norm.param$norm.method == 'rank.unit'){
      feat.rank <- apply(feat.red, 2, rank, ties.method='average')
      feat.norm <- apply(feat.rank, 2, FUN=function(x){x/sqrt(sum(x^2))})
    } else if (norm.param$norm.method == 'clr'){
      stopifnot(!is.null(norm.param$log.n0) && !is.null(norm.param$geometric.mean) && all(names(norm.param$geometric.mean) == row.names(feat.red)))
      feat.log <- feat.red + norm.param$log.n0
      feat.norm <- log(feat.log/norm.param$geometric.mean)
    } else if (norm.param$norm.method == 'rank.std'){
      stopifnot(!is.null(norm.param$feat.mean) && !is.null(norm.param$feat.adj.sd) &&
                all(names(norm.param$feat.mean) == row.names(feat.red)) && all(names(norm.param$feat.adj.s) == row.names(feat.red)))
      feat.rank <- apply(feat.red, 2, rank, ties.method='average')
      feat.norm <- (feat.rank - norm.param$feat.mean)/norm.param$feat.adj.s
    } else if (norm.param$norm.method == 'log.std'){
      stopifnot(!is.null(norm.param$log.n0) && !is.null(norm.param$feat.mean) && !is.null(norm.param$feat.adj.sd) &&
                all(names(norm.param$feat.mean) == row.names(feat.red)) &&
                all(names(norm.param$feat.adj.s) == row.names(feat.red)))
      feat.log <- log10(feat.red + norm.param$log.n0)
      feat.norm <- (feat.log - norm.param$feat.mean)/norm.param$feat.adj.sd
    } else if (norm.param$norm.method == 'log.unit'){
      stopifnot(!is.null(norm.param$log.n0) && !is.null(norm.param$norm.margin) && norm.param$norm.margin %in% c(1,2,3) &&
                ((norm.param$norm.margin == 1 && !is.null(norm.param$feat.norm.denom) && all(names(norm.param$feat.norm.denom) == row.names(feat.red))) ||
                 (norm.param$norm.margin == 2 && !is.null(norm.param$norm.fun)) ||
                 (norm.param$norm.margin == 3 && !is.null(norm.param$global.norm))))
      feat.log <- log10(feat.red + norm.param$log.n0)
      if (norm.param$norm.margin == 1){
        feat.norm <- feat.log/norm.param$feat.norm.denom
      } else if (norm.param$norm.margin == 2){
        feat.norm <- apply(feat.log, 2, FUN=norm.param$norm.fun)
      } else if (norm.param$norm.margin == 3){
        feat.norm <- feat.log/norm.param$global.norm
      }
    }
    if(verbose>1) cat('+ feature sparsity after normalization: ', 100*mean(feat.norm == 0), '%\n', sep='')
    stopifnot(!any(is.na(feat.norm)))

  }
  siamcat@phyloseq@otu_table <- otu_table(feat.norm, taxa_are_rows = T)
  e.time <- proc.time()[3]
  if(verbose>1) cat("+ finished normalize.features in",e.time-s.time,"s\n")
  if(verbose==1)cat("Features normalized successfully.\n")
  return(siamcat)
}
