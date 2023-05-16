#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Perform feature normalization
#'
#' @description This function performs feature normalization according to 
#' user-specified parameters.
#'
#' @usage normalize.features(siamcat, norm.method = c("rank.unit", "rank.std",
#' "log.std", "log.unit", "log.clr", "std", "pass"), 
#' norm.param = list(log.n0 = 1e-06, sd.min.q = 0.1, n.p = 2, norm.margin = 1),
#' feature.type='filtered', verbose = 1)
#'
#' @param siamcat an object of class \link{siamcat-class}
#'
#' @param norm.method string, normalization method, can be one of these:
#' \code{c('rank.unit', 'rank.std', 'log.std', 'log.unit', 'log.clr','std', 
#' 'pass')}
#'
#' @param norm.param list, specifying the parameters of the different 
#' normalization methods, see Details for more information
#'
#' @param feature.type string, on which type of features should the function
#' work? Can be either \code{"original"}, \code{"filtered"}, or
#' \code{"normalized"}. Please only change this paramter if you know what
#' you are doing!
#'
#' @param verbose integer, control output: \code{0} for no output at all,
#' \code{1} for only information about progress and success, \code{2} for
#' normal level of information and \code{3} for full debug information,
#' defaults to \code{1}
#'
#' @section Implemented methods:
#' There are seven different normalization methods available, which 
#' might need additional parameters, which are passed via the \code{norm.param}
#' list:
#' \itemize{
#' \item \code{'rank.unit'} - converts features to ranks and normalizes 
#' each column (=sample) by the square root of the sum of ranks
#' This method does not require additional parameters.
#' \item \code{'rank.std'} - converts features to ranks and applies z-score 
#' standardization. 
#' This method requires \code{sd.min.q} (minimum quantile of the standard 
#' deviation to be added to all features in order to avoid underestimation of 
#' standard deviation) as additional parameter.
#' \item \code{'log.clr'} - centered log-ratio transformation.
#' This methods requires a pseudocount (\code{log.n0}) before 
#' log-transformation.
#' \item \code{'log.std'} - log-transforms features and applies z-score 
#' standardization.
#' This method requires both a pseudocount (\code{log.n0}) and \code{sd.min.q}
#' \item \code{'log.unit'} - log-transforms features and normalizes by 
#' features or samples with different norms.
#' This method requires a pseudocount (\code{log.n0}) and then additionally the
#' parameters \code{norm.maring} (margin over which to normalize, similarly to 
#' the \code{apply}-syntax: Allowed values are \code{1} for normalization 
#' over features, \code{2} over samples, and \code{3} for normalization 
#' by the global maximum) and the parameter \code{n.p} (vector norm to be 
#' used, can be either \code{1} for \code{x/sum(x)} or \code{2} for 
#' \code{x/sqrt(sum(x^2))}).
#' \item \code{'std'} - z-score standardization without any other 
#' transformation
#' This method only requires the \code{sd.min.q} parameter
#' \item \code{'pass'} - pass-through normalization will not change 
#' the features}
#'
#' @section Frozen normalization:
#' The function additionally allows to perform a frozen normalization on a
#' different dataset. After normalizing the first dataset, the \code{norm_feat}
#' slot in the SIAMCAT object contains all parameters of the normalization,
#' which you can access via the \link{norm_params} accessor.
#'
#' In order to perform a frozen normalization of a new dataset, you can run the
#' function supplying the normalization parameters as argument to
#' \code{norm.param}:
#' \code{norm.param=norm_params(siamcat_reference)}. See also the example below.
#'
#' @keywords SIAMCAT normalize.features
#'
#' @export
#'
#' @encoding UTF-8
#'
#' @return an object of class \link{siamcat-class} with normalized features
#'
#' @examples
#' # Example data
#' data(siamcat_example)
#'
#' # Simple example
#' siamcat_norm <- normalize.features(siamcat_example,
#'     norm.method='rank.unit')
#'
#' # log.unit example
#' siamcat_norm <- normalize.features(siamcat_example,
#'     norm.method='log.unit',
#'     norm.param=list(log.n0=1e-05, n.p=1, norm.margin=1))
#'
#' # log.std example
#' siamcat_norm <- normalize.features(siamcat_example,
#'     norm.method='log.std',
#'     norm.param=list(log.n0=1e-05, sd.min.q=.1))
#'
#' # Frozen normalization
#' # normalize the object siamcat with the same parameters as used in 
#' # siamcat_reference
#' # 
#' # this is not run
#' # siamcat_norm <- normalize.features(siamcat,
#' #   norm.param=norm_params(siamcat_reference))
normalize.features <- function(siamcat,
    norm.method = c("rank.unit", "rank.std",
        "log.std", "log.unit",
        "log.clr", "std", "pass"),
    norm.param = list(
        log.n0 = 1e-06,
        sd.min.q = 0.1,
        n.p = 2,
        norm.margin = 1
    ),
    feature.type='filtered',
    verbose = 1) {
    if (verbose > 1)
        message("+ starting normalize.features")
    s.time <- proc.time()[3]

    if (!feature.type %in% c('original', 'filtered', 'normalized')){
        stop("Unrecognised feature type, exiting...\n")
    }

    # get the right features
    if (feature.type == 'original'){
        feat <- get.orig_feat.matrix(siamcat)
        if (verbose > 1) message('+ normalizing original features')
    } else if (feature.type == 'filtered'){
        if (is.null(filt_feat(siamcat, verbose=0))){
            stop('Features have not yet been filtered, exiting...\n')
        }
        feat <- get.filt_feat.matrix(siamcat)
    } else if (feature.type == 'normalized'){
        if (is.null(norm_feat(siamcat, verbose=0))){
            stop('Features have not yet been normalized, exiting...\n')
        }
        msg <- paste('Normalizing features that have already been normalized!',
            '\nPlease note that some functionalities may not work properly')
        warning(msg)
        feat <- get.norm_feat.matrix(siamcat)
    }

    if (is.null(norm.param$norm.method)) {
        # de novo normalization
        if (verbose > 1){
            msg <- paste( "+++ performing de novo normalization using the ",
                norm.method, " method")
            message(msg)
        }
        ### remove features with missing values
        keep.idx <- rowSums(is.na(feat)) == 0
        if (any(!keep.idx)) {
            feat.red.na <- feat[keep.idx,]
            if (verbose > 1) {
                msg <- paste0("+++ removed ", nrow(feat.red.na) - nrow(feat),
                        " features with missing values (retaining ",
                        nrow(feat.red.na), ")")
                message(msg)
            }
        } else {
            feat.red.na <- feat
        }
        ### remove features with zero sd
        keep.idx.sd <- (rowSds(feat.red.na) == 0)
        if (any(keep.idx.sd)) {
            feat.red <- feat.red.na[!keep.idx.sd,]
            if (verbose > 1) {
                msg <- paste0("+++ removed ", nrow(feat.red.na)-nrow(feat.red),
                        " features with no variation across samples
                        (retaining ", nrow(feat.red), ")")
                message(msg)
            }
        } else {
            feat.red <- feat.red.na
        }
        ## check if one of the allowed norm.methods have been supplied
        if (length(norm.method) != 1){
            stop('Please supply only a single normalization method! Exiting...')
        }
        if (!norm.method %in%
            c("rank.unit", "rank.std", "log.std", "log.unit",
                "log.clr", "std", "pass")) {
            stop("Unknown normalization method! Exiting...")
        }

        ### keep track of normalization parameters
        par <- list()
        par$norm.method <- norm.method
        par$retained.feat <- rownames(feat.red)
        if (verbose > 2)
            message("+++ checking if parameters are compatible with each other")
        # check if the right set of normalization parameters have been
        # supplied for the chosen norm.method
        if (norm.method %in% c("rank.std", "std") &&
                is.null(norm.param$sd.min.q)) {
            stop(
                "The rank.std method requires the parameter sd.min.q, which is
                not supplied. Exiting ..."
            )
        }
        if (startsWith(norm.method, "log")){
            if (any(feat.red < 0)){
                stop('Can not perform log-transform on ',
                    'negative data. Exiting...')
            }
            if (is.null(norm.param$log.n0)){
                warning(
                    "Pseudo-count before log-transformation not supplied!
                    Estimating it as 5% percentile...\n"
                )
                norm.param$log.n0 <-
                    quantile(feat.red[feat.red != 0], 0.05)
            }
        }
        if (norm.method == "log.std" &&
                (is.null(norm.param$sd.min.q))) {
            msg <- paste0("The log.std method requires the parameter",
                " sd.min.q, which is not supplied. Exiting ...")
            stop(msg)
        }
        if (norm.method == "log.unit" &&
                (is.null(norm.param$n.p) ||
                        is.null(norm.param$norm.margin))) {
            msg <- paste0("The log.std method requires the parameters n.p and", 
                "norm.margin, which are not supplied. Exiting ...")
            stop(msg)
        }

        ### keep track of normalization parameters
        par$log.n0 <- norm.param$log.n0
        par$n.p <- norm.param$n.p
        par$norm.margin <- norm.param$norm.margin

        if (verbose > 1){
            msg <- paste0("+ feature sparsity before normalization: ",
                formatC(100 * mean(feat.red == 0), digits = 4), "%")
            message(msg)
        }
        # normalization
        if (verbose > 2)
            message("+++ performing normalization")
        if (norm.method == "rank.unit") {
            feat.rank <- colRanks(feat.red,
                preserveShape = TRUE,
                ties.method = 'average')
            stopifnot(!any(is.na(feat.rank)))
            feat.norm <- t(t(feat.rank) / sqrt(colSums(feat.rank ^ 2)))
            dimnames(feat.norm) <- dimnames(feat.red)
        } else if (norm.method == "log.clr") {
            gm <- exp(colMeans(log(feat.red + norm.param$log.n0)))
            feat.norm <- log(t(t((feat.red + norm.param$log.n0)) / gm))
        } else if (norm.method == "rank.std") {
            feat.rank <- colRanks(feat.red,
                preserveShape = TRUE,
                ties.method = 'average')
            dimnames(feat.rank) <- dimnames(feat.red)
            m <- rowMeans(feat.rank)
            s <- rowSds(feat.rank)
            q <- quantile(s, norm.param$sd.min.q, names = FALSE)
            stopifnot(q > 0)
            feat.norm <- (feat.rank - m) / (s + q)
            par$feat.mean <- m
            par$feat.adj.sd <- s + q
        } else if (norm.method == "std") {
            m <- rowMeans(feat.red)
            s <- rowSds(feat.red)
            q <- quantile(s, norm.param$sd.min.q, names = FALSE)
            stopifnot(q > 0)
            feat.norm <- (feat.red - m) / (s + q)
            names(m) <- rownames(feat.red)
            names(s) <- rownames(feat.red)
            par$feat.mean <- m
            par$feat.adj.sd <- s + q
        } else if (norm.method == "log.std") {
            feat.log <- log10(feat.red + norm.param$log.n0)
            m <- rowMeans(feat.log)
            s <- rowSds(feat.log)
            q <- quantile(s, norm.param$sd.min.q, names = FALSE)
            stopifnot(q > 0)
            feat.norm <- (feat.log - m) / (s + q)
            names(m) <- rownames(feat.log)
            names(s) <- rownames(feat.log)
            par$feat.mean <- m
            par$feat.adj.sd <- s + q
        } else if (norm.method == "log.unit") {
            feat.log <- log10(feat.red + norm.param$log.n0)
            if (norm.param$n.p == 1) {
                norm.fun <- function(x) {
                    x / sum(x)
                }
                par$norm.fun <- norm.fun
                par$feat.norm.denom <- rowSums(feat.log)
            } else if (norm.param$n.p == 2) {
                norm.fun <- function(x) {
                    x / sqrt(sum(x ^ 2))
                }
                par$norm.fun <- norm.fun
                par$feat.norm.denom <- sqrt(rowSums(feat.log ^ 2))
            } else {
                stop("Unknown vector norm, must be either 1 or 2. Exiting...")
            }
            if (norm.param$norm.margin == 1) {
                feat.norm <- t(apply(feat.log, 1, FUN = norm.fun))
            } else if (norm.param$norm.margin == 2) {
                feat.norm <- apply(feat.log, 2, FUN = norm.fun)
            } else if (norm.param$norm.margin == 3) {
                par$global.norm <- max(feat.log)
                feat.norm <- feat.log / max(feat.log)
            } else {
                stop(
                    "Unknown margin for normalization, must be either
                    1 (for features), 2 (for samples), or 3 (global).
                    Exiting..."
                )
            }
        } else if (norm.method == "pass"){
            feat.norm <- feat.red
        }
        if (verbose > 1){
            msg <- paste("+++ feature sparsity after normalization: ",
                formatC(100 * mean(feat.norm == 0), digits = 4), "%")
            message(msg)
        }
        stopifnot(!any(is.na(feat.norm)))
    } else {
        # frozen normalization
        if (verbose > 1){
            msg <- paste("+ performing frozen ", norm.param$norm.method,
                " normalization using the supplied parameters")
            message(msg)
        }
        # check if all retained.feat in norm.params are also found in features
        stopifnot(all(norm.param$retained.feat %in% row.names(feat)))
        feat.red <- feat[norm.param$retained.feat,]

        if (verbose > 1){
            msg <- paste0("+ feature sparsity before normalization: ",
                formatC(100 * mean(feat.red == 0), digits = 4), "%")
            message(msg)
        }

        # normalization
        if (verbose > 2)
            message("+++ performing normalization")
        if (norm.param$norm.method == "rank.unit") {
            feat.rank <- colRanks(feat.red,
                preserveShape = TRUE,
                ties.method = 'average')
            stopifnot(!any(is.na(feat.rank)))
            feat.norm <- t(t(feat.rank) / sqrt(colSums(feat.rank ^ 2)))
            dimnames(feat.norm) <- dimnames(feat.red)

        } else if (norm.param$norm.method == "log.clr") {
            gm <- exp(colMeans(log(feat.red + norm.param$log.n0)))
            feat.norm <- log(t(t((feat.red + norm.param$log.n0)) / gm))

        } else if (norm.param$norm.method == "rank.std") {
            stopifnot(
                !is.null(norm.param$feat.mean) &&
                    !is.null(norm.param$feat.adj.sd) &&
                    all(
                        names(norm.param$feat.mean) ==
                            row.names(feat.red)
                    ) && all(
                        names(norm.param$feat.adj.s) ==
                            row.names(feat.red)
                    )
            )
            feat.rank <- colRanks(feat.red,
                preserveShape = TRUE,
                ties.method = 'average')
            dimnames(feat.rank) <- dimnames(feat.red)
            feat.norm <- (feat.rank - norm.param$feat.mean) /
                norm.param$feat.adj.s

        } else if (norm.param$norm.method == "std") {
            stopifnot(
                !is.null(norm.param$feat.mean) &&
                    !is.null(norm.param$feat.adj.sd) &&
                    all(
                        names(norm.param$feat.mean) == row.names(feat.red)
                    ) &&
                    all(
                        names(norm.param$feat.adj.s) == row.names(feat.red)
                    )
            )
            feat.norm <- (feat.red - norm.param$feat.mean) /
                norm.param$feat.adj.sd
        } else if (norm.param$norm.method == "log.std") {
            stopifnot(
                !is.null(norm.param$log.n0) &&
                    !is.null(norm.param$feat.mean) &&
                    !is.null(norm.param$feat.adj.sd) &&
                    all(
                        names(norm.param$feat.mean) == row.names(feat.red)
                    ) &&
                    all(
                        names(norm.param$feat.adj.s) == row.names(feat.red)
                    )
            )
            feat.log <- log10(feat.red + norm.param$log.n0)
            feat.norm <- (feat.log - norm.param$feat.mean) /
                norm.param$feat.adj.sd

        } else if (norm.param$norm.method == "log.unit") {
            stopifnot(
                !is.null(norm.param$log.n0) &&
                    !is.null(norm.param$norm.margin) &&
                    norm.param$norm.margin %in%
                    c(1, 2, 3) && ((
                        norm.param$norm.margin == 1 &&
                            !is.null(norm.param$feat.norm.denom) &&
                            all(
                                names(norm.param$feat.norm.denom) ==
                                    row.names(feat.red)
                            )
                    ) || (
                        norm.param$norm.margin == 2 &&
                            !is.null(norm.param$norm.fun)
                    ) || (
                        norm.param$norm.margin ==
                            3 && !is.null(norm.param$global.norm)
                    )
                    )
            )
            feat.log <- log10(feat.red + norm.param$log.n0)
            if (norm.param$norm.margin == 1) {
                feat.norm <- feat.log / norm.param$feat.norm.denom
            } else if (norm.param$norm.margin == 2) {
                feat.norm <- apply(feat.log, 2, FUN = norm.param$norm.fun)
            } else if (norm.param$norm.margin == 3) {
                feat.norm <- feat.log / norm.param$global.norm
            }
        } else if (norm.param$norm.method == "pass"){
            feat.norm <- feat.red
        }
        if (verbose > 1){
            msg <- paste0("+ feature sparsity after normalization: ",
                formatC(100 * mean(feat.norm == 0), digits = 4), "%")
            message(msg)
        }
        stopifnot(!any(is.na(feat.norm)))
        par <- norm.param
    }

    norm_feat(siamcat) <- list(
        norm.feat=feat.norm,
        norm.param=par)

    e.time <- proc.time()[3]
    if (verbose > 1){
        msg <- paste("+ finished normalize.features in", 
            formatC(e.time - s.time, digits = 3), "s")
        message(msg)
    }
    if (verbose == 1)
        message("Features normalized successfully.")
    return(siamcat)
}
