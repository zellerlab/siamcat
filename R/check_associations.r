#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Calculate associations between features and labels
#'
#' @description This function computes different measures of association 
#' between features and the label and stores the results in the 
#' \code{association} slot of the SIAMCAT object
#'
#' @usage check.associations(siamcat, formula="feat~label", test='wilcoxon', 
#' alpha=0.05, mult.corr="fdr", log.n0=1e-06, pr.cutoff=1e-06, 
#' probs.fc=seq(.1, .9, .05), paired=NULL, feature.type='filtered', 
#' verbose = 1)
#'
#' @param siamcat object of class \link{siamcat-class}
#' 
#' @param formula string, formula used for testing, see Details for more 
#' information, defaults to \code{"feat~label"}
#' 
#' @param test string, statistical test used for the association testing, can
#' be either \code{'wilcoxon'} or \code{'lm'}, see Details for more 
#' information, defaults to \code{'wilcoxon'}
#'
#' @param alpha float, significance level, defaults to \code{0.05}
#'
#' @param mult.corr string, multiple hypothesis correction method, see 
#' \code{\link[stats]{p.adjust}}, defaults to \code{"fdr"}
#'
#' @param log.n0 float, pseudo-count to be added before log-transformation of 
#' the data, defaults to \code{1e-06}. Will be ignored if 
#' \code{feature.type} is \code{"normalized"}.
#'
#' @param pr.cutoff float, cutoff for the prevalence computation, defaults to 
#' \code{1e-06}
#'
#' @param probs.fc numeric vector, quantiles used to calculate the generalized
#' fold change between groups, see Details for more information,  
#' defaults to \code{seq(.1, .9, .05)}
#' 
#' @param paired character, column name of the meta-variable containing 
#' information for a paired test, defaults to \code{NULL}
#' 
#' @param feature.type string, on which type of features should the function
#' work? Can be either \code{c("original", "filtered", or "normalized")}.
#' Please only change this parameter if you know what you are doing!
#'
#' If \code{feature.type} is \code{"normalized"}, the normalized abundances
#' will not be log10-transformed.
#' 
#' @param verbose integer, control output: \code{0} for no output at all, 
#' \code{1} for only information about progress and success, \code{2} for 
#' normal level of information and \code{3} for full debug information, 
#' defaults to \code{1}
#'
#' @keywords SIAMCAT check.associations
#' 
#' @section Statistical testing:
#' The function uses the Wilcoxon test as default statistical test for binary 
#' classification problems. Alternatively, a simple linear model (as 
#' implemented in \link[stats]{lm}) can be used as well. For regression 
#' problems, the function defaults to the linear model.
#' 
#' @section Effect sizes:
#' The function calculates several measures for the effect size of the 
#' assocations between microbial features and the label. For binary 
#' classification problems, these associations are: \itemize{
#' \item AUROC (area under the Receiver Operating Characteristics curve) as a 
#' non-parametric measure of enrichment,
#' \item the generalized fold change (gFC), a pseudo-fold change which is 
#' calculated as geometric mean of the differences between quantiles across 
#' both groups,
#' \item prevalence shift (difference in prevalence between the two groups).}
#' For regression problems, the effect sizes are: \itemize{
#' \item Spearman correlation between the feature and the label.}
#' 
#' @section Confounder-corrected testing:
#' To correct for possible confounders while testing for association, the 
#' function uses linear mixed effect models as implemented in the 
#' \link{lmerTest} package. To do so, the test formula needs to be adjusted 
#' to include the confounder. For example, when correcting for the metadata 
#' information \code{Sex}, the formula would be: 
#' \code{'feat~label+(1|Sex)'} (see also the example below).
#' 
#' Please note that modifying the formula parameter in this function might
#' lead to unexpected results!
#'
#' @section Paired testing:
#' For paired testing, e.g. when the same patient has been sampled before and
#' after an intervention, the `paired` parameter can be supplied to the 
#' function. This indicated a column in the metadata table that holds the 
#' information about pairing.
#' 
#' @return object of class \link{siamcat-class} with the slot 
#' \code{associations} filled
#'
#' @export
#' 
#' @encoding UTF-8
#'
#' @examples
#' # Example data
#' data(siamcat_example)
#'
#' # Simple example
#' siamcat_example <- check.associations(siamcat_example)
#'
#'
#' # Confounder-corrected testing (corrected for Sex)
#' #
#' # this is not run during checks
#' # siamcat_example <- check.associations(siamcat_example, 
#' #     formula='feat~label+(1|Sex)', test='lm')
#' 
#' # Paired testing
#' #
#' # this is not run during checks
#' # siamcat_paired <- check.associations(siamcat_paired, 
#' #     paired='Individual_ID')
check.associations <- function(siamcat, formula="feat~label", test='wilcoxon',
    alpha=0.05, mult.corr="fdr", log.n0=1e-06, pr.cutoff=1e-06,
    probs.fc=seq(.1, .9, .05), paired=NULL, feature.type='filtered',
    verbose = 1) {

        if (verbose > 1)
            message("+ starting check.associations")
        s.time <- proc.time()[3]

        # check mult.corr
        if (!mult.corr %in% c("holm", "hochberg", "hommel", "bonferroni",
                                "BH", "BY", "fdr", "none")) {
            stop("Unknown multiple testing correction method: '", mult.corr,
                "'. Exiting!\n  Must of one of c('holm', 'hochberg', ",
                "'hommel', 'bonferroni', 'BH', 'BY', 'fdr', none')")
        }
        # check test
        if (!test %in% c('wilcoxon', 'lm')){
            stop("Unknown testing method: '", test,
                "'. Exiting!\n  Must of one of c('wilcoxon', 'lm')")
        }
        # check label
        label <- label(siamcat)
        if (label$type == 'TEST'){
            stop('Can not check assocations for a',
                ' SIAMCAT object with TEST label! Exiting...')
        }
        if (label$type=='CONTINUOUS' & test == 'wilcoxon'){
            warning("Cannot test a SIAMCAT object with regression label using",
                    " the Wilcoxon test. Changing test to 'lm'!")
            test <- 'lm'
        }
        meta <- meta(siamcat)
        # check feature type
        if (!feature.type %in% c('original', 'filtered', 'normalized')){
            stop("Unrecognised feature type, exiting...\n")
        }
        # get features
        if (feature.type == 'original'){
            feat <- get.orig_feat.matrix(siamcat)
        } else if (feature.type == 'filtered'){
            if (is.null(filt_feat(siamcat, verbose=0))){
                stop('Features have not yet been filtered, exiting...\n')
            }
            feat <- get.filt_feat.matrix(siamcat)
        } else if (feature.type == 'normalized'){
            if (is.null(norm_feat(siamcat, verbose=0))){
                stop('Features have not yet been normalized, exiting...\n')
            }
            feat <- get.norm_feat.matrix(siamcat)
        }
        if (any(is.na(feat))){
            stop('Features contain NAs. Exiting...')
        }
        if ((any(colSums(feat) > 1.01) | any(feat < -0.01)) &
            feature.type != 'normalized'){
            stop('This function expects compositional data. Exiting...')
        }
        # check paired
        if (!is.null(paired)){
            if (label$type!='BINARY'){
                stop('Paired testing is only supported for binary labels!')
            }
            if (!paired %in% colnames(meta)){
                msg <- paste0("Column with pairing information not present in",
                    " the metadata. Exiting...")
                stop(msg)
            }
            if (is.null(meta)){
                stop("Metadata is needed for paired testing!")
            }
            if (test != 'wilcoxon'){
              stop("Paired information can only be used with the Wilcoxon test")
            }
            # check that each entry in "paired" has two samples associated with
            # a different label, filter out the rest
            groups <- unique(meta[[paired]])
            if (verbose > 2) message("+++ Starting with ",
                                length(groups), " pairings")
            groups.red <- groups[vapply(groups, FUN = function(x){
                temp <- label$label[rownames(meta[meta[[paired]]==x,])]
                if (length(unique(temp))!=2){
                    return(FALSE)
                } else if (length(temp)!=2){
                    return(FALSE)
                }
                return(TRUE)
            }, FUN.VALUE = logical(1))]
            if (length(groups.red) > 5){
                if (verbose > 2) {
                    message("+++ Keeping ", length(groups.red),
                            " pairings with exactly two samples!")
                }
            } else {
                msg <- paste0("Per pairing, exactly 2 samples with different",
                    " label are needed! Only ", length(groups.red),
                    " pairing fulfill this requirement.")
                stop(msg)
            }
            meta.red <- meta[meta[[paired]] %in% groups.red,]
            feat <- feat[,rownames(meta.red)]
            label$label <- label$label[rownames(meta.red)]
            meta <- meta.red
        }
        # check the formula
        if (test=='wilcoxon' & formula!='feat~label'){
            warning("For confounder-corrected tests, SIAMCAT uses linear",
                    " mixed effect models. Changing test to 'lm'!")
            test <- 'lm'
        }

        param.list <- list(formula=formula, alpha=alpha, mult.corr=mult.corr,
                            log.n0=log.n0, pr.cutoff=pr.cutoff,
                            test=test, feature.type=feature.type,
                            paired=paired, probs.fc=probs.fc)
        # check if results are already available
        if (!is.null(associations(siamcat, verbose=0))){
            old.params <- assoc_param(siamcat)
            check <- any(
                all.equal(param.list[-which(names(param.list)=='alpha')],
                            old.params[-which(names(old.params)=='alpha')])
                == TRUE)
            check <- all(check, nrow(associations(siamcat)) == nrow(feat))
            check <- all(check,
                        all(rownames(associations(siamcat)) == rownames(feat)))

            if (check){
                message("+ Enrichments have already been calculated!")
                associations(siamcat) <- list(
                    assoc.results=associations(siamcat),
                    assoc.param=param.list)
                res <- associations(siamcat)
            } else {
                res <- analyze.markers(feat, meta, label, param.list)
                associations(siamcat) <- list(
                    assoc.results=res,
                    assoc.param=param.list)
            }
        } else {
            res <- analyze.markers(feat, meta, label, param.list)
            associations(siamcat) <- list(
                assoc.results=res,
                assoc.param=param.list)
        }

        if (verbose > 1){
            msg <- paste('+++ found', sum(res$p.adj < alpha, na.rm = TRUE),
                'significant associations at a significance level <', alpha)
            message(msg)
        }

        e.time <- proc.time()[3]
        if (verbose > 1){
            msg <- paste("+ finished check.associations in",
                formatC(e.time - s.time, digits = 3), "s")
            message(msg)
        }
        return(siamcat)
}


# ##############################################################################
#' @keywords internal
analyze.markers <- function(feat, meta, label, param.list){
    if (label$type=='CONTINUOUS'){
        res <- analyze.continuous.markes(feat, meta, label, param.list)
    } else if (label$type=='BINARY'){
        res <- analyze.binary.markers(feat, meta, label, param.list)
    }
    return(res)
}

# ##############################################################################
### maker analysis for two-class data
#     calculate p-value with Wilcoxon
#     fold change as normalized absolute difference between quantiles
#     prevalence shift
#     single marker AUC
#' @keywords internal
analyze.binary.markers <- function(feat, meta, label, param.list) {


    if (param.list$feature.type=='normalized'){
        take.log <- FALSE
    } else {
        take.log <- TRUE
    }

    if (any(feat[feat != 0] < param.list$log.n0) & take.log==TRUE){
        cnt <- length(which(feat[feat!=0] < param.list$log.n0))
        percentage <- (cnt/length(feat[feat!=0]))*100
        if (percentage >= 5){
            msg <- paste0('### Some values (',cnt, ' or ',
                formatC(percentage, digits=2),
                '% of non-zero entries',
                ') are smaller than the given detection limit!')
            message(msg)
        }
    }

    ############################################################################
    ### Calculate wilcoxon, pseudo-FC, prevalence shift, and AUC for all feats
    ############################################################################

    positive.label <- max(label$info)
    negative.label <- min(label$info)

    feat <- feat[,names(label$label)]

    pb <- progress_bar$new(total = nrow(feat))

    if (is.null(meta)){
        df.temp <- data.frame(label=label$label)
    } else {
        df.temp <- data.frame(lapply(names(meta),
                    FUN = function(x){meta@.Data[[which(names(meta)==x)]]}))
        colnames(df.temp) <- names(meta)
        rownames(df.temp) <- rownames(meta)
        df.temp$label <- label$label
    }
    effect.size <- t(vapply(rownames(feat), FUN = function(x){

        df.temp$feat <- feat[x,rownames(df.temp)]
        if (!is.null(param.list$paired)){
            df.temp <- df.temp[sort(df.temp[[param.list$paired]],
                                    index.return=TRUE)$ix,]
        }
        x.pos <- df.temp[df.temp$label==positive.label,'feat']
        x.neg <- df.temp[df.temp$label==negative.label,'feat']

        # prevalence
        pr.n <- mean(x.neg >= param.list$pr.cutoff)
        pr.p <- mean(x.pos >=  param.list$pr.cutoff)
        pr.shift <- c(pr.p - pr.n, pr.n, pr.p)

        # AUC
        temp <- roc(cases = x.pos, controls=x.neg, ci = TRUE, direction = '<')
        aucs <- c(temp$ci)

        # gFC
        if (take.log){
            x.pos <- log10(x.pos + param.list$log.n0)
            x.neg <- log10(x.neg + param.list$log.n0)
        }
        if (is.null(param.list$paired)){
            q.p <- quantile(x.pos, probs = param.list$probs.fc)
            q.n <- quantile(x.neg, probs = param.list$probs.fc)
            fc <- sum(q.p - q.n) / length(q.p)
        } else {
            fc <- mean(x.pos-x.neg)
        }

        # p.val
        if (param.list$test=='wilcoxon'){
            if (is.null(param.list$paired)) {
              p.val <- wilcox.test(formula=as.formula(param.list$formula),
                                   data=df.temp, exact=FALSE)$p.value
            } else {
              p.val <- wilcox.test(x.pos, x.neg, paired=TRUE, 
                                   exact=FALSE)$p.value
            }
        } else if (param.list$test=='lm'){
            if (take.log) df.temp$feat <- log10(df.temp$feat +
                                                    param.list$log.n0)
            if (!grepl("1\\|", param.list$formula)){
                fit <- lm(formula = as.formula(param.list$formula),
                            data=df.temp)
                res <- coefficients(summary(fit))
                p.val <- res['label',4]
            } else {
                fit <- suppressMessages(
                    lmerTest::lmer(formula=as.formula(param.list$formula),
                                    data=df.temp))
                res <- coefficients(summary(fit))
                p.val <- res['label',5]
            }
        }

        pb$tick()
        return(c('fc' = fc, 'p.val' = p.val, 'auc' = aucs[2],
                    'auc.ci.l' = aucs[1], 'auc.ci.h' = aucs[3],
                    'pr.shift' = pr.shift[1], 'pr.n' = pr.shift[2],
                    'pr.p' = pr.shift[3]))
    }, FUN.VALUE = double(8)))
    effect.size <- as.data.frame(effect.size)
    ### Apply multi-hypothesis testing correction

    if (param.list$mult.corr == 'none') {
        warning('No multiple hypothesis testing performed.')
        effect.size$p.adj <- effect.size$p.val
    } else {
        effect.size$p.adj <-
            p.adjust(effect.size$p.val, method = param.list$mult.corr)
    }

    return("effect.size" = effect.size)
}

# ##############################################################################
### maker analysis for regression
#     calculate p-value with LM
#     spearman correlation
#' @keywords internal
analyze.continuous.markes <- function(feat, meta, label, param.list) {

    if (param.list$feature.type=='normalized'){
        take.log <- FALSE
    } else {
        take.log <- TRUE
    }

    if (any(feat[feat != 0] < param.list$log.n0) & take.log==TRUE){
        cnt <- length(which(feat[feat!=0] < param.list$log.n0))
        percentage <- (cnt/length(feat[feat!=0]))*100
        if (percentage >= 5){
            msg <- paste0('### Some values (',cnt, ' or ',
                            formatC(percentage, digits=2),
                            '% of non-zero entries',
                            ') are smaller than the given detection limit!')
            message(msg)
        }
    }

    ############################################################################
    ### Calculate p value with LM, spearman's rho
    ############################################################################

    feat <- feat[,names(label$label)]

    pb <- progress_bar$new(total = nrow(feat))

    if (is.null(meta)){
        df.temp <- data.frame(label=label$label)
    } else {
        df.temp <- data.frame(lapply(names(meta),
            FUN = function(x){meta@.Data[[which(names(meta)==x)]]}))
        colnames(df.temp) <- names(meta)
        rownames(df.temp) <- rownames(meta)
        df.temp$label <- label$label[rownames(meta)]
    }


    effect.size <- t(vapply(rownames(feat), FUN = function(x){

        df.temp$feat <- feat[x,rownames(df.temp)]

        cor.sp <- cor(df.temp$feat, df.temp$label, method='spearman')

        # p.val
        if (take.log) df.temp$feat <- log10(df.temp$feat + param.list$log.n0)
        if (param.list$formula=='feat~label'){
            fit <- lm(formula = as.formula(param.list$formula),
                        data=df.temp)
            res <- coefficients(summary(fit))
            p.val <- res[2,4]
            fc <- res[2,1]
        } else {
            fit <- suppressMessages(
                lmer(formula=as.formula(param.list$formula),
                                data=df.temp))
            res <- coefficients(summary(fit))
            p.val <- res[2,5]
            fc <- res[2,1]
        }


        pb$tick()
        return(c('fc' = fc, 'p.val' = p.val, 'spearman' = cor.sp))
    }, FUN.VALUE = double(3)))
    effect.size <- as.data.frame(effect.size)
    ### Apply multi-hypothesis testing correction

    if (param.list$mult.corr == 'none') {
        warning('No multiple hypothesis testing performed.')
        effect.size$p.adj <- effect.size$p.val
    } else {
        effect.size$p.adj <-
            p.adjust(effect.size$p.val, method = param.list$mult.corr)
    }

    return("effect.size" = effect.size)
}
