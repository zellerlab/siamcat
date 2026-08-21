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
#' select.feats=NULL, verbose = 1)
#'
#' @param siamcat object of class \link{siamcat-class}
#' 
#' @param formula string, formula used for testing, see Details for more 
#' information, defaults to \code{"feat~label"}
#' 
#' @param test string, statistical test used for the association testing, can
#' be either \code{'wilcoxon'}, \code{'lm'}, or \code{'lmer'} see Details for more 
#' information, defaults to \code{NULL}, which uses the Wilcoxon test for binary labels
#' and a linear model for continuous labels
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
#' @param select.feats character vector, names of a subset of features to
#' test for association, defaults to \code{NULL}, in which case all features
#' of the given \code{feature.type} are tested. If supplied, all entries must
#' match feature names present in the chosen feature matrix, otherwise the
#' function will stop with an error. This is useful to restrict testing to
#' a pre-specified set of features (e.g. features of interest) without
#' having to re-filter or re-normalize the whole SIAMCAT object.
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
#' implemented in \link[stats]{lm}) or a mixed-effect model (\link[lme4]{lmer}) can be used as well. For regression 
#' problems, the function defaults to the linear model.
#' 
#' @section Effect sizes:
#' The function calculates several measures for the effect size of the 
#' assocations between microbial features and the label. For binary 
#' classification problems, these associations are: \itemize{
#' \item AUROC (area under the Receiver Operating Characteristics curve) as a 
#' non-parametric measure of enrichment,
#' \item the generalized fold change (gFC) for the untransformed abundances, a pseudo-fold change which is 
#' calculated as geometric mean of the differences between quantiles across 
#' both groups (for paired tests the real fold change is used),
#' \item rank-biserial correlation, a non-parametric effect size metric 
#' \item prevalence shift (difference in prevalence between the two groups).}
#' \item effect size from the linear or mixed effect model (beta)
#' For regression problems, the effect sizes are: \itemize{
#' \item Spearman correlation between the feature and the label.}
#' \item effect size from the linear or mixed effect model (beta)
#' 
#' @section Confounder-corrected testing:
#' To correct for possible confounders while testing for association, the 
#' function uses linear or linear mixed effect models as implemented in the 
#' \link{stats} and \link{lme4} packages. To do so, the test formula needs to be adjusted 
#' to include the confounder. For example, when correcting for the metadata 
#' information \code{Sex} as a fixed effect, the formula would be: 
#' \code{'feat~label+Sex'} (see also the example below).
#' To treat sex as a random effect instead:
#' \code{'feat~label+(1|Sex)'} (see also the example below).
#' 
#' Please note that modifying the formula parameter in this function might
#' lead to unexpected results!
#'
#' @section Paired testing:
#' For paired testing, e.g. when the same patient has been sampled before and
#' after an intervention, the `paired` parameter can be supplied to the 
#' function. This indicated a column in the metadata table that holds the 
#' information about pairing. Note: this is applicable only for the Wilcoxon test.
#'
#' @section Selecting a feature subset:
#' By default, \code{check.associations} tests all features of the chosen
#' \code{feature.type}. If you only want to test association for a subset
#' of features (for example a pre-defined panel of taxa/genes), pass a
#' character vector of feature names via \code{select.feats}. All names
#' must be found among the rownames of the feature matrix corresponding to
#' \code{feature.type}, otherwise the function stops with an error.
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
#' #     formula='feat~label+Sex', test='lm')
#' 
#' # Paired testing
#' #
#' # this is not run during checks
#' # siamcat_paired <- check.associations(siamcat_paired, 
#' #     paired='Individual_ID')
#'
#' # Testing only a subset of features
#' #
#' # this is not run during checks
#' # siamcat_example <- check.associations(siamcat_example,
#' #     select.feats=c('Feature1', 'Feature2', 'Feature3'))
check.associations <- function(siamcat, formula="feat~label",
    test=NULL, alpha=0.05, mult.corr="fdr", log.n0=1e-06, pr.cutoff=1e-06,
    probs.fc=seq(.1, .9, .05), paired=NULL, feature.type='normalized',
    select.feats=NULL,
    verbose = 1) {
        canonical_formula_obj <- as.formula("feat ~ label")

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
        
        # check formulas
        formula_obj <- as.formula(formula)
        if (!'label' %in% attr(terms(formula_obj), "term.labels")) {
            stop("The formula must contain 'label' as a term.")
        }
        formula_obj <- as.formula(formula)
        random_effects_present <- !is.null(reformulas::findbars(formula_obj))

        if (!is.null(test)){
            if (length(test) != 1) {
                stop("Test must be a string of length 1 or NULL")
            }
        }
        
        # check label
        label <- label(siamcat)
        if (label$type == 'TEST'){
            stop('Can not check assocations for a',
                ' SIAMCAT object with TEST label! Exiting...')
        }
        if (label$type=='CONTINUOUS' && !is.null(test)){
            if (test == 'wilcoxon'){
                stop("Cannot test a SIAMCAT object with regression label using the Wilcoxon test.")
            }
        }
        
        # set NULL test
        if (is.null(test)) {
            if (label$type=='CONTINUOUS' && random_effects_present){
                test <- 'lmer'
            } else if (label$type=='CONTINUOUS' && !random_effects_present){
                test <- 'lm'
            } else if (label$type=='BINARY' && formula_obj == canonical_formula_obj){
                test <- 'wilcoxon'
            } else if (label$type=='BINARY' && (formula_obj != canonical_formula_obj & !random_effects_present)){
                test <- 'lm'
            } else if (label$type=='BINARY' && (formula_obj != canonical_formula_obj & random_effects_present)){
                test <- 'lmer'
            } else {
                stop("An error occurred in determining the type of test to perform, please raise an issue with the developpers of this package")
            }
            message("Setting test to ", test)
        }
        
        # check test
        allowed_tests <- c('wilcoxon', 'lm', 'lmer')
        if (!test %in% allowed_tests){
            stop("Unknown testing method: '", test,
                "'. Exiting!\n  Must of one of ", paste(allowed_tests, collapse=", "))
        } 
        if (test == "wilcoxon" && formula_obj != canonical_formula_obj) {
            stop("wilcoxon test does not support the use of covariates.")
        }
        if (test == "lm" && random_effects_present) {
            stop("lm test cannot be used with random effects in the formula")
        }
        if (test == "lmer" && !random_effects_present) {
            stop("lmer test cannot be used without random effects in the formula")
        }
        
        meta <- meta(siamcat)
        
        # check feature type
        if (!feature.type %in% c('original', 'filtered', 'normalized')){
            stop("Unrecognised feature type, exiting...\n")
        }
        
        # get features
        feat_orig <- get.orig_feat.matrix(siamcat)
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

        # check and apply select.feats
        if (!is.null(select.feats)){
            if (!is.character(select.feats)){
                stop("select.feats must be a character vector of feature names.")
            }
            if (length(select.feats) == 0){
                stop("select.feats must contain at least one feature name.")
            }
            if (any(duplicated(select.feats))){
                warning("Duplicated entries in select.feats detected, ",
                    "keeping only unique feature names.")
                select.feats <- unique(select.feats)
            }
            missing.feats <- setdiff(select.feats, rownames(feat))
            if (length(missing.feats) > 0){
                msg <- paste0("The following features in select.feats are ",
                    "not present among the '", feature.type,
                    "' features: ", paste(head(missing.feats, 10),
                    collapse=", "),
                    if (length(missing.feats) > 10) ", ..." else "",
                    ". Exiting...")
                stop(msg)
            }
            if (verbose > 1){
                message("+ restricting association testing to ",
                    length(select.feats), " out of ", nrow(feat),
                    " features via select.feats")
            }
            feat <- feat[select.feats, , drop=FALSE]
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

        param.list <- list(
            formula=formula, alpha=alpha, mult.corr=mult.corr,
            log.n0=log.n0, pr.cutoff=pr.cutoff,
            test=test, feature.type=feature.type,
            paired=paired, probs.fc=probs.fc,
            select.feats=select.feats
        )

        # if only alpha changed no need to rerun, just update the param.list
        if (!is.null(associations(siamcat, verbose=0))){
            old.params <- assoc_param(siamcat)
            check <- isTRUE(
                all.equal(
                    param.list[-which(names(param.list)=='alpha')],
                    old.params[-which(names(old.params)=='alpha')]
                )
            )
            check <- all(check, nrow(associations(siamcat)) == nrow(feat))
            if (check) {
                check <- all(check,
                        all(rownames(associations(siamcat)) == rownames(feat)))
            }

            if (check){
                message("+ Enrichments have already been calculated!")
                associations(siamcat) <- list(
                    assoc.results=associations(siamcat),
                    assoc.param=param.list)
                res <- associations(siamcat)
            } else {
                res <- analyze.markers(feat, feat_orig, meta, label, param.list)
                associations(siamcat) <- list(
                    assoc.results=res,
                    assoc.param=param.list)
            }
        } else {
            res <- analyze.markers(feat, feat_orig, meta, label, param.list)
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
analyze.markers <- function(feat, feat_orig, meta, label, param.list){
    # warn if the pseudocount is too large
    if (any(feat[feat != 0] < param.list$log.n0) & param.list$feature.type != 'normalized'){
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

    feat <- feat[,names(label$label)]
    
    if (is.null(meta)){
        df.temp <- data.frame(label=label$label)
    } else {
        if ("label" %in% names(meta)) stop("A column called 'label' in the metadata is not allowed.")
        df.temp <- data.frame(lapply(names(meta), FUN = function(x){meta@.Data[[which(names(meta)==x)]]}))
        colnames(df.temp) <- names(meta)
        rownames(df.temp) <- rownames(meta)
        df.temp$label <- label$label[rownames(df.temp)]
    }

    if (label$type=='CONTINUOUS'){
        ret <- analyze.continuous.markers(df.temp, feat, feat_orig, meta, label, param.list)
    } else if (label$type=='BINARY'){
        ret <- analyze.binary.markers(df.temp, feat, feat_orig, meta, label, param.list)
    }
    
    ret <- as.data.frame(ret)
    
    ### Apply multi-hypothesis testing correction
    if (param.list$mult.corr == 'none') {
        warning('No multiple hypothesis testing performed.')
        ret$p.adj <- ret$p.val
    } else {
        ret$p.adj <- p.adjust(ret$p.val, method = param.list$mult.corr)
    }

    return(ret)
}

###############################################################################
### maker analysis for two-class data
#' @keywords internal
analyze.binary.markers <- function(df.temp, feat, feat_orig, meta, label, param.list) {
    pb <- progress_bar$new(total = nrow(feat))
    positive.label <- max(label$info)
    negative.label <- min(label$info)
    formula_obj <- as.formula(param.list$formula)
    formula_null_obj <- update(formula_obj, . ~ . - label)

    ret <- t(vapply(rownames(feat), FUN = function(xname){
        df.temp$feat <- feat[xname,rownames(df.temp)]
        df.temp$feat_orig <- feat_orig[xname,rownames(df.temp)]
        df.temp$feat_log <- log10(df.temp$feat_orig + param.list$log.n0)

        if (param.list$feature.type != 'normalized') df.temp$feat <- log10(df.temp$feat + param.list$log.n0)
        
        # ensure correct pairing order
        if (!is.null(param.list$paired)){
            df.temp <- df.temp[order(df.temp[[param.list$paired]]),]
            pos_pairing <- df.temp[df.temp$label==positive.label, param.list$paired]
            neg_pairing <- df.temp[df.temp$label==negative.label, param.list$paired]
            if (
                !(
                    all(pos_pairing == neg_pairing) &&
                    length(unique(pos_pairing)) == length(pos_pairing) &&
                    length(unique(neg_pairing)) == length(neg_pairing)
                )
            ){
                stop("Pairing mismatch detected, please raise an issue with the maintainers of this package.")
            }
        }

        x.pos_orig <- df.temp[df.temp$label==positive.label,'feat_orig']
        x.neg_orig <- df.temp[df.temp$label==negative.label,'feat_orig']
        x.pos_log <- df.temp[df.temp$label==positive.label,'feat_log']
        x.neg_log <- df.temp[df.temp$label==negative.label,'feat_log']
        x.pos <- df.temp[df.temp$label==positive.label,'feat']
        x.neg <- df.temp[df.temp$label==negative.label,'feat']
        
        # prevalence
        pr.all <- mean(df.temp[, 'feat_orig'] > param.list$pr.cutoff)
        pr.p <- mean(x.pos_orig >= param.list$pr.cutoff)
        pr.n <- mean(x.neg_orig >= param.list$pr.cutoff)
        pr.shift <- pr.p - pr.n

        # AUC
        temp <- roc(cases = x.pos, controls = x.neg, ci = TRUE, direction = '<')
        aucs <- c(temp$ci) # strip attributes

        # FC
        if (is.null(param.list$paired)){
            q.p <- quantile(x.pos_log, probs = param.list$probs.fc)
            q.n <- quantile(x.neg_log, probs = param.list$probs.fc)
            fc_log10 <- mean(q.p - q.n)
        } else {
            fc_log10 <- mean(x.pos_log - x.neg_log)
        }
        
        # rank biserial correlation
        rank_biserial <- effectsize::rank_biserial(x.pos, x.neg, paired=!is.null(param.list$paired))$r_rank_biserial

        # p.val
        if (param.list$test=='wilcoxon'){
            beta <- NA
            p.val <- wilcox.test(x.pos, x.neg, paired=!is.null(param.list$paired), exact=FALSE)$p.value
        } else {
            if (param.list$test=='lm'){
                fit <- lm(formula=formula_obj, data=df.temp)
                fit_null <- lm(formula=formula_null_obj, data=df.temp)
                beta <- coef(fit)[['label']]
                p.val <- anova(fit_null, fit, test="LRT")[2, "Pr(>Chi)"]
            } else if (param.list$test == "lmer"){
                fit <- suppressMessages(
                    lme4::lmer(formula=formula_obj, data=df.temp, REML = FALSE)
                )
                fit_null <- suppressMessages(
                    lme4::lmer(formula=formula_null_obj, data=df.temp, REML = FALSE)
                )
                beta <- lme4::fixef(fit)[['label']]
                p.val <- anova(fit_null, fit, test="LRT")[2, "Pr(>Chisq)"]
            } else {
                stop("Unrecognised test, please raise an issue with the package developper.")
            }
        }

        pb$tick()
        return(c(
                'fc.log10' = fc_log10, 'p.val' = p.val,
                'beta' = beta, 'auc' = aucs[2], 
                'auc.ci.l' = aucs[1], 'auc.ci.h' = aucs[3],
                'pr.shift' = pr.shift, 'pr.n' = pr.n,
                'pr.p' = pr.p, 'pr.all' = pr.all,
                'rank.biserial' = rank_biserial
        ))
    }, FUN.VALUE = double(11)))

    # names are the xnames beacuse vapply sets those from the input
    return(ret)
}

# ##############################################################################
### maker analysis for regression
#' @keywords internal
analyze.continuous.markers <- function(df.temp, feat, feat_orig, meta, label, param.list) {
    pb <- progress_bar$new(total = nrow(feat))
    formula_obj <- as.formula(param.list$formula)
    formula_null_obj <- update(formula_obj, . ~ . - label)

    ret <- t(vapply(rownames(feat), FUN = function(xname){
        df.temp$feat <- feat[xname,rownames(df.temp)]
        df.temp$feat_orig <- feat_orig[xname,rownames(df.temp)]
        
        if (param.list$feature.type != 'normalized') df.temp$feat <- log10(df.temp$feat + param.list$log.n0)

        # prevalence and spearman
        pr.all <- mean(df.temp[, 'feat_orig'] > param.list$pr.cutoff)
        cor.sp <- cor(df.temp$feat, df.temp$label, method='spearman')
        
        if (param.list$test=='lm'){
            fit <- lm(formula=formula_obj, data=df.temp)
            fit_null <- lm(formula=formula_null_obj, data=df.temp)
            beta <- coef(fit)[['label']]
        } else if (param.list$test == "lmer"){
            fit <- suppressMessages(
                lme4::lmer(formula=formula_obj, data=df.temp, REML = FALSE)
            )
            fit_null <- suppressMessages(
                lme4::lmer(formula=formula_null_obj, data=df.temp, REML = FALSE)
            )
            beta <- lme4::fixef(fit)[['label']]
        } else {
            stop("Unrecognised test, please raise an issue with the package developper.")
        }

        p.val <- anova(fit_null, fit, test="LRT")[2, "Pr(>Chi)"]

        pb$tick()
        return(c(
            'p.val' = p.val, 'spearman' = cor.sp, 
            'beta' = beta, 'pr.all' = pr.all
        ))
    }, FUN.VALUE = double(4)))

    return(ret)
}