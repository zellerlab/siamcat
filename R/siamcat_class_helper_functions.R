#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Filter the label of a SIMACAT object
#'
#' @description This functions filters the label in a SIAMCAT object
#'
#' @param siamcat an object of class \link{siamcat-class}
#'
#' @param ids vector, can contain either names or indices of samples to be
#' retained
#'
#' @param verbose integer, control output: \code{0} for no output at all,
#'     \code{1} for only information about progress and success, \code{2} for
#'     normal level of information and \code{3} for full debug information,
#'     defaults to \code{1}
#'
#' @keywords filter.label
#'
#' @export filter.label
#'
#' @details This function filters the label contained in a SIAMCAT object,
#' based on the provided \code{ids}. The IDs can be either sample names or
#' indices to be retained.
#'
#' Predominantly for internal use...
#'
#' \strong{Please note:} It makes sense to run \link{validate.data} after
#' filtering the label.
#'
#' @return siamcat an object of class \link{siamcat-class}
#'
#' @examples
#' data(siamcat_example)
#'
#' # simple working example
#' siamcat_filtered <- filter.label(siamcat_example, ids=c(1:20))

filter.label <- function(siamcat, ids, verbose = 1) {
    label_old <- label(siamcat)
    # check ids
    if (all(is.numeric(ids))){
        if (any(ids <= 0) | any(ids > length(label_old$label))){
            stop('Provided IDs do not match the existing label! Exiting...\n')
        }
    } else if (all(is.character(ids))){
        if (!all(ids %in% names(label_old$label))){
            stop("Could not find all IDs in the label! Exiting...\n")
        }
    } else {
        stop("Please provide only one kind of ID",
            "(either sample names or indices)\n")
    }

    labels_new <- list(
        label = label_old$label[ids],
        info = label_old$info,
        type = label_old$type
    )

    if (verbose > 1)
        message(paste(
            "+ Keeping labels of",
            length(labels_new$label),
            "sample(s)."
        ))

    label(siamcat) <- labels_new
    return(siamcat)
}

# based on https://github.com/joey711/phyloseq/blob/master/R/show-methods.R
#' @title Show method for siamcat class object
#' @description Show method for siamcat class object
#' @rdname show-methods
#' @return none
#' @keywords internal
setMethod("show", "siamcat", function(object) {
    cat("siamcat-class object", fill = TRUE)


    # if it is a SIAMCAT.v1 model, print a warning
    if (is(models(object, verbose=0)[[1]], "WrappedModel")){
        message("Attention:\n",
                "This SIAMCAT object seems to have been constructed with ",
                "version 1.x, based on 'mlr'.\nYour current SIAMCAT version ",
                "has been upgraded to use 'mlr3' internally.\nPlease consider ",
                "re-training your SIAMCAT object or downgrading your SIAMCAT ",
                "version in order to continue.")
        stop("Exiting...")
    }
    # Label object
    if (!is.null(label(object)))
        label <- label(object)
        type <- label$type
        if (type == 'TEST'){
            n <- length(label$label)
            cat(paste("label()                Label object:        ",
                "Test label for", n, "samples"), fill=TRUE)
        } else if (type == 'BINARY'){
            p.lab <- names(which(label$info == max(label$info)))
            n.lab <- setdiff(names(label$info), p.lab)
            p.n <- length(which(label$label == max(label$info)))
            n.n <- length(which(label$label == min(label$info)))
            cat(paste("label()                Label object:        ",
                n.n, n.lab, "and", p.n, p.lab, "samples", sep = " "),
                fill = TRUE)
        } else if (type=='CONTINUOUS'){
            cat(paste("label()                Label object:        ",
                    length(label$label), "samples, ranging from",
                    sprintf(fmt='%.2f', label$info[1]), "to",
                    sprintf(fmt='%.2f', label$info[2]), sep = " "),
                fill = TRUE)
        }

    # filtered features
    if (!is.null(filt_feat(object, verbose=0))){
        temp <- filt_params(object)
        filtering.methods <- vapply(temp, FUN=function(x){x$filter.method},
            FUN.VALUE = character(1))
        if (length(filtering.methods) <= 3){
            cat(paste("filt_feat()            Filtered features:   ",
                nrow(get.filt_feat.matrix(object)), 'features after',
                paste(filtering.methods, collapse=', '), 'filtering'),
                    fill=TRUE)
        } else {
            cat(paste("filt_feat()            Filtered features:   ",
                nrow(get.filt_feat.matrix(object)), 'features after',
                length(filtering.methods), 'filtering steps'), fill=TRUE)
        }

    }

    # assocations testing
    if (!is.null(associations(object, verbose=0))){
        alpha <- assoc_param(object)$alpha
        temp.df <- associations(object)
        cat(paste("associations()         Associations:         Results from",
                "association testing\n                          ",
                "                 ",
                "with", length(which(temp.df$p.adj < alpha)),
                'significant features at alpha', alpha,
                sep=' '), fill=TRUE)
    }

    # normalizes features
    if (!is.null(norm_feat(object, verbose=0))) {
        cat(paste("norm_feat()            Normalized features: ",
        nrow(get.norm_feat.matrix(object)), "features normalized",
        "using", norm_params(object)$norm.method, sep = " "), fill = TRUE)
    }

    # data split
    if (!is.null(data_split(object, verbose=0))) {
        cat(paste("data_split()           Data split:          ",
            data_split(object)$num.resample, "cv rounds with",
            data_split(object)$num.folds, "folds", sep = " "), fill = TRUE)
    }

    # model list
    if (!is.null(model_type(object, verbose=0))) {
        cat(paste("model_list()           Model list:          ",
            length(models(object)) , model_type(object), "models", sep = " "),
            fill = TRUE)
    }
    # model list
    if (!is.null(feature_weights(object, verbose=0))) {
        cat(paste("feature_weights()      Feature weights:     ",
            "Summary of feature weights [ see also weight_matrix() ]",
            sep = " "),
            fill = TRUE)
    }

    # predictions
    if (!is.null(pred_matrix(object, verbose=0))) {
        cat(paste("pred_matrix()          Prediction matrix:   ",
            "Predictions for",
            nrow(pred_matrix(object)), "samples from",
            ncol(pred_matrix(object)), "cv rounds", sep = " "), fill = TRUE)
    }

    # evaluation data
    if (!is.null(eval_data(object, verbose=0))) {
        if ('auroc' %in% names(eval_data(object))){
            cat(paste( "eval_data()            Evaluation data:     ",
                        "Average AUC:", round(eval_data(object)$auroc, 3),
                        sep = " "), fill = TRUE)
        } else if ('r2' %in% names(eval_data(object))){
            cat(paste("eval_data()            Evaluation data:     ",
                    "Average R-squared:", round(eval_data(object)$r2, 3),
                    sep = " "), fill = TRUE)
        }
    }

    # print otu_table (always there).
    cat("\ncontains phyloseq-class experiment-level object @phyloseq:",
        fill = TRUE)
    cat(paste("phyloseq@otu_table()   OTU Table:            [ ",
        ntaxa(otu_table(physeq(object))), " taxa and ",
        nsamples(otu_table(physeq(object))), " samples ]", sep = ""),
        fill = TRUE)

    # print Sample Data if there
    if (!is.null(sample_data(physeq(object), FALSE))) {
        cat(paste("phyloseq@sam_data()    Sample Data:          [ ",
            nrow(sample_data(physeq(object))), " samples by ",
            ncol(sample_data(physeq(object))), " sample variables ]",
            sep = ""), fill = TRUE)
    }

    # print tax Tab if there
    if (!is.null(tax_table(physeq(object), FALSE))) {
        cat(paste("phyloseq@tax_table()   Taxonomy Table:       [ ",
            nrow(tax_table(physeq(object))), " taxa by ",
            ncol(tax_table(physeq(object))), " taxonomic ranks ]", sep = ""),
            fill = TRUE)
    }

    # print tree if there
    if (!is.null(phy_tree(physeq(object), FALSE))) {
        cat(paste("phyloseq@phy_tree()    Phylogenetic Tree:    [ ",
            ntaxa(phy_tree(physeq(object))), " tips and ",
            phy_tree(physeq(object))$Nnode, " internal nodes ]", sep = ""),
            fill = TRUE)
    }

    # print refseq summary if there
    if (!is.null(refseq(physeq(object), FALSE))) {
        cat(paste("phyloseq@refseq()      ", class(refseq(physeq(object)))[1],
            ": [ ", ntaxa(refseq(physeq(object))), " reference sequences ]",
            sep = ""), fill = TRUE)
    }
})
