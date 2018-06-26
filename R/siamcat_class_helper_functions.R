#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' Reset features in siamcat@phylose@otu_table to those in siamcat@orig_feat
#' @title reset.features
#' @name reset.features
#' @description Function reset features in siamcat@phylose@otu_table to those in
#' siamcat@orig_feat in an object of class \link{siamcat-class}
#' @param siamcat an object of class \link{siamcat-class}t
#' @return A new \link{siamcat-class} object
#' @export
#' @examples
#'     data(siamcat_example)
#'     siamcat_example <- reset.features(siamcat_example)
reset.features <- function(siamcat) {
    features(siamcat) <- orig_feat(siamcat)
    return(siamcat)
}

#' @title Filter samples from \code{siamcat@label}
#' @description This functions filters \code{siamcat@label}.
#' @param siamcat an object of class \link{siamcat-class}
#' @param ids names of samples to be left in the \code{siamcat@label}
#' @param verbose control output: \code{0} for no output at all, \code{1}
#' for more information about progress and success, defaults to \code{1}
#' @keywords filter.label
#' @export
#' @return siamcat an object of class \link{siamcat-class}
#' @examples
#'
#'     data(siamcat_example)
#'     # simple working example
#'     siamcat_filtered <- filter.label(siamcat_example, ids=c(1:10))
#'
filter.label <- function(siamcat, ids, verbose = 1) {
    label_old <- label(siamcat)
    labels_new <- list(
        label = label_old$label[ids],
        header = label_old$header,
        info = label_old$info,
        positive.lab = label_old$positive.lab,
        negative.lab = label_old$negative.lab,
        n.lab = label_old$n.lab,
        p.lab = label_old$p.lab
    )
    labels_new$n.idx <- labels_new$label == labels_new$negative.lab
    labels_new$p.idx <- labels_new$label == labels_new$positive.lab
    
    if (verbose > 0)
        message(paste(
            "Keeping labels of",
            length(labels_new$label),
            "sample(s)."
        ))
    
    label(siamcat) <- label(labels_new)
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
    if (!is.null(label(object)))
        cat(
            paste(
                "label()                label:           ",
                sum(label(object)$n.idx),
                label(object)$n.lab,
                "and",
                sum(label(object)$p.idx),
                label(object)$p.lab,
                "samples",
                sep = " "
            ),
            fill = TRUE
        )
    if (length(norm_param(object))) {
        cat(
            paste(
                "norm_param()           norm_param:       Features normalized
                using",
                norm_param(object)$norm.method,
                sep = " "
            ),
            fill = TRUE
        )
    }
    if (length(data_split(object)$num.folds)) {
        cat(
            paste(
                "data_split()            data_split:       ",
                data_split(object)$num.resample,
                "cv rounds with",
                data_split(object)$num.folds,
                "folds",
                sep = " "
            ),
            fill = TRUE
        )
    }
    if (length(model_type(object))) {
        cat(paste(
            "model_list()            model_list:       ",
            length(models(object)) , model_type(object),
                "models", sep = " "), fill = TRUE)
    }
    if (nrow(pred_matrix(object))) {
        cat(
            paste(
                "pred_matrix()           pred_matrix:       Predictions for",
                nrow(pred_matrix(object)),
                "samples from",
                ncol(pred_matrix(object)),
                "cv rounds",
                sep = " "
            ),
            fill = TRUE
        )
    }
    if (length(eval_data(object))) {
        cat(paste(
            "eval_data()             eval_data:         Average AUC:",
            round(eval_data(object)$auc.average[[1]], 3),
            sep = " "
        ),
            fill = TRUE)
    }
    
    # print otu_table (always there).
    cat("\ncontains phyloseq-class experiment-level object @phyloseq:",
        fill = TRUE)
    cat(
        paste(
            "phyloseq@otu_table()   OTU Table:         [ ",
            ntaxa(otu_table(physeq(object))),
            " taxa and ",
            nsamples(otu_table(physeq(object))),
            " samples ]",
            sep = ""
        ),
        fill = TRUE
    )
    
    # print Sample Data if there
    if (!is.null(sample_data(physeq(object), FALSE))) {
        cat(
            paste(
                "phyloseq@sam_data()    Sample Data:       [ ",
                dim(sample_data(physeq(object)))[1],
                " samples by ",
                dim(sample_data(physeq(object)))[2],
                " sample variables ]",
                sep = ""
            ),
            fill = TRUE
        )
    }
    
    # print tax Tab if there
    if (!is.null(tax_table(physeq(object), FALSE))) {
        cat(
            paste(
                "phyloseq@tax_table()   Taxonomy Table:    [ ",
                dim(tax_table(physeq(object)))[1],
                " taxa by ",
                dim(tax_table(physeq(object)))[2],
                " taxonomic ranks ]",
                sep = ""
            ),
            fill = TRUE
        )
    }
    
    # print tree if there
    if (!is.null(phy_tree(physeq(object), FALSE))) {
        cat(
            paste(
                "phyloseq@phy_tree()    Phylogenetic Tree: [ ",
                ntaxa(phy_tree(physeq(object))),
                " tips and ",
                phy_tree(physeq(object))$Nnode,
                " internal nodes ]",
                sep = ""
            ),
            fill = TRUE
        )
    }
    
    # print refseq summary if there
    if (!is.null(refseq(physeq(object), FALSE))) {
        cat(
            paste(
                "phyloseq@refseq()      ",
                class(refseq(physeq(object)))[1],
                ":
                [ ",
                ntaxa(refseq(physeq(object))),
                " reference sequences ]",
                sep = ""
            ),
            fill = TRUE
        )
    }
})
