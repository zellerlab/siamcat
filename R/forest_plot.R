#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL

#' @title Visualize associations between features and classes as forest plot
#'
#' @description This function creates a forest plot to visualize the
#' association between features and the label
#'
#' @usage forest.plot(siamcat)
#'
#' @param siamcat object of class \link{siamcat-class}
#'
#' @param alpha float, override significance threshold from siamcat object
#' 
#' @param effect.metric string, effect metric to use for the x-axis.
#' If \code{NULL}, the default will be "beta" for linear models and "rank.biserial" for Wilcoxon tests.
#' Valid options for classification tasks are "fc.log10", "beta" (only for linear models), "rank.biserial", and "auc".
#' Valid options for regression tasks are "beta", "spearman".
#'
#' @param fn.plot string, filename for the plot (any extension supported by
#' ggsave is allowed). If \code{fn.plot} is \code{NULL}, the plot will only
#' be returned as a ggplot object.
#' 
#' @ param max.feat integer, maximum number of features to plot. Defaults to 30.
#' 
#' @ param max.pos.feat integer, maximum number of positively associated features to plot. Defaults to \code{NULL}.
#' If set, this parameter overrides max.feat.
#' 
#' @ param max.neg.feat integer, maximum number of negatively associated features to plot. Defaults to \code{NULL}.
#' If set, this parameter overrides max.feat.
#'
#' @param color.scheme valid R color scheme or vector of valid R colors (must
#' be of length 2 for positive and negative associations),
#' defaults to \code{NULL}, which uses the standard muted red and blue palette.
#'
#' @param font.size integer, base font size for the plot
#' 
#' @param feat.renamer function to use to rename the features for plotting.
#' Defaults to \code{NULL} where no renaming is applied.
#' 
#' @param title string or boolean or NULL, title for the plot.
#' If \code{NULL}, the label name will be used.
#' If \code{FALSE}, no title will be added.
#' 
#' @param width numeric, width of the plot in inches, defaults to 7
#' 
#' @param height numeric, height of the plot in inches, defaults to 6
#'
#' @return Returns the ggplot plot object
#'
#' @keywords SIAMCAT forest.plot
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
#' forest.plot(siamcat_example, fn.plot = "./forest.pdf")
forest.plot <- function(
    siamcat, alpha = NULL, effect.metric = NULL, fn.plot = NULL,
    color.scheme = NULL, max.feat = 30, max.pos.feat = NULL, max.neg.feat = NULL,
    feat.renamer = NULL, font.size = 14, title=NULL, width = 7, height = 6
) {
    if (is.null(color.scheme)) color.scheme <- pos_neg_neut_palette[1:2]
    associations <- associations(siamcat, verbose = 0)
    if (is.null(associations)) {
        stop(
            "SIAMCAT object does not contain association testing results! ",
            "Exiting..."
        )
    }
    associations$label <- rownames(associations)
    if (!is.null(feat.renamer)) {
        if (!is.function(feat.renamer)) stop("feat.renamer must be a function")
        newlabels <- feat.renamer(associations$label)
        if (!is.character(newlabels) && length(newlabels) == nrow(associations)) {
            stop ("feat.renamer must return a character vector of the same length as its input")
        }
        if (any(duplicated(newlabels))) {
            message("Duplicated values in feature names:")
            message(newlabels[duplicated(newlabels)])
            stop("feat.renamer must not return duplicate values")
        }
        associations$label <- newlabels
    }
    assoc.param <- assoc_param(siamcat)
    assoc.param$alpha <- ifelse(is.null(alpha), assoc.param$alpha, alpha)

    switch(
        as.character(length(color.scheme)),
        "1" = {
            tryCatch(
                {
                    col <- brewer.pal(2, color.scheme)
                },
                error = function(e) {
                    stop(
                        "color.scheme contains 1 element, so it is interpreted ",
                        "as a ColorBrewer palette, but the palette name is invalid."
                    )
                }
            )
        }, "2" = {
            tryCatch(col2rgb(color.scheme), error = function(e) {
                stop(
                    "color.scheme contains 2 elements, so it is interpreted ",
                    "as containing individual colors, but the color names are",
                    " invalid."
                )
            })
            col <- color.scheme
        },
        {
            stop(
                "color.scheme must contain 2 R color values or",
                " the name of 1 ColorBrewer palette."
            )
        }
    )

    if (is.null(effect.metric) && assoc.param$test %in% c("lm", "lmer")) {
        effect.metric <- "beta"
    } else if (is.null(effect.metric) && assoc.param$test == "wilcoxon") {
        effect.metric <- "rank.biserial"
    } else if (
        !(
            (label(siamcat)$type == "BINARY" && assoc.param$test == "wilcoxon" && effect.metric %in% c("fc.log10", "auc", "rank.biserial")) ||
            (label(siamcat)$type == "BINARY" && assoc.param$test %in% c("lm", "lmer") && effect.metric %in% c("fc.log10", "beta", "auc", "rank.biserial")) ||
            (label(siamcat)$type == "CONTINUOUS" && effect.metric %in% c("beta", "spearman"))
        )
    ) {
        stop(sprintf("effect.metric %s is invalid for label type %s and test %s", effect.metric, label(siamcat)$type, assoc.param$test))
    }

    switch(effect.metric,
        "fc.log10" = {
            xlab <- bquote(log[10] ~ "fold change")
            associations$eff <- associations$fc
            eff.midpoint <- 0
        },
        "beta" = {
            xlab <- bquote("Effect size (" ~ beta ~ ")")
            associations$eff <- associations$beta
            eff.midpoint <- 0
        },
        "rank.biserial" = {
            xlab <- "Rank biserial correlation"
            associations$eff <- associations$rank.biserial
            eff.midpoint <- 0
        },
        "spearman" = {
            xlab <- bquote(Spearman~rho)
            associations$eff <- associations$spearman
            eff.midpoint <- 0
        },
        "auc" = {
            xlab <- "Area under the ROC curve"
            associations$eff <- associations$auc
            eff.midpoint <- 0.5
        },
        {stop(sprintf("effect.metric %s is invalid", effect.metric))}
    )
    
    associations <- associations[associations$p.adj < assoc.param$alpha,]
    if (nrow(associations) == 0) stop("No significant associations to show.")
    associations <- associations[order(associations$p.adj),]

    switch(label(siamcat)$type,
        "CONTINUOUS" = {
            fill_lab <- "Effect direction"
            associations$class <- ifelse(associations$eff > eff.midpoint, "Positive", "Negative")
            names(col) <- c("Positive", "Negative")
        },
        "BINARY" = {
            fill_lab <- "Enrichment"
            control.label <- names(label(siamcat)$info)[label(siamcat)$info == -1]
            case.label <- names(label(siamcat)$info)[label(siamcat)$info == 1]
            names(col) <- c(case.label, control.label)
            associations$class <- ifelse(associations$eff > eff.midpoint, case.label, control.label)
        },
        {
            stop("label type is invalid. Please raise an issue with the package developers.")
        }
    )

    if (is.null(max.feat) && (is.null(max.pos.feat) || is.null(max.neg.feat))) {
        stop("max.feat cannot be NULL if either max.pos.feat and max.neg.feat are also NULL.")
        if (length(max.feat) != 1) stop("max.feat must be a single integer value.")
    }
    if (is.null(max.pos.feat) && is.null(max.neg.feat)) {
        associations <- head(associations, max.feat)
    } else {
        if (!is.null(max.pos.feat) && length(max.pos.feat) != 1) stop("max.pos.feat must be a single integer value.")
        if (!is.null(max.neg.feat) && length(max.neg.feat) != 1) stop("max.neg.feat must be a single integer value.")
        if (is.null(max.pos.feat)) {
            max.pos.feat <- max(max.feat - max.neg.feat, 0)
        } else if (is.null(max.neg.feat)) {
            max.neg.feat <- max(max.feat - max.pos.feat, 0)
        }
        associations <- rbind(
            head(associations[associations$eff > eff.midpoint, ], max.pos.feat),
            head(associations[associations$eff < eff.midpoint, ], max.neg.feat)
        )
    }

    associations$label <- factor(associations$label, levels = associations$label[order(associations$eff)])

    plot <- ggplot(
        associations,
        aes(x = eff, y = label, fill = class)
    ) +
        geom_col() +
        scale_fill_manual(
            values = col, breaks = names(col),
            guide = guide_legend(override.aes = list(size = 6))
        ) +
        labs(x = xlab, fill = fill_lab) +
        theme_siamcat_forest(font.size)


    if (!isFALSE(title)) {
        if (!is.character(title)) {
            title <- label(siamcat)$name
        }
        plot <- plot + ggtitle(title)
    }

    # save the plot
    if (!is.null(fn.plot)) {
        ggsave(
            fn.plot, plot, bg = "white", height = height, width = width
        )
    }

    return(plot)
}
