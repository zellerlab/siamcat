#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Visualize associations between features and classes as volcano plot
#'
#' @description This function creates a volcano plot to visualize the
#' association between features and the label
#'
#' @usage volcano.plot(siamcat)
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
#' @param color.scheme valid R color scheme or vector of valid R colors (must
#' be of length 3 for positive, negative, and non-significant associations),
#' defaults to \code{c('red', 'blue', 'gray')}
#'
#' @param annotate integer, number of features to annotate with the name
#'
#' @param annot.size integer, size of annotation text
#'
#' @param annot.y.exp float, percent vertical expansion of the plot area
#' to accommodate annotations
#' 
#' @param annot.force float, passed to geom_text_repel
#' 
#' @param annot.force_pull float, passed to geom_text_repel
#'
#' @param font.size integer, base font size for the plot
#'
#' @return Returns the ggplot plot object
#'
#' @keywords SIAMCAT volcano.plot
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
#' volcano.plot(siamcat_example, fn.plot = "./volcano.pdf")
volcano.plot <- function(
    siamcat, alpha = NULL, effect.metric = NULL, fn.plot = NULL,
    color.scheme = c("red", "blue", "gray"), annotate = 3,
    annot.size = 4, annot.y.exp = 0.2, annot.force=20, annot.force_pull=0, font.size = 14
) {
    associations <- associations(siamcat, verbose = 0)
    if (is.null(associations)) {
        stop(
            "SIAMCAT object does not contain association testing results! ",
            "Exiting..."
        )
    }
    associations$label <- rownames(associations)
    assoc.param <- assoc_param(siamcat)
    assoc.param$alpha <- ifelse(is.null(alpha), assoc.param$alpha, alpha)

    switch(
        as.character(length(color.scheme)),
        "1" = {
            tryCatch(
                {
                    col <- brewer.pal(3, color.scheme)[c(1, 3, 2)]
                },
                error = function(e) {
                    stop(
                        "color.scheme contains 1 element, so it is interpreted ",
                        "as a ColorBrewer palette, but the palette name is invalid."
                    )
                }
            )
        }, "3" = {
            tryCatch(col2rgb(color.scheme), error = function(e) {
                stop(
                    "color.scheme contains 3 elements, so it is interpreted ",
                    "as containing individual colors, but the color names are",
                    " invalid."
                )
            })
            col <- color.scheme
        },
        {
            stop(
                "color.scheme must contain 3 R color values or",
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

    ns.label <- "n. s."
    switch(label(siamcat)$type,
        "CONTINUOUS" = {
            fill_lab <- "Effect direction"
            associations$class <- ifelse(
                associations$p.adj < assoc.param$alpha,
                ifelse(associations$eff > eff.midpoint, "Positive", "Negative"),
                ns.label
            )
            names(col) <- c("Positive", "Negative", ns.label)
        },
        "BINARY" = {
            fill_lab <- "Enrichment"
            control.label <- names(label(siamcat)$info)[label(siamcat)$info == -1]
            case.label <- names(label(siamcat)$info)[label(siamcat)$info == 1]
            names(col) <- c(case.label, control.label, ns.label)
            associations$class <- ifelse(
                associations$p.adj < assoc.param$alpha,
                ifelse(associations$eff > eff.midpoint, case.label, control.label),
                ns.label
            )
        },
        {
            stop("label type is invalid. Please raise an issue with the package developers.")
        }
    )

    associations.to.label <- associations[associations$p.adj < assoc.param$alpha,]
    associations.to.label <- associations.to.label[order(associations.to.label$p.adj),]
    associations.to.label <- rbind(
        head(associations.to.label[associations.to.label$eff > eff.midpoint, ], annotate),
        head(associations.to.label[associations.to.label$eff < eff.midpoint, ], annotate)
    )

    mult.corr.label <- c(
        "holm" = "q",
        "hochberg" = "q",
        "hommel" = "q",
        "bonferroni" = "q",
        "BH" = "q",
        "BY" = "q",
        "fdr" = "q",
        "none" = "p"
    )[[assoc.param$mult.corr]]

    plot <- ggplot(
        associations,
        aes(
            x = eff, y = -log10(p.adj),
            size = pr.all, fill = class, label = label
        )
    ) +
        geom_point(shape = 21, alpha = 0.3) +
        ggrepel::geom_text_repel(
            data = associations.to.label,
            ylim = c(
                max(-log10(associations$p.adj)),
                max(-log10(associations$p.adj)) * (1 + annot.y.exp)
            ),
            show.legend = FALSE, size = annot.size, force = annot.force, force_pull = annot.force_pull
        ) +
        geom_hline(
            yintercept = -log10(assoc.param$alpha),
            color = "gray", lty = "dashed", lwd = 0.5
        ) +
        scale_fill_manual(
            values = col, breaks = names(col),
            guide = guide_legend(override.aes = list(size = 6))
        ) +
        scale_size(guide = guide_legend(reverse = TRUE)) +
        labs(
            x = xlab, y = bquote(-log[10](italic(.(mult.corr.label)))),
            size = "Prevalence", fill = fill_lab
        ) +
        theme_siamcat(font.size)

    # compute maximum y (hard to do a priori because of ggrepel)
    # this is for shifting the annotation of relative units to the plot area
    y_span <- diff(layer_scales(plot)$y$range$range)
    plot <- plot + annotate(
        "text",
        x = Inf, y = -log10(assoc.param$alpha) + 0.02*y_span,
        label = deparse(bquote(alpha ~ "=" ~ .(assoc.param$alpha))),
        color = "gray40", parse = TRUE, hjust = 1, vjust = 0
    )

    if (nrow(associations.to.label) != 0) {
        plot <- plot + scale_y_continuous(
            expand = expansion(mult = c(0.05, annot.y.exp))
        )
    }

    # save the plot
    if (!is.null(fn.plot)) {
        ggsave(
            fn.plot, plot, bg = "white", height = 6, width = 7
        )
    }

    return(plot)
}
