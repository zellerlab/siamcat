#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#'@title Check and visualize associations between features and classes
#'
#'@description This function calculates for each feature a pseudo-fold change
#'     (geometrical mean of the difference between quantiles) between the
#'     different classes found in labels.
#'
#'     Significance of the differences is computed for each feature using a
#'     Wilcoxon test followed by multiple hypothesis testing correction.
#'
#'     Additionally, the Area Under the Receiver Operating Characteristic Curve
#'     (AU-ROC) and a prevalence shift are computed for the features found to be
#'     associated with the two different classes at a user-specified
#'     significance level \code{alpha}.
#'
#'     Finally, the function produces a plot of the top \code{max.show}
#'     associated features, showing the distribution of the log10-transformed
#'     abundances for both classes, and user-selected panels for the effect
#'     (AU-ROC, Prevalence Shift, and Fold Change)
#'
#'@usage check.associations(siamcat,fn.plot,color.scheme = "RdYlBu",
#'     alpha =0.05,mult.corr = "fdr", sort.by = "fc",detect.lim = 1e-06,
#'     pr.cutoff = 1e-6, max.show = 50, plot.type = "quantile.box",
#'     panels = c("fc","auroc"),verbose = 1)
#'
#'@param siamcat object of class \link{siamcat-class}
#'
#'@param fn.plot filename for the pdf-plot
#'
#'@param color.scheme valid R color scheme or vector of valid R colors (must be
#'     of the same length as the number of classes), defaults to \code{'RdYlBu'}
#'
#'@param alpha float, significance level, defaults to \code{0.05}
#'
#'@param mult.corr multiple hypothesis correction method, see
#'     \code{\link[stats]{p.adjust}}, defaults to \code{"fdr"}
#'
#'@param sort.by string, sort features by p-value (\code{"p.val"}), by fold
#'     change (\code{"fc"}) or by prevalence shift (\code{"pr.shift"}),
#'     defaults to \code{"fc"}
#'
#'@param detect.lim float, pseudocount to be added before log-transformation of
#'     the data, defaults to \code{1e-06}
#'
#'@param pr.cutoff float, cutoff for the prevalence computation, defaults to
#'     \code{1e-06}
#'
#'@param max.show integer, how many associated features should be shown,
#'     defaults to \code{50}
#'
#'@param plot.type string, specify how the abundance should be plotted, must be
#'     one of these: \code{c("bean", "box", "quantile.box", "quantile.rect")},
#'     defaults to \code{"quantile.box"}
#'
#'@param panels vector, name of the panels to be plotted next to the log10-
#'     transformed abundances, possible entries are \code{c("fc", "auroc",
#'     "prevalence")}, defaults to \code{c("fc", "auroc")}
#'
#'@param verbose control output: \code{0} for no output at all, \code{1} for
#'     only information about progress and success, \code{2} for normal level of
#'     information and \code{3} for full debug information, defaults to \code{1}
#'
#'@return Does not return anything, but produces an association plot
#'
#'@keywords SIAMCAT check.associations
#'
#'@export
#'
#'@examples
#'# Example data
#'     data(siamcat_example)
#'# since the whole pipeline has been run in the example data, exchange the
#'# normalized features with the original features
#'     siamcat_example <- reset.features(siamcat_example)
#'
#'# Simple example
#'     check.associations(siamcat_example, './assoc_plot.pdf')
#'
#'# Plot associations as bean plot
#'     check.associations(siamcat_example, './assoc_plot_bean.pdf',
#'     plot.type='bean')
#'
#'# Plot assocations as box plot
#'# Additionally, sort by p-value instead of by fold change
#'     check.associations(siamcat_example, './assoc_plot_fc.pdf',
#'     plot.type='box', sort.by='p.val')
#'
#'# Custom colors
#'     check.associations(siamcat_example, './assoc_plot_blue_yellow.pdf',
#'     plot.type='box', color.scheme=c('cornflowerblue', '#ffc125'))

check.associations <-
    function(siamcat,
        fn.plot,
        color.scheme = "RdYlBu",
        alpha = 0.05,
        mult.corr = "fdr",
        sort.by = "fc",
        detect.lim = 1e-06,
        pr.cutoff = 1e-6,
        max.show = 50,
        plot.type = "quantile.box",
        panels = c("fc", "auroc"),
        verbose = 1) {
# check panel and plot.type parameter
        if (verbose > 1)
            message("+ starting check.associations")
        s.time <- proc.time()[3]

        if (!all(panels %in% c("fc", "auroc", "prevalence"))) {
            stop("Unknown panel-type selected!")
        }
        if (length(panels) > 3) {
            warning(
                "Plot layout is not suited for more than 3 panels.
                Continuing with first three panels."
            )
            panels <- panels[seq_len(3)]
        }
        if ((!plot.type %in%
                c("bean", "box", "quantile.box", "quantile.rect")) ||
                length(plot.type) != 1) {
            warning("Plot type has not been specified properly! Continue with q
                uantile.box.")
            plot.type <- "quantile.box"
        }
# either give n_classes colors or color palette
        col <- check.color.scheme(color.scheme, label(siamcat))

        feat <- get.features.matrix(siamcat)
        label <- label(siamcat)

### Calculate different effect sizes
        if (verbose > 2)
            message("+++ analysing features\n")
        result.list <- analyse.binary.marker(
            feat = feat,
            label = label,
            detect.lim = detect.lim,
            colors = col,
            pr.cutoff = pr.cutoff,
            mult.corr = mult.corr,
            alpha = alpha,
            max.show = max.show,
            sort.by = sort.by,
            probs.fc = seq(.1, .9, .05),
            verbose = verbose
        )

###
        effect.size <- result.list$effect.size
        truncated <- result.list$truncated
        detect.lim <- result.list$detect.lim
        feat.red    <- result.list$feat.red
        feat.red.log <- log10(feat.red + detect.lim)

#############################################################################
### generate plots with significant associations between features and labels

# make plot matrix dependent on panels parameters
        if (verbose > 2)
            message("+++ preparing plotting layout")
        if (length(panels) == 3) {
            layout.mat <- cbind(2, 1, t(seq(3, length.out = length(panels))))
            widths <- c(0.5, 0.1, rep(0.4 / 3, length(panels)))
        } else {
            layout.mat <- cbind(2, 1, t(seq(3, length.out = length(panels))))
            widths <- c(0.5, 0.1, rep(0.2, length(panels)))
        }
        pdf(fn.plot,
            paper = 'special',
            height = 8.27,
            width = 11.69) # format:A4 landscape

        layout(mat = layout.mat, widths = widths)

#############################################################################
# PANEL 2: P-VALUES
# print p-values in second panel of the plot
        associations.pvals.plot(p.vals = effect.size$p.adj,
            alpha = alpha,
            verbose = verbose)

#############################################################################
# PANEL 1: DATA
# prepare margins
        associations.margins.plot(species_names = row.names(feat.red),
            verbose = verbose)


        if (verbose > 2)
            message("+++ plotting results")
        if (plot.type == "bean") {
            associations.bean.plot(feat.red.log,
                label,
                col = col,
                verbose = verbose)
        } else if (plot.type == "box") {
            associations.box.plot(feat.red.log,
                label,
                col = col,
                verbose = verbose)
        } else if (plot.type == "quantile.box") {
            associations.quantile.box.plot(feat.red.log,
                label,
                col = col,
                verbose = verbose)
        } else if (plot.type == "quantile.rect") {
            associations.quantile.rect.plot(feat.red.log,
                label,
                col = col,
                verbose = verbose)
        }

# plot title
        if (!truncated) {
            title(main = 'Differentially abundant features',
                xlab = 'Abundance (log10-scale)')
        } else {
            title(
                main = paste(
                    'Differentially abundant features\nshowing top',
                    max.show,
                    'features'
                ),
                xlab = 'Abundance (log10-scale)'
            )
        }

#############################################################################
# OTHER PANELS
        for (p in panels) {
            if (p == "fc") {
                associations.fcs.plot(
                    fc.all = effect.size$fc,
                    binary.cols = effect.size$bcol,
                    verbose = verbose
                )
            } else if (p == "prevalence") {
                associations.pr.shift.plot(
                    pr.shifts = effect.size[,c('pr.n', 'pr.p')],
                    col = col,
                    verbose = verbose
                )
            } else if (p == "auroc") {
                associations.aucs.plot(
                    aucs = effect.size[, c('auc', 'auc.ci.l', 'auc.ci.h')],
                    binary.cols = effect.size$bcol,
                    verbose = verbose
                )
            }
        }

# close pdf device
        tmp <- dev.off()
        e.time <- proc.time()[3]
        if (verbose > 1)
            message(paste(
                "+ finished check.associations in",
                formatC(e.time - s.time, digits = 3),
                "s"
            ))
        if (verbose == 1)
            message(paste(
                "Plotted associations between features and label
                successfully to:",
                fn.plot
            ))
    }


### one function for each type of plot
# bean plot
associations.bean.plot <-
    function(data.mat, label, col, verbose = 1) {
        if (verbose > 2)
            message("+ starting associations.bean.plot")
# create data.frame in format for beanplot
        bean.data <- data.frame(data = c(data.mat))
        bean.data$factor <- c(vapply(
            label$label,
            FUN = function(x) {
                paste(rownames(data.mat), names(label$info$class.descr[match(x,
                    label$info$class.descr)]))
            },
            FUN.VALUE = character(nrow(data.mat)),
            USE.NAMES = FALSE
        ))
# ensure correct ordering by converting to a factor
        bean.data$factor <- factor(bean.data$factor,
            levels = paste(rep(rownames(data.mat), each = 2),
                            c(label$n.lab, label$p.lab)))

        mn <- as.integer(c(min(bean.data$data)))
        mx <- as.integer(c(max(bean.data$data)))

        plot(
            NULL,
            xlab = '',
            ylab = '',
            xaxs = 'i',
            yaxs = 'i',
            axes = FALSE,
            xlim = c(mn - 1.5, mx + 1),
            ylim = c(0.45, nrow(data.mat) + 0.6),
            type = 'n'
        )
        ticks <- mn:mx
        for (v in ticks) {
            abline(v = v,
                lty = 3,
                col = 'lightgrey')
        }
        tick.labels <- formatC(10 ^ ticks, format = 'E', digits = 0)
        axis(
            side = 1,
            at = ticks,
            labels = tick.labels,
            cex.axis = 0.7
        )

        beanplot(
            data ~ factor,
            data = bean.data,
            side = "both",
            bw = "nrd0",
            col = list(col[1], col[2]),
            horizontal = TRUE,
            names = c(""),
            show.names = FALSE,
            beanlines = "median",
            maxstripline = 0.2,
            what = c(FALSE, TRUE, TRUE, FALSE),
            axes = FALSE,
            add = TRUE
        )


        legend(
            'topright',
            legend = c(label$p.lab, label$n.lab),
            fill = rev(col),
            bty = 'n'
        )
        associations.labels.plot(rownames(data.mat),
            plot.type = 'bean',
            verbose = verbose)
        if (verbose > 2)
            message("+ finished associations.bean.plot")
    }

# box plot
associations.box.plot <-
    function(data.mat, label, col, verbose = 1) {
        if (verbose > 2)
            message("+ starting associations.box.plot")
        box.colors <- rep(c(col[1], col[2]), nrow(data.mat))

# create data.frame in format for beanplot
        plot.data <- data.frame(data = c(data.mat))
        plot.data$factor <- c(vapply(
            label$label,
            FUN = function(x) {
                paste(rownames(data.mat), names(label$info$class.descr[match(x,
                    label$info$class.descr)]))
            },
            FUN.VALUE = character(nrow(data.mat)),
            USE.NAMES = FALSE
        ))
# ensure correct ordering by converting to a factor
        plot.data$factor <- factor(plot.data$factor,
            levels = paste(rep(rownames(data.mat), each = 2),
                            c(label$n.lab, label$p.lab)))

        mn <- as.integer(c(min(data.mat)))
        mx <- as.integer(c(max(data.mat)))

        plot(
            NULL,
            xlab = '',
            ylab = '',
            xaxs = 'i',
            yaxs = 'i',
            axes = FALSE,
            xlim = c(mn - 0.2, mx + 1),
            ylim = c(+0.5, nrow(data.mat) * 2 + 0.5),
            type = 'n'
        )
        ticks <- mn:mx
        for (v in ticks) {
            abline(v = v,
                lty = 3,
                col = 'lightgrey')
        }
        boxplot(
            data ~ factor,
            data = plot.data,
            horizontal = TRUE,
            names = c(""),
            show.names = FALSE,
            col = box.colors,
            axes = FALSE,
            outcol = c(col[1], col[2]),
            add = TRUE
        )


        tick.labels <- formatC(10 ^ ticks, format = 'E', digits = 0)
        axis(
            side = 1,
            at = ticks,
            labels = tick.labels,
            cex.axis = 0.7
        )
        legend(
            'topright',
            legend = c(label$p.lab, label$n.lab),
            fill = rev(col),
            bty = 'n'
        )
        associations.labels.plot(row.names(data.mat),
            plot.type = 'box',
            verbose = verbose)
        if (verbose > 2)
            message("+ finished associations.box.plot")
    }

# quantile.box plot
#' @keywords internal
associations.quantile.box.plot <- function(data.mat, label, col,
    verbose = 1) {
    if (verbose > 2)
        message("+ starting associations.quantile.box.plot")
    pos.col <- col[2]
    neg.col <- col[1]

    n.spec <- nrow(data.mat)
    p.m = min(data.mat, na.rm = TRUE)
    plot(
        rep(p.m, n.spec),
        seq_len(n.spec),
        xlab = '',
        ylab = '',
        yaxs = 'i',
        axes = FALSE,
        xlim = c(p.m, 0),
        ylim = c(0.5, n.spec + 0.5),
        frame.plot = FALSE,
        type = 'n'
    )
    for (v in seq(p.m, -1, 1)) {
        abline(v = v,
            lty = 3,
            col = 'lightgrey')
    }

    tck = floor(p.m):0
    axis(
        1,
        tck,
        formatC(10 ^ tck, format = 'E', digits = 0),
        las = 1,
        cex.axis = 0.7
    )

# get quantiles
    quant.probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
    quantiles.pos = rowQuantiles(data.mat[, label$p.idx], probs = quant.probs,
        na.rm = TRUE)
    quantiles.neg = rowQuantiles(data.mat[, label$n.idx], probs = quant.probs,
        na.rm = TRUE)

# inter-quartile range
    associations.quantiles.plot(quantiles.pos, up = TRUE, pos.col)
    associations.quantiles.plot(quantiles.neg, up = FALSE, neg.col)


# scatter plot on top
    for (i in seq_len(n.spec)) {
        pos.col.t <- change.transparency(pos.col)
        neg.col.t <- change.transparency(neg.col)

        points(
            data.mat[i, label$p.idx],
            rep(i + 0.15,
                sum(label$p.idx)) + rnorm(sum(label$p.idx), sd = 0.03),
            pch = 16,
            cex = 0.6,
            col = pos.col.t
        )
        points(
            data.mat[i, label$n.idx],
            rep(i - 0.15,
                sum(label$n.idx)) + rnorm(sum(label$n.idx), sd = 0.03),
            pch = 16,
            cex = 0.6,
            col = neg.col.t
        )
    }
    legend(
        'topright',
        legend = c(label$p.lab, label$n.lab),
        fill = rev(col),
        bty = 'n'
    )
    associations.labels.plot(row.names(data.mat),
        plot.type = 'quantile.box',
        verbose = verbose)
    if (verbose > 2)
        message("+ finished associations.quantile.box.plot")
}

# quantile.rect plot
#'@keywords internal
associations.quantile.rect.plot <-
    function(data.mat, label, col, verbose = 1) {
        if (verbose > 2)
            message("+ starting associations.quantile.rect.plot")
        n.spec <- nrow(data.mat)
        quant.probs <- seq(from = 0.1, to = 0.9, by = 0.1)

        quantiles.pos = rowQuantiles(data.mat[, label$p.idx],
                                        probs = quant.probs,
                                        na.rm = TRUE)
        quantiles.neg = rowQuantiles(data.mat[, label$n.idx],
                                        probs = quant.probs,
                                        na.rm = TRUE)

        p.mn <- min(data.mat, na.rm = TRUE)
        p.mx <- max(data.mat, na.rm = TRUE)

        plot(
            rep(p.mn, n.spec),
            seq_len(n.spec),
            xlab = '',
            ylab = '',
            yaxs = 'i',
            axes = FALSE,
            xlim = c(p.mn, p.mx),
            ylim = c(0, n.spec),
            frame.plot = FALSE,
            type = 'n'
        )
        for (v in seq(p.mn, 0, 1)) {
            abline(v = v,
                lty = 3,
                col = 'lightgrey')
        }

        tck = floor(p.mn):0
        axis(
            1,
            tck,
            formatC(10 ^ tck, format = 'E', digits = 0),
            las = 1,
            cex.axis = 0.7
        )

# create different tints of the colours
        colors.p <-
            rev(create.tints(vec = seq(0, 1, length.out = 4),
            colour = col[2]))
        colors.n <-
            rev(create.tints(vec = seq(0, 1, length.out = 4),
            colour = col[1]))

        associations.quantile.rect.sub.plot(quantiles.pos, up = TRUE, colors.p)
        associations.quantile.rect.sub.plot(quantiles.neg, up = FALSE, colors.n)
        associations.quantile.median.sub.plot(quantiles.pos, up = TRUE)
        associations.quantile.median.sub.plot(quantiles.neg, up = FALSE)

        mtext(
            'Quantiles',
            3,
            line = 0,
            at = 1,
            adj = 1.675,
            padj = 0.45,
            las = 1,
            cex = 0.7
        )
        legend(
            -1.75,
            n.spec,
            legend = c(
                "40%-60%",
                "30%-70%",
                "20%-80%",
                "10%-90%",
                "median",
                "",
                "",
                "",
                "",
                ""
            ),
            bty = 'n',
            cex = 1,
            fill = c(rev(colors.p), 'white', rev(colors.n), 'white'),
            lwd <- 1.3,
            ncol = 2,
            border = c(
                "black",
                "black",
                "black",
                "black",
                "white",
                "black",
                "black",
                "black",
                "black",
                "white"
            )
        )
        legend(
            -1.675,
            n.spec,
            legend = c("", "", "", "", ""),
            bty = 'n',
            lty = c(0, 0, 0, 0, 0),
    # cap legend size for diamond (should look symmetric to other symbols)
            pch = 18,
            cex = 1,
            pt.cex = c(0, 0, 0, 0, min(35 / n.spec, 2.25))
        )
        associations.labels.plot(rownames(data.mat),
            plot.type = 'quantile.rect',
            verbose = verbose)
        if (verbose > 2)
            message("+ finished associations.quantile.rect.plot")
}

### Prepare margins for the first plots make left margin as big as the
### longest label or maximally 20.1 lines
#'@keywords internal
associations.margins.plot <-
    function(species_names, p.label, verbose = 1) {
        if (verbose > 2)
            message("+ starting associations.margins.plot")
        cex.org <- par()$cex
        par(mar = c(5.1, 18, 4.1, 1.1), cex = 1)
        temp = par()$mai
        cex.labels <- min(.7, (((
            par()$pin[2] / length(species_names)
        ) * .6) /
                max(
                    strheight(species_names, units = 'inches')
                )))
        max_name <- max(strwidth(species_names, units = 'inches',
            cex = cex.labels)) + temp[4]
        temp[2] <- min(temp[2], max_name)
        par(mai = temp, cex = cex.org)
        if (verbose > 2)
            message("+ finished associations.margins.plot")
    }

### Plot single feature AUCs in single panel
#'@keywords internal
associations.aucs.plot <- function(aucs, binary.cols, verbose = 1) {
    if (verbose > 2)
        message("+ starting associations.aucs.plot")
# set margins
    par(mar = c(5.1, 0, 4.1, 1.6))
# plot background
    plot(
        NULL,
        xlab = '',
        ylab = '',
        xaxs = 'i',
        yaxs = 'i',
        axes = FALSE,
        xlim = c(0, 1),
        ylim = c(0.5, nrow(aucs) + 0.5),
        type = 'n'
    )
    ticks             <- seq(0, 1.0, length.out = 5)
    tick.labels <- formatC(ticks, digits = 2)
# plot gridlines
    for (v in ticks) {
        abline(v = v,
            lty = 3,
            col = 'lightgrey')
    }
# make thicker line at .5
    abline(v = .5, lty = 1, col = 'lightgrey')
# plot single feature aucs
    for (i in seq_len(nrow(aucs))) {
        segments(
            x0 = aucs[i, 2],
            x1 = aucs[i, 3],
            y0 = i,
            col = 'lightgrey',
            lwd = 1.5
        )
        points(aucs[i, 1], i, pch = 18, col = binary.cols[i])
        points(aucs[i, 1],
            i,
            pch = 5,
            col = 'black',
            cex = 0.9)
    }

# Title and axis label
    axis(
        side = 1,
        at = ticks,
        labels = tick.labels,
        cex.axis = 0.7
    )
    title(main = 'Feature AUCs', xlab = 'AU-ROC')
    if (verbose > 2)
        message("+ finished associations.aucs.plot")
}

### Plot fold changes in single panel
#'@keywords internal
associations.fcs.plot <-
    function(fc.all, binary.cols,    verbose = 1) {
        if (verbose > 2)
            message("+ starting associations.fcs.plot")
# margins
        par(mar = c(5.1, 0, 4.1, 1.6))
# get minimum and maximum fcs
        mx <- max(ceiling(abs(
            range(fc.all, na.rm = TRUE, finite = TRUE)
        )))
        mn        <- -mx
# plot background
        plot(
            NULL,
            xlab = '',
            ylab = '',
            xaxs = 'i',
            yaxs = 'i',
            axes = FALSE,
            xlim = c(mn, mx),
            ylim = c(0.2, length(fc.all) + 0.2),
            type = 'n'
        )
        grid(NULL, NA, lty = 3, col = 'lightgrey')
# plot bars
        barplot(
            fc.all,
            horiz = TRUE,
            width = 0.6,
            space = 2 / 3,
            col = binary.cols,
            axes = FALSE,
            add = TRUE,
            names.arg = FALSE
        )
# gridlines and axes labels
        ticks <- seq(from = mn,
            to = mx,
            length.out = 5)
        tick.labels <- formatC(ticks, digits = 2)
        axis(
            side = 1,
            at = ticks,
            labels = tick.labels,
            cex.axis = 0.7
        )
        title(main = 'Fold change', xlab = 'Generalized fold change')
        if (verbose > 2)
            message("+ finished associations.fcs.plot")
    }

### Plot prevalence shifts in single panel
associations.pr.shift.plot <-
    function(pr.shifts, col, verbose = 1) {
        if (verbose > 2)
            message("+ starting associations.pr.shift.plot")
# margins
        par(mar = c(5.1, 0, 4.1, 1.6))

# plot background
        plot(
            NULL,
            xlab = '',
            ylab = '',
            xaxs = 'i',
            yaxs = 'i',
            axes = FALSE,
            xlim = c(0, 1),
            ylim = c(0.2, nrow(pr.shifts) + 0.2),
            type = 'n'
        )
# gridlines and axes labels
        ticks <- seq(from = 0,
            to = 1,
            length.out = 5)
        for (v in ticks) {
            abline(v = v,
                lty = 3,
                col = 'lightgrey')
        }
        tick.labels <- formatC(ticks * 100, digits = 3)
        axis(
            side = 1,
            at = ticks,
            labels = tick.labels,
            cex.axis = 0.7
        )
# plot bars
        row.names(pr.shifts) <- NULL
        barplot(
            t(pr.shifts),
            horiz = TRUE,
            axes = FALSE,
            add = TRUE,
            space = c(0, 4 / 3),
            beside = TRUE,
            width = .3,
            col = c(col[1], col[2])
        )        
        title(main = 'Prevalence shift', xlab = 'Prevalence [%]')
        if (verbose > 2)
            message("+ finished associations.pr.shift.plot")
    }

# p-vals
#'@keywords internal
associations.pvals.plot <- function(p.vals, alpha,    verbose = 1) {
    if (verbose > 2)
        message("+ starting associations.pvals.plot")
# margins
    par(mar = c(5.1, .0, 4.1, 1.6))
    p.vals.log <- -log10(p.vals)
# get minimum and maximum
    mx        <-
        max(ceiling(abs(
            range(p.vals.log, na.rm = TRUE, finite = TRUE)
        )))
    mn        <- 0
    p.vals.log[is.infinite(p.vals.log)] <- mx
# plot background
    plot(
        NULL,
        xlab = '',
        ylab = '',
        xaxs = 'i',
        yaxs = 'i',
        axes = FALSE,
        xlim = c(mn, mx),
        ylim = c(0.2, length(p.vals) + 0.2),
        type = 'n'
    )
    grid(NULL, NA, lty = 3, col = 'lightgrey')
# plot bars
    barplot(
        p.vals.log,
        horiz = TRUE,
        width = 0.6,
        space = 2 / 3,
        col = 'lightgrey',
        axes = FALSE,
        add = TRUE,
        names.arg = FALSE
    )
# gridlines and axes labels
    ticks <- seq(from = mn, to = mx)
    abline(v = -log10(alpha),
        lty = 1,
        col = 'red')
    tick.labels <- formatC(ticks, digits = 2)
    axis(
        side = 1,
        at = ticks,
        labels = tick.labels,
        cex.axis = 0.7
    )
    title(main = 'Significance', xlab = '-log10(adj. p value)')
    if (verbose > 2)
        message("+ finished associations.pvals.plot")
}


# check if a string is a valid r color reprensentation
# from stackoverflow: Check if character string is a valid color representation
# https://stackoverflow.com/questions/13289009
#'@keywords internal
is.color <- function(x) {
    vapply(
        x,
        FUN = function(z) {
            tryCatch(
                is.matrix(col2rgb(z)),
                error = function(e)
                    FALSE
            )
        },
        FUN.VALUE = logical(1)
    )
}

### check the user-supplied color scheme for validity
### color scheme may either be a single RColorBrewer palette or a vector of
### the same length as the number of classes containing interpretable colors
### as strings
#'@keywords internal
check.color.scheme <- function(color.scheme, label, verbose = 1) {
    if (verbose > 2)
        message("+ starting check.color.scheme")
    n.classes = ifelse(label$info$type == 'BINARY', 2,
        length(unique(label$label)))

    if (length(color.scheme) == 1 &&
            class(color.scheme) == 'character') {
        if (n.classes == 2) {
    # if color scheme and binary label, make colors as before
            if (!color.scheme %in% row.names(brewer.pal.info)) {
                warning(
                    "Not a valid RColorBrewer palette name, defaulting to
                    RdBu.\n See brewer.pal.info for more information about
                    RColorBrewer palettes."
                )
                color.scheme <- 'RdYlBu'
            }
            colors <-
                rev(colorRampPalette(brewer.pal(
                    brewer.pal.info[color.scheme,
                        'maxcolors'], color.scheme
                ))(2))
        } else {
    # if color scheme and multiclass label, make colors either directly out
    # of the palette (if n.classes smaller than maxcolors) or like before
            if (!color.scheme %in% row.names(brewer.pal.info)) {
                warning(
                    "Not a valid RColorBrewer palette name, defaulting to
                    Set3.\n See brewer.pal.info for more information about
                    RColorBrewer palettes."
                )
                color.scheme <- 'Set3'
            }
    # if color scheme and multiclass label, check that the palette is not
    # divergent or sequential, but qualitative. Only issue warning.
            if (brewer.pal.info[color.scheme, 'category'] != 'qual')
                warning("Using a divergent or sequential color palette for
                        multiclass data.")
            if (n.classes <= brewer.pal.info[color.scheme, 'maxcolors']) {
                colors <- brewer.pal(n.classes, color.scheme)
            } else {
                warning("The data contains more classes than the color.palette
                    provides.")
                colors <-
                    rev(colorRampPalette(brewer.pal(
                        brewer.pal.info[color.scheme,
                            'maxcolors'], color.scheme
                    ))(n.classes))
            }
        }
        } else if (length(color.scheme == n.classes) &&
                all(is.color(color.scheme))) {
    # if colors, check that all strings are real colors and check that
    # the same length as n classes
    # convert color names to hex representation
            colors <-
                vapply(
                    color.scheme,
                    FUN = function(x) {
                        rgb(t(col2rgb(x)),
                            maxColorValue = 255)
                    },
                    FUN.VALUE = character(1),
                    USE.NAMES = FALSE
                )
        } else {
            stop("Supplied colors do not match the number of classes or are no
                valid colors")
        }
# add transparency
    colors <- vapply(
        colors,
        FUN = function(x) {
            paste0(x, '85')
        },
        FUN.VALUE = character(1),
        USE.NAMES = FALSE
    )
    if (verbose > 2)
        message("+ finished check.color.scheme")
    return(colors)
        }

#'@keywords internal
associations.labels.plot <-
    function(labels, plot.type,    verbose = 1) {
        if (verbose > 2)
            message("+ starting associations.labels.plot")
        adj <- rep(0, length(labels))
        if (plot.type == 'quantile.rect')
            adj <- rep(-0.5, length(labels))
        if (plot.type == 'box')
            adj <- -0.5 + seq_along(labels)
        cex.org <- par()$cex
        par(cex = 1)
        cex.labels <- min(.7, (((
            par()$pin[2] / length(labels)
        ) * .6) /
                max(strheight(labels, units = 'inches'))))
        for (i in seq_along(labels)) {
            mtext(
                labels[i],
                2,
                line = 0,
                at = i + adj[i],
                las = 1,
                cex = cex.labels
            )
        }
        par(cex = cex.org)
        if (verbose > 2)
            message("+ finished associations.labels.plot")
    }


### maker analysis for two-class data
#     calculate p-value with Wilcoxon
#     fold change as normalized absolute difference between quantiles
#     prevalence shift
#     single marker AUC
#'@keywords internal
analyse.binary.marker <- function(feat,
    label,
    detect.lim,
    colors,
    pr.cutoff,
    mult.corr,
    alpha,
    max.show,
    sort.by,
    probs.fc = seq(.1, .9, .05),
    verbose = 1) {
    if (verbose > 1)
        message("+ starting analyse.binary.marker")
    s.time <- proc.time()[3]
##############################################################################
### Calculate wilcoxon, pseudo-FC, prevalence shift, and AUC for each feature
##############################################################################
    if (verbose > 1)
        message('+++ calculating effect size for each feature.')
    if (is.null(detect.lim)) {
        warning(
            "Pseudo-count before log-transformation not supplied! Estimating it
            as 5% percentile.\n"
        )
        detect.lim <- quantile(feat[feat != 0], 0.05)
    }
    if (verbose)
        pb = txtProgressBar(max = nrow(feat), style = 3)
    effect.size <- data.frame(t(apply(
        feat,
        1,
        FUN = function(x) {
    # pseudo-fold change as differential quantile area
            q.p <-
                quantile(log10(x[label$p.idx] + detect.lim), probs = probs.fc)
            q.n <-
                quantile(log10(x[label$n.idx] + detect.lim), probs = probs.fc)
            fc <- sum(q.p - q.n) / length(q.p)

    # wilcoxon
            p.val <-
                wilcox.test(x[label$n.idx], x[label$p.idx],
                            exact = FALSE)$p.value

    # AU-ROC
            temp    <-
                roc(
                    predictor = x,
                    response = label$label,
                    ci = TRUE,
                    direction = '<'
                )
            aucs <- c(temp$ci)

    # prevalence shift
            temp.n <-
                sum(x[label$n.idx] >= pr.cutoff) / sum(label$n.idx)
            temp.p <-
                sum(x[label$p.idx] >= pr.cutoff) / sum(label$p.idx)
            pr.shift <- c(temp.p - temp.n, temp.n, temp.p)
            if (verbose)
                setTxtProgressBar(pb, (pb$getVal() + 1))
            return(
                c(
                    'fc' = fc,
                    'p.val' = p.val,
                    'auc' = aucs[2],
                    'auc.ci.l' = aucs[1],
                    'auc.ci.h' = aucs[3],
                    'pr.shift' = pr.shift[1],
                    'pr.n' = pr.shift[2],
                    'pr.p' = pr.shift[3]
                )
            )
        }
    )))

    effect.size$bcol <-
        ifelse(effect.size[, 'auc'] >= 0.5, colors[2], colors[1])

### Apply multi-hypothesis testing correction
    if (!tolower(mult.corr) %in%
        c('none', 'bonferroni', 'holm', 'fdr', 'bhy')) {
        stop(
            "! Unknown multiple testing correction method:', mult.corr,'
            Stopping!\n Must of one of c('none','bonferroni', 'holm','fdr',
            'bhy')"
        )
    }
    if (mult.corr == 'none') {
        warning('WARNING: No multiple hypothesis testing performed.')
        effect.size$p.adj <- effect.size$p.val
    } else {
        effect.size$p.adj <-
            p.adjust(effect.size$p.val, method = tolower(mult.corr))
    }

    if (verbose > 1)
        message(
            paste(
                '+++ found',
                sum(effect.size$p.adj < alpha,
                    na.rm = TRUE),
                'significant associations at a significance level <',
                alpha
            )
        )
    idx <- which(effect.size$p.adj < alpha)

    if (length(idx) == 0) {
        stop('No significant associations found. Stopping.\n')
    } else if (length(idx) < 5) {
        warning('Less than 5 associations found. Consider changing your alpha
            value')
    }

    idx <- idx[order(effect.size$p.adj[idx], decreasing = TRUE)]

# # truncated the list for the following plots
    truncated = FALSE
    if (length(idx) >= max.show) {
        truncated = TRUE
        idx <- idx[(length(idx) - max.show + 1):length(idx)]
        if (verbose > 1)
            message(
                paste(
                    '+++ truncating the list of significant
                    associations to the top',
                    max.show
                )
            )
    }

### Sort features
    if (verbose > 2)
        message('+++ sorting features')
    if (!sort.by %in% c('fc', 'p.val', 'pr.shift')) {
        if (verbose > 1)
            message(paste(
                '+++ Unknown sorting option:',
                sort.by,
                '. Instead order by fold change.'
            ))
        sort.by <- 'fc'
    }
    if (sort.by == 'fc') {
        fc.sign <- ifelse(effect.size[idx, 'fc'] == 0, 1,
            sign(effect.size[idx, 'fc']))
        p.adj.log <- -log10(effect.size$p.adj[idx]) * fc.sign
        idx <- idx[order(p.adj.log, decreasing = FALSE)]
    } else if (sort.by == 'p.val') {
        idx <- idx[order(effect.size$p.adj[idx], decreasing = TRUE)]
    } else if (sort.by == 'pr.shift') {
        pr.sign <- ifelse(effect.size[idx, 'pr.shift'] == 0, 1,
            sign(effect.size[idx, 'pr.shift']))
        p.adj.log <- -log10(effect.size$p.adj[idx]) * pr.sign
        idx <- idx[order(p.adj.log, decreasing = FALSE)]
    }
    e.time <- proc.time()[3]
    if (verbose > 1)
        message(paste(
            "+ finished analyse.binary.marker in",
            formatC(e.time - s.time, digits = 3),
            "s"
        ))
    return(
        list(
            "effect.size" = effect.size[idx, ],
            "feat.red" = feat[idx, , drop = FALSE],
            "truncated" = truncated,
            "detect.lim" = detect.lim
        )
    )
    }

#'@keywords internal
change.transparency <- function(col.name) {
    if (nchar(col.name) > 7) {
# adjust alpha channel by reducing transparency
        a = substr(col.name, nchar(col.name) - 1, nchar(col.name))
        a = 1 - (1 - as.numeric(paste('0x', a, sep = '')) / 255) / 2
        new.col = gsub('..$', toupper(as.hexmode(round(a * 255))), col.name)
    } else {
        new.col <- col.name
    }
    return(new.col)
}

#'@keywords internal
associations.quantiles.plot <- function(quantiles, up = TRUE, col) {
    n.spec <- nrow(quantiles)
    adj.y0 <- ifelse(up, 0, 0.3)
    adj.y1 <- ifelse(up, 0.3, 0)
# box
    rect(quantiles[, 2],
        seq_len(n.spec) - adj.y0,
        quantiles[, 4],
        seq_len(n.spec) + adj.y1,
        col = col)
# 90% interval
    segments(quantiles[, 1], seq_len(n.spec), quantiles[, 5], seq_len(n.spec))
    segments(
        quantiles[, 1],
        y0 = seq_len(n.spec) - adj.y0 / 3 * 2,
        y1 = seq_len(n.spec) + adj.y1 / 3 * 2
    )
    segments(
        quantiles[, 5],
        y0 = seq_len(n.spec) - adj.y0 / 3 * 2,
        y1 = seq_len(n.spec) + adj.y1 / 3 * 2
    )
# median
    segments(
        quantiles[, 3],
        y0 = seq_len(n.spec) - adj.y0,
        y1 = seq_len(n.spec) + adj.y1,
        lwd = 3
    )
}

#'@keywords internal
create.tints <- function(colour, vec) {
    new.cols <-
        vapply(
            vec,
            FUN = function(x) {
                rgb(matrix(col2rgb(colour) / 255 +
                        (1 - col2rgb(colour) / 255) * x, ncol = 3))
            },
            FUN.VALUE = character(1)
        )
    return(new.cols)
}

#'@keywords internal
associations.quantile.rect.sub.plot <-
    function(quantiles, up = TRUE, colors) {
        n.spec <- nrow(quantiles)
        adj.y0 <- ifelse(up, 0, 0.3)
        adj.y1 <- ifelse(up, 0.3, 0)
        for (i in seq_len(ncol(quantiles) / 2)) {
            rect(
                quantiles[, i],
                (0.5:n.spec) - adj.y0,
                quantiles[, ncol(quantiles) + 1 - i],
                (0.5:n.spec) + adj.y1,
                col = colors[i],
                border = c("black"),
                lwd = 0.9
            )
        }
}

#'@keywords internal
associations.quantile.median.sub.plot <-
    function(quantiles, up = TRUE) {
        n.spec <- nrow(quantiles)
        adj.y <- ifelse(up, 0.15,-0.15)
        points(
            quantiles[, ceiling(ncol(quantiles) / 2)],
            y = (0.5:n.spec) + adj.y,
            pch = 18,
            cex = min(35 / n.spec, 4)
        )
    }
