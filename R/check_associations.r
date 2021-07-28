#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#'@title Check and visualize associations between features and classes
#'
#'@description This function computes different measures of association between
#'     features and the label and visualizes the results
#'
#'@usage check.associations(siamcat, fn.plot=NULL, color.scheme = "RdYlBu",
#'     alpha =0.05, mult.corr = "fdr", sort.by = "fc",
#'     detect.lim = 1e-06, pr.cutoff = 1e-6, max.show = 50,
#'     plot.type = "quantile.box",
#'     panels = c("fc","auroc"), prompt = TRUE,
#'     feature.type = 'filtered', paired=NULL, verbose = 1)
#'
#'@param siamcat object of class \link{siamcat-class}
#'
#'@param fn.plot string, filename for the pdf-plot. If \code{fn.plot} is
#'     \code{NULL}, the plot will be produced in the active graphics device.
#'
#'@param color.scheme valid R color scheme or vector of valid R colors (must be
#'     of the same length as the number of classes), defaults to \code{'RdYlBu'}
#'
#'@param alpha float, significance level, defaults to \code{0.05}
#'
#'@param mult.corr string, multiple hypothesis correction method, see
#'     \code{\link[stats]{p.adjust}}, defaults to \code{"fdr"}
#'
#'@param sort.by string, sort features by p-value (\code{"p.val"}), by fold
#'     change (\code{"fc"}) or by prevalence shift (\code{"pr.shift"}),
#'     defaults to \code{"fc"}
#'
#'@param detect.lim float, pseudocount to be added before log-transformation of
#'     the data, defaults to \code{1e-06}. Will be ignored if
#'     \code{feature.type} is \code{"normalized"}.
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
#'@param panels vector, name of the panels to be plotted next to the
#'     abundances, possible entries are \code{c("fc", "auroc",
#'     "prevalence")}, defaults to \code{c("fc", "auroc")}
#'
#'@param prompt boolean, turn on/off prompting user input when not plotting
#'      into a pdf-file, defaults to TRUE
#'
#'@param feature.type string, on which type of features should the function
#' work? Can be either \code{c()"original", "filtered", or "normalized")}.
#' Please only change this paramter if you know what you are doing!
#'
#' If \code{feature.type} is \code{"normalized"}, the normalized abundances
#' will not be log10-transformed.
#' 
#'@param paired character, column name of the meta-variable containing 
#' information for a paired test
#'      
#' @param verbose integer, control output: \code{0} for no output at all,
#'     \code{1} for only information about progress and success, \code{2} for
#'     normal level of information and \code{3} for full debug information,
#'     defaults to \code{1}
#'
#'@return object of class \link{siamcat-class} with the slot
#'     \code{associations} filled
#'
#'@keywords SIAMCAT check.associations
#'
#'@details For each feature, this function calculates different measures of
#'     association between the feature and the label. In detail, these
#'     associations are: \itemize{
#'     \item Significance as computed by a Wilcoxon test followed by multiple
#'     hypothesis testing correction.
#'     \item AUROC (Area Under the Receiver Operating Characteristics Curve)
#'     as a non-parameteric measure of enrichment (corresponds to the effect
#'     size of the Wilcoxon test).
#'     \item The generalized Fold Change (gFC) is a pseudo fold change
#'     which is calculated as geometric mean of the differences between the
#'     quantiles for the different classes found in the label.
#'     \item The prevalence shift between the two different classes found in
#'     the label.
#'     }
#'
#'     Finally, the function produces a plot of the top \code{max.show}
#'     associated features at a user-specified significance level \code{alpha},
#'     showing the distribution of the log10-transformed abundances for both
#'     classes, and user-selected panels for the effect (AU-ROC, Prevalence
#'     Shift, and Fold Change).
#'@export
#'
#'@examples
#' # Example data
#' data(siamcat_example)
#'
#' # Simple example
#' siamcat_example <- check.associations(siamcat_example,
#'     fn.plot='./assoc_plot.pdf')
#'
#' # Plot associations as box plot
#' siamcat_example <- check.associations(siamcat_example,
#'     fn.plot='./assoc_plot_box.pdf', plot.type='box')
#'
#' # Additionally, sort by p-value instead of by fold change
#' siamcat_example <- check.associations(siamcat_example,
#'     fn.plot='./assoc_plot_fc.pdf', plot.type='box', sort.by='p.val')
#'
#' # Custom colors
#' siamcat_example <- check.associations(siamcat_example,
#'     fn.plot='./assoc_plot_blue_yellow.pdf', plot.type='box',
#'     color.scheme=c('cornflowerblue', '#ffc125'))

check.associations <- function(siamcat, fn.plot=NULL, color.scheme = "RdYlBu",
    alpha = 0.05, mult.corr = "fdr", sort.by = "fc", detect.lim = 1e-06,
    pr.cutoff = 1e-6, max.show = 50, plot.type = "quantile.box",
    panels = c("fc", "auroc"), prompt=TRUE, feature.type='filtered',
    paired=NULL, verbose = 1) {

        if (verbose > 1)
            message("+ starting check.associations")
        s.time <- proc.time()[3]

        # check panel and plot.type parameter
        if (!all(panels %in% c("fc", "auroc", "prevalence"))) {
            stop("Unknown panel-type selected!")
        }
        panels <- unique(panels)
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
            warning("Plot type has not been specified properly! Continue with",
                    " quantile.box.")
            plot.type <- "quantile.box"
        }
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
        # check label
        label <- label(siamcat)
        meta <- meta(siamcat)
        if (label$type == 'TEST'){
            stop('Can not check assocations for a',
            ' SIAMCAT object with TEST label! Exiting...')
        }
        # check paired information 
        if (!is.null(paired)){
            if (!paired %in% colnames(meta)){
                stop(paste0("Column with pairing information not present in",
                            " the metadata. Exiting..."))
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
                stop(paste0("Per pairing, exactly 2 samples with different",
                            " label are needed! Only ", length(groups.red), 
                            " pairing fulfill this requirement."))
            }
            meta.red <- meta[meta[[paired]] %in% groups.red,]
            pairing.info <- meta.red[[paired]]
            names(pairing.info) <- rownames(meta.red)
            feat <- feat[,names(pairing.info)]
            label$label <- label$label[names(pairing.info)]
        } else {
            pairing.info <- NULL
        }
        
        # check fn.plot
        if (is.null(fn.plot)) {
            message(paste0('### WARNING: Not plotting to a pdf-file.\n',
                '### The plot is optimized for landscape DIN-A4 (or similar) ',
                'layout.\n### Please make sure that your plotting region is',
                ' large enough!!!\n### Use at your own risk...'))
            if (prompt == TRUE){
                continue <- askYesNo('Are you sure that you want to continue?',
                    default = TRUE,
                    prompts = getOption("askYesNo",
                        gettext(c("Yes", "No", "Cancel"))))
            } else {
                continue <- TRUE
            }
            if (!continue || is.na(continue)){
                opt <- options(show.error.messages = FALSE)
                on.exit(options(opt))
                stop('Exiting...')
            }
            par.old <- par(no.readonly=TRUE)
        }
        # either give n_classes colors or color palette
        col <- check.color.scheme(color.scheme, label)


        ### Calculate different effect sizes
        if (verbose > 2)
            message("+++ analysing features\n")
        probs.fc <- seq(.1, .9, .05)
        if (is.null(associations(siamcat, verbose=0))){
            result.list <- analyse.binary.marker(
                feat = feat,
                label = label,
                detect.lim = detect.lim,
                colors = col,
                pr.cutoff = pr.cutoff,
                mult.corr = mult.corr,
                alpha = alpha,
                probs.fc = probs.fc,
                take.log=ifelse(feature.type == 'normalized', FALSE, TRUE),
                pairing.info = pairing.info,
                verbose = verbose
            )
            # update siamcat
            associations(siamcat) <- list(
                assoc.results=result.list$effect.size,
                assoc.param=list(detect.lim=result.list$detect.lim,
                    pr.cutoff=pr.cutoff, probs.fc=probs.fc,
                    mult.corr=mult.corr, alpha=alpha,
                    feature.type=feature.type, paired=!is.null(paired)))
        } else {
            # if already existing, check parameters
            old.params <- assoc_param(siamcat)
            new.params <- list(detect.lim=detect.lim,
                pr.cutoff=pr.cutoff, probs.fc=probs.fc,
                mult.corr=mult.corr, alpha=alpha,
                feature.type=feature.type, paired=!is.null(paired))
            check <- any(all.equal(new.params, old.params) == TRUE)
            check <- all(check, nrow(associations(siamcat)) == nrow(feat))
            check <- all(check,
                all(rownames(associations(siamcat)) == rownames(feat)))
            # if the same, don't compute again but rather use the old resutls
            if (check){
                result.list <- list()
                result.list$effect.size <- associations(siamcat)
                result.list$detect.lim <- assoc_param(siamcat)$detect.lim
            } else {
                result.list <- analyse.binary.marker(
                    feat = feat,
                    label = label,
                    detect.lim = detect.lim,
                    colors = col,
                    pr.cutoff = pr.cutoff,
                    mult.corr = mult.corr,
                    alpha = alpha,
                    probs.fc = probs.fc,
                    take.log=ifelse(feature.type == 'normalized', FALSE, TRUE),
                    pairing.info = pairing.info,
                    verbose = verbose
                )
                # update siamcat
                associations(siamcat) <- list(
                    assoc.results=result.list$effect.size,
                    assoc.param=list(detect.lim=result.list$detect.lim,
                        pr.cutoff=pr.cutoff, probs.fc=probs.fc,
                        mult.corr=mult.corr, alpha=alpha,
                        feature.type=feature.type, paired=!is.null(paired)))
            }
        }

        ########################################################################
        # extract relevant info for plotting
        temp <- get.plotting.idx(result.list$effect.size, alpha=alpha,
            sort.by=sort.by, max.show=max.show, verbose=verbose)
        if (is.null(temp)){
            return(siamcat)
        }
        
        effect.size <- result.list$effect.size[temp$idx, , drop=FALSE]
        truncated <- temp$truncated
        detect.lim <- result.list$detect.lim
        feat.red    <- feat[temp$idx, , drop=FALSE]

        if (feature.type == 'normalized'){
            feat.plot <- feat.red
        } else {
            feat.red.log <- log10(feat.red + detect.lim)
            feat.plot <- feat.red.log
        }
        ########################################################################
        ### generate plots with significant associations between
        ##      features and labels

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
        if (!is.null(fn.plot)) {
            pdf(fn.plot,
                paper = 'special',
                height = 8.27,
                width = 11.69) # format:A4 landscape
            }

        layout(mat = layout.mat, widths = widths)

        ########################################################################
        # PANEL 2: P-VALUES
        # print p-values in second panel of the plot
        associations.pvals.plot(p.vals = effect.size$p.adj,
            alpha = alpha, mult.corr = mult.corr,
            verbose = verbose)

        ########################################################################
        # PANEL 1: DATA
        # prepare margins
        associations.margins.plot(species_names = row.names(feat.red),
            verbose = verbose)


        if (verbose > 2)
            message("+++ plotting results")
        if (plot.type == "bean") {
            associations.bean.plot(feat.plot,
                label,
                col = col,
                take.log=ifelse(feature.type == 'normalized', FALSE, TRUE),
                verbose = verbose)
        } else if (plot.type == "box") {
            associations.box.plot(feat.plot,
                label,
                take.log=ifelse(feature.type == 'normalized', FALSE, TRUE),
                col = col,
                verbose = verbose)
        } else if (plot.type == "quantile.box") {
            associations.quantile.box.plot(feat.plot,
                label,
                col = col,
                take.log=ifelse(feature.type == 'normalized', FALSE, TRUE),
                verbose = verbose)
        } else if (plot.type == "quantile.rect") {
            associations.quantile.rect.plot(feat.plot,
                label,
                col = col,
                take.log=ifelse(feature.type == 'normalized', FALSE, TRUE),
                verbose = verbose)
        }

        # plot title
        xlab <- ifelse(feature.type=='normalized',
                    'Normalized abundance', 'Abundance (log10-scale)')
        if (!truncated) {
            title(main = 'Differentially abundant features', xlab = xlab)
        } else {
            title(
                main = paste(
                    'Differentially abundant features\nshowing top',
                    max.show,
                    'features'
                ), xlab = xlab)
        }

        ########################################################################
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
        if (!is.null(fn.plot)) {
            tmp <- dev.off()
        } else {
            par(par.old)
        }
        e.time <- proc.time()[3]
        if (verbose > 1)
            message(paste(
                "+ finished check.associations in",
                formatC(e.time - s.time, digits = 3),
                "s"
            ))
        if (verbose == 1 & !is.null(fn.plot))
            message(paste(
                "Plotted associations between features and label",
                "successfully to:", fn.plot
            ))
        return(siamcat)
    }

# ##############################################################################
### AUC
#' @keywords internal
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
    ticks <- seq(0, 1.0, length.out = 5)
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

# ##############################################################################
### FC
#' @keywords internal
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
        mn <- -mx
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

# ##############################################################################
### PREVALENCE
#' @keywords internal
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

# ##############################################################################
# P-VALUES
#' @keywords internal
associations.pvals.plot <- function(p.vals, alpha,  mult.corr,  verbose = 1) {
    if (verbose > 2)
        message("+ starting associations.pvals.plot")
    # margins
    par(mar = c(5.1, .0, 4.1, 1.6))
    p.vals.log <- -log10(p.vals)
    # get minimum and maximum
    mx <-
        max(ceiling(abs(
            range(p.vals.log, na.rm = TRUE, finite = TRUE)
        )))
    mn <- 0
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
    if (mult.corr != 'none'){
        title(main = 'Significance', xlab = '-log10(adj. p value)')
    } else {
        title(main='Significance', xlab='-log10(p value)')
    }

    if (verbose > 2)
        message("+ finished associations.pvals.plot")
}

# ##############################################################################
# COLOR
# check if a string is a valid r color reprensentation
# from stackoverflow: Check if character string is a valid color representation
# https://stackoverflow.com/questions/13289009
#' @keywords internal
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
#' @keywords internal
check.color.scheme <- function(color.scheme, label, verbose = 1) {
    if (verbose > 2)
        message("+ starting check.color.scheme")
    n.classes = ifelse(label$type == 'BINARY', 2,
        length(unique(label$label)))

    if (length(color.scheme) == 1 &&
            is.character(color.scheme)) {
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

#' @keywords internal
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

#' @keywords internal
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

# ##############################################################################
# UTILITY FUNCTIONS
### Prepare margins for the first plots make left margin as big as the
### longest label or maximally 20.1 lines
#' @keywords internal
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

#' @keywords internal
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

#' @keywords internal
associations.quantiles.plot <- function(quantiles, up = TRUE, col) {
    n.spec <- nrow(quantiles)
    adj.y0 <- ifelse(up, 0, 0.3)
    adj.y1 <- ifelse(up, 0.3, 0)
    #  box
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

# ##############################################################################
# BEAN PLOT
#' @keywords internal
associations.bean.plot <-
    function(data.mat, label, col, take.log=TRUE, verbose = 1) {
        if (verbose > 2)
            message("+ starting associations.bean.plot")

        p.label <- max(label$info)
        n.label <- min(label$info)

        # create data.frame in format for beanplot
        bean.data <- data.frame(data = c(data.mat))
        bean.data$factor <- c(vapply(
            label$label,
            FUN = function(x) {
                paste(rownames(data.mat),
                names(label$info)[match(x, label$info)])
            },
            FUN.VALUE = character(nrow(data.mat)),
            USE.NAMES = FALSE
        ))
        # ensure correct ordering by converting to a factor
        bean.data$factor <- factor(bean.data$factor,
            levels = paste(rep(rownames(data.mat), each = 2),
                            names(label$info[order(label$info)])))

        mn <- floor(c(min(bean.data$data)))
        mx <- ceiling(c(max(bean.data$data)))

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
        if (take.log){
            tick.labels <- formatC(10 ^ ticks, format = 'E', digits = 0)
            axis(side = 1,at = ticks,labels = tick.labels,cex.axis = 0.7)
        } else {
            axis(1, ticks, cex.axis=0.7)
        }

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
            legend = c(names(which(label$info == p.label)),
                names(which(label$info == n.label))),
            fill = rev(col),
            bty = 'n'
        )
        associations.labels.plot(rownames(data.mat),
            plot.type = 'bean',
            verbose = verbose)
        if (verbose > 2)
            message("+ finished associations.bean.plot")
    }

# ##############################################################################
# BOX PLOT
#' @keywords internal
associations.box.plot <-
    function(data.mat, label, col, take.log=TRUE, verbose = 1) {
        if (verbose > 2)
            message("+ starting associations.box.plot")
        box.colors <- rep(c(col[1], col[2]), nrow(data.mat))

        p.label <- max(label$info)
        n.label <- min(label$info)

        # create data.frame in format for beanplot
        plot.data <- data.frame(data = c(data.mat))
        plot.data$factor <- c(vapply(
            label$label,
            FUN = function(x) {
                paste(rownames(data.mat),
                names(label$info)[match(x, label$info)])
            },
            FUN.VALUE = character(nrow(data.mat)),
            USE.NAMES = FALSE
        ))

        # ensure correct ordering by converting to a factor
        plot.data$factor <- factor(plot.data$factor,
            levels = paste(rep(rownames(data.mat), each = 2),
                            names(label$info[order(label$info)])))

        mn <- floor(c(min(data.mat)))
        mx <- ceiling(c(max(data.mat)))

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


        if (take.log){
            tick.labels <- formatC(10 ^ ticks, format = 'E', digits = 0)
            axis(side = 1,at = ticks,labels = tick.labels,cex.axis = 0.7)
        } else {
            axis(1, ticks, cex.axis=0.7)
        }
        legend(
            'topright',
            legend = c(names(which(label$info == p.label)),
                names(which(label$info == n.label))),
            fill = rev(col),
            bty = 'n'
        )
        associations.labels.plot(row.names(data.mat),
            plot.type = 'box',
            verbose = verbose)
        if (verbose > 2)
            message("+ finished associations.box.plot")
    }

# ##############################################################################
# QUANTILE BOX PLOT
#' @keywords internal
associations.quantile.box.plot <- function(data.mat, label, take.log=TRUE, col,
    verbose = 1) {
    if (verbose > 2)
        message("+ starting associations.quantile.box.plot")
    pos.col <- col[2]
    neg.col <- col[1]
    p.label <- max(label$info)
    n.label <- min(label$info)
    p.idx <- which(label$label == p.label)
    n.idx <- which(label$label == n.label)
    p.n <- length(which(label$label == p.label))
    n.n <- length(which(label$label == n.label))

    n.spec <- nrow(data.mat)
    if (take.log){
        p.min <- floor(min(data.mat, na.rm = TRUE))
        p.max <- 0
    } else {
        p.min <- floor(min(data.mat, na.rm = TRUE))
        p.max <- ceiling(max(data.mat, na.rm = TRUE))
    }

    plot(
        rep(p.min, n.spec),
        seq_len(n.spec),
        xlab = '',
        ylab = '',
        yaxs = 'i',
        axes = FALSE,
        xlim = c(p.min, p.max),
        ylim = c(0.5, n.spec + 0.5),
        frame.plot = FALSE,
        type = 'n'
    )
    for (v in seq(p.min, p.max, 1)) {
        abline(v = v,
            lty = 3,
            col = 'lightgrey')
    }

    tck = p.min:p.max
    if (take.log){
        axis(1, tck, formatC(10 ^ tck, format='E', digits=0),
            las=1, cex.axis=0.7)
    } else {
        axis(1, tck, las=1, cex.axis=0.7)
    }

    # get quantiles
    quant.probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
    quantiles.pos = rowQuantiles(data.mat[, p.idx, drop=FALSE],
        probs = quant.probs, na.rm = TRUE, drop=FALSE)
    quantiles.neg = rowQuantiles(data.mat[, n.idx, drop=FALSE],
        probs = quant.probs, na.rm = TRUE, drop=FALSE)

        # inter-quartile range
    associations.quantiles.plot(quantiles.pos, up = TRUE, pos.col)
    associations.quantiles.plot(quantiles.neg, up = FALSE, neg.col)


    # scatter plot on top
    for (i in seq_len(n.spec)) {
        pos.col.t <- change.transparency(pos.col)
        neg.col.t <- change.transparency(neg.col)

        points(
            data.mat[i, p.idx],
            rep(i + 0.15, p.n) + rnorm(p.n, sd = 0.03),
            pch = 16,
            cex = 0.6,
            col = pos.col.t
        )
        points(
            data.mat[i, n.idx],
            rep(i - 0.15, n.n) + rnorm(n.n, sd = 0.03),
            pch = 16,
            cex = 0.6,
            col = neg.col.t
        )
    }
    legend(
        'topright',
        legend = c(names(which(label$info == p.label)),
            names(which(label$info == n.label))),
        fill = rev(col),
        bty = 'n'
    )
    associations.labels.plot(row.names(data.mat),
        plot.type = 'quantile.box',
        verbose = verbose)
    if (verbose > 2)
        message("+ finished associations.quantile.box.plot")
}

# ##############################################################################
# QUANTILE RECT PLOT
#' @keywords internal
associations.quantile.rect.plot <-
    function(data.mat, label, col, take.log=TRUE, verbose = 1) {
        if (verbose > 2)
            message("+ starting associations.quantile.rect.plot")
        n.spec <- nrow(data.mat)
        quant.probs <- seq(from = 0.1, to = 0.9, by = 0.1)

        p.label <- max(label$info)
        n.label <- min(label$info)
        p.idx <- which(label$label == p.label)
        n.idx <- which(label$label == n.label)

        quantiles.pos = rowQuantiles(data.mat[, p.idx, drop=FALSE],
                                        probs = quant.probs,
                                        na.rm = TRUE, drop=FALSE)
        quantiles.neg = rowQuantiles(data.mat[, n.idx, drop=FALSE],
                                        probs = quant.probs,
                                        na.rm = TRUE, drop=FALSE)

        if (take.log){
            p.min <- floor(min(data.mat, na.rm = TRUE))
            p.max <- 0
        } else {
            p.min <- floor(min(data.mat, na.rm = TRUE))
            p.max <- ceiling(max(data.mat, na.rm = TRUE))
        }

        plot(
            rep(p.min, n.spec),
            seq_len(n.spec),
            xlab = '',
            ylab = '',
            yaxs = 'i',
            axes = FALSE,
            xlim = c(p.min, p.max),
            ylim = c(0, n.spec),
            frame.plot = FALSE,
            type = 'n'
        )
        for (v in seq(p.min, p.max, 1)) {
            abline(v = v,
                lty = 3,
                col = 'lightgrey')
        }

        tck = p.min:p.max
        if (take.log){
            axis(1, tck, formatC(10^tck, format='E', digits=0),
                las=1, cex.axis=0.7)
        } else {
            axis(1, tck, las=1, cex.axis=0.7)
        }

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

        legend(0.3*p.min, n.spec,
            legend = c("Quantiles", "40%-60%", "30%-70%", "20%-80%", "10%-90%",
                "median", "", "", "", "", ""),
            bty = 'n', cex = 1, fill = c('white', rev(colors.p), 'white',
                'white', rev(colors.n), 'white'),
            lwd = 1.3, ncol = 2, border = c("white", "black", "black",
                "black", "black", "white", "white", "black", "black", "black",
                "black", "white"))
        legend(0.3*p.min + abs(0.016*p.min), n.spec,
            legend = c("", "", "", "", "", ""), bty = 'n',
            lty = c(0, 0, 0, 0, 0, 0),
            # cap legend size for diamond (should look
            #   symmetric to other symbols)
            pch = 18, cex = 1,
            pt.cex = c(0, 0, 0, 0, 0, min(35 / n.spec, 2.25)))
        legend('bottomright',
            legend = c(names(which(label$info == max(label$info))),
                names(which(label$info == min(label$info)))),
            fill = rev(col), bty = 'n')
        associations.labels.plot(rownames(data.mat),
            plot.type = 'quantile.rect',
            verbose = verbose)
        if (verbose > 2)
            message("+ finished associations.quantile.rect.plot")
}

#' @keywords internal
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

#' @keywords internal
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

# ##############################################################################
### maker analysis for two-class data
#     calculate p-value with Wilcoxon
#     fold change as normalized absolute difference between quantiles
#     prevalence shift
#     single marker AUC
#' @keywords internal
analyse.binary.marker <- function(feat, label, detect.lim, colors,
    pr.cutoff, mult.corr, alpha, max.show, sort.by, probs.fc = seq(.1, .9, .05),
    take.log=TRUE, pairing.info=NULL, verbose = 1) {
    if (verbose > 1)
        message("+ starting analyse.binary.marker")
    s.time <- proc.time()[3]
    ############################################################################
    ### Calculate wilcoxon, pseudo-FC, prevalence shift, and AUC for all feats
    ############################################################################
    if (verbose > 1)
        message('+++ calculating effect size for each feature.')
    if (is.null(detect.lim) & take.log==TRUE) {
        warning(
            "Pseudo-count before log-transformation not supplied! Estimating it
            as 5% percentile.\n"
        )
        detect.lim <- quantile(feat[feat != 0], 0.05)
    }
    if (any(feat[feat != 0] < detect.lim) & take.log==TRUE){
        cnt <- length(which(feat[feat!=0] < detect.lim))
        percentage <- (cnt/length(feat[feat!=0]))*100
        if (percentage >= 5){
            warning(paste0('### Some values (',cnt, ' or ',
            formatC(percentage, digits=2),
            '% of non-zero entries',
            ') are smaller than the given detection limit!'))
        }
    }

    positive.label <- max(label$info)
    negative.label <- min(label$info)

    if (is.null(pairing.info)){
        feat <- feat[,names(label$label)]
    } else {
        feat <- feat[,names(pairing.info)]
    }
    
    if (verbose)
        pb <- progress_bar$new(total = nrow(feat))
    
    effect.size <- data.frame(t(apply(feat, 1, FUN = function(x) {
        # pseudo-fold change as differential quantile area
        if (is.null(pairing.info)){
            x.pos <- x[which(label$label==positive.label)]
            x.neg <- x[which(label$label==negative.label)]
            if (take.log == TRUE){
                q.p <- quantile(log10(x.pos + detect.lim), probs = probs.fc)
                q.n <- quantile(log10(x.neg + detect.lim), probs = probs.fc)
            } else {
                q.p <- quantile(x.pos, probs = probs.fc)
                q.n <- quantile(x.neg, probs = probs.fc)
            }
            fc <- sum(q.p - q.n) / length(q.p)

            # wilcoxon
            p.val <- wilcox.test(x.pos, x.neg, exact = FALSE)$p.value
            
            # AU-ROC
            temp <- roc(cases = x.pos, controls=x.neg, ci = TRUE,
                        direction = '<')
            aucs <- c(temp$ci)
            
            # prevalence shift
            temp.n <- mean(x.neg >= pr.cutoff)
            temp.p <- mean(x.pos >= pr.cutoff)
            pr.shift <- c(temp.p - temp.n, temp.n, temp.p)
        } else {
            pairing.info <- sort(pairing.info)
            x.pos <- x[names(which(label$label[names(pairing.info)] == 
                positive.label))]
            x.neg <- x[names(which(label$label[names(pairing.info)] == 
                negative.label))]
            if (take.log == TRUE){
                fc <- mean(log10(x.pos + detect.lim) - 
                    log10(x.neg + detect.lim))
            } else {
                fc <- mean(x.pos-x.neg)
            }
            p.val <- wilcox.test(x.pos, x.neg, exact=FALSE, 
                paired=TRUE)$p.value
            # AU-ROC
            temp <- roc(cases = x.pos, 
                        controls = x.neg,
                        ci = TRUE,
                        direction = '<')
            aucs <- c(temp$ci)
            
            # prevalence shift
            temp.n <- mean(x.neg >= pr.cutoff)
            temp.p <- mean(x.pos >= pr.cutoff)
            pr.shift <- c(temp.p - temp.n, temp.n, temp.p)
            
        }

        if (verbose)
            pb$tick()
        return(c('fc' = fc, 'p.val' = p.val, 'auc' = aucs[2],
            'auc.ci.l' = aucs[1], 'auc.ci.h' = aucs[3],
            'pr.shift' = pr.shift[1], 'pr.n' = pr.shift[2],
            'pr.p' = pr.shift[3]))
        }
    )))

    effect.size$bcol <-
        ifelse(effect.size[, 'auc'] >= 0.5, colors[2], colors[1])

    ### Apply multi-hypothesis testing correction
    if (!mult.corr %in% c("holm", "hochberg", "hommel", "bonferroni",
        "BH", "BY", "fdr", "none")) {
        stop("Unknown multiple testing correction method: '", mult.corr, 
            "'. Exiting!\n  Must of one of c('holm', 'hochberg', 'hommel', ", 
            "'bonferroni', 'BH', 'BY', 'fdr', none')")
    }
    if (mult.corr == 'none') {
        warning('WARNING: No multiple hypothesis testing performed.')
        effect.size$p.adj <- effect.size$p.val
    } else {
        effect.size$p.adj <-
            p.adjust(effect.size$p.val, method = mult.corr)
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

    e.time <- proc.time()[3]
    if (verbose > 1)
        message(paste(
            "+ finished analyse.binary.marker in",
            formatC(e.time - s.time, digits = 3),
            "s"
        ))
    return(
        list(
            "effect.size" = effect.size,
            "detect.lim" = detect.lim
        )
    )
}

#' @keywords internal
get.plotting.idx <- function(df.results, alpha, sort.by, max.show, verbose){

    idx <- which(df.results$p.adj < alpha)

    if (length(idx) == 0) {
        warning(paste0('No significant associations found.',
        ' No plot will be produced.\n'))
        return(NULL)
    } else if (length(idx) < 5) {
        warning(paste0('Less than 5 associations found. Consider',
        ' changing your alpha value.'))
    }

    idx <- idx[order(df.results$p.adj[idx], decreasing = TRUE)]

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
    if (!sort.by %in% c('fc', 'p.val', 'pr.shift', 'auc')) {
            message(paste0(
                '+++ Unknown sorting option: ',
                sort.by,
                '. Instead order by fold change.'
            ))
        sort.by <- 'fc'
    }
    if (sort.by == 'fc') {
        idx <- idx[order(df.results$fc[idx], decreasing = FALSE)]
    } else if (sort.by == 'p.val') {
        idx <- idx[order(df.results$p.adj[idx], decreasing = TRUE)]
    } else if (sort.by == 'pr.shift') {
        idx <- idx[order(df.results$pr.shift[idx], decreasing = FALSE)]
    } else if (sort.by == 'auc'){
        idx <- idx[order(df.results$auc[idx], decreasing = FALSE)]
    }
    return(list('idx'=idx,
                'truncated'=truncated))
}
