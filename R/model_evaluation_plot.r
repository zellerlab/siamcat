#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Model Evaluation Plot
#' @description Produces two plots for model evaluation. The first plot shows
#'     the Receiver Operating Characteristic (ROC)-curves, the other the
#'     Precision-recall (PR)-curves for the different cross-validation
#'     repetitions.
#' @usage model.evaluation.plot(..., fn.plot = NULL,
#'     colours=NULL, verbose = 1)
#' @param ... one or more object of class \link{siamcat-class}, can be named
#' @param fn.plot string, filename for the pdf-plot
#' @param colours colour specification for the different \link{siamcat-class}-
#'     objects, defaults to \code{NULL} which will cause the colours to be
#'     picked from the \code{'Set1'} palette
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'     for only information about progress and success, \code{2} for normal
#'     level of information and \code{3} for full debug information,
#'     defaults to \code{1}
#' @keywords SIAMCAT model.evaluation.plot
#' @export
#' @return Does not return anything, but produces the model evaluation plot.
#' @examples
#' data(siamcat_example)
#'
#' # simple working example
#' model.evaluation.plot(siamcat_example, fn.plot='./eval.pdf')
#'
#' # plot several named SIAMCAT object
#' # (although we use only one example object here)
#' model.evaluation.plot('Example_1'=siamcat_example,
#'     'Example_2'=siamcat_example, colours=c('red', 'blue'),
#'     fn.plot='./eval.pdf')
model.evaluation.plot <- function(..., fn.plot=NULL, colours = NULL,
        verbose = 1) {
    if (verbose > 1)
        message("+ starting model.evaluation.plot")
    s.time <- proc.time()[3]

    if(!is.null(fn.plot)) pdf(fn.plot, onefile = TRUE)

    if (verbose > 2)
        message("+ plotting ROC")
    if (is.null(fn.plot)) par(ask=TRUE)
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    plot(
        NULL,
        xlim = c(0, 1),
        ylim = c(0, 1),
        xlab = "False positive rate",
        ylab = "True positive rate",
        type = "n"
    )
    title(paste("ROC curve for the model", sep = " "))
    abline(a = 0, b = 1, lty = 3)

    args <- list(...)

    if (length(args) > 1) {
        # checks
        if (!all(vapply(args, class,
            FUN.VALUE = character(1)) == 'siamcat')){
                stop("Please supply only SIAMCAT objects. Exiting...")
        }
        if (any(vapply(args,
            FUN=function(x){is.null(eval_data(x, verbose=0))},
            FUN.VALUE = logical(1)))){
                stop("Not all SIAMCAT objects have evaluation data. Exiting...")
        }

        n <- length(args)
        if (is.null(colours)) {
            if (n > 9) {
                colours <- colorRampPalette(brewer.pal(9, 'Set1'))(n)
                warning(paste0('Consider plotting fewer',
                                ' ROC-Curves into the same plot...'))
            } else if (n == 2) {
                colours <- brewer.pal(3, 'Set1')[rev(seq_len(2))]
            } else {
                colours <- brewer.pal(n, 'Set1')
            }
        }
        stopifnot(length(colours) == n)

        # ROC
        legend.val <- c()
        # plot each roc curve for each eval data object
        for (i in seq_along(args)) {
            legend.val <- c(legend.val,
                            as.numeric(single.roc.plot(args[[i]],
                                                        colours[i],
                                                        verbose=verbose)))
        }
        if (!is.null(names(args))) {
            legend('bottomright',
                    legend= paste0(names(args),
                                    ' AUC: ' ,
                                    format(legend.val, digits=3)),
                    col=colours, lty=1, lwd=2, cex=0.8, y.intersp=1.5)
        } else {
            legend('bottomright',
                    legend= paste0('AUC: ' ,
                                    format(legend.val, digits=3)),
                    col=colours, lty=1, lwd=2, cex=0.8, y.intersp=1.5)
        }

        # PR
        # precision recall curve
        if (verbose > 2)
            message("+ plotting PRC")
        plot(
            NULL,
            xlim = c(0, 1),
            ylim = c(0, 1),
            xlab = "Recall",
            ylab = "Precision",
            type = "n"
        )
        title(paste("Precision-recall curve for the model", sep = " "))
        legend.val <- c()
        # plot each roc curve for each eval data object
        for (i in seq_along(args)) {
            legend.val <- c(legend.val,
                            as.numeric(single.pr.plot(args[[i]],
                                                        colours[i],
                                                        verbose=verbose)))
        }
        if (!is.null(names(args))) {
            legend('bottomright',
                    legend= paste0(names(args),
                                    ' AUC: ' ,
                                    format(legend.val, digits=3)),
                    col=colours, lty=1, lwd=2, cex=0.8, y.intersp=1.5)
        } else {
            legend('bottomright',
                    legend= paste0('AUC: ' ,
                                    format(legend.val, digits=3)),
                    col=colours, lty=1, lwd=2, cex=0.8, y.intersp=1.5)
        }

    } else if (length(args) == 1) {
        # checks
        if (!all(class(args[[1]]) == 'siamcat'))
            stop('Please supply a SIAMCAT object. Exiting...')
        if(is.null(eval_data(args[[1]]))){
            stop('SIAMCAT object has no evaluation data. Exiting...')
        }

        # ROC
        if (is.null(colours)) colours <- 'black'
        auroc <- single.roc.plot(args[[1]], colours, verbose=verbose)
        if (data_split(args[[1]])$num.resample == 1) {
            text(0.7, 0.1, paste("AUC:", format(auroc, digits = 3)))
        } else {
            text(0.7, 0.1, paste("Mean-prediction AUC:",
                    format(auroc, digits = 3)))
        }

        # PR
        if (verbose > 2)
            message("+ plotting PRC")
        plot(
            NULL,
            xlim = c(0, 1),
            ylim = c(0, 1),
            xlab = "Recall",
            ylab = "Precision",
            type = "n"
        )
        title(paste("Precision-recall curve for the model", sep = " "))
        label <- label(args[[1]])
        abline(h = mean(label$label == max(label$info)),
            lty = 3)
        auprc <- single.pr.plot(args[[1]], colours, verbose=verbose)
        if (data_split(args[[1]])$num.resample == 1) {
            text(0.7, 0.1, paste("AUC:", format(auprc, digits = 3)))
        } else {
            text(0.7, 0.1, paste("Mean AUC:", format(auprc, digits = 3)))
        }
    } else {
        stop('No SIAMCAT object supplied. Exiting...')
    }

    if(!is.null(fn.plot)) tmp <- dev.off()
    if (is.null(fn.plot)) par(ask=FALSE)
    e.time <- proc.time()[3]
    if (verbose > 1)
        message(paste(
            "+ finished model.evaluation.plot in",
            formatC(e.time - s.time, digits = 3),
            "s"
        ))
    if (verbose == 1 & !is.null(fn.plot))
        message(paste(
            "Plotted evaluation of predictions successfully to:",
            fn.plot
        ))

}

single.pr.plot <- function(siamcat, colour, verbose) {

    eval.data <- eval_data(siamcat)

    # pr curves for resampling
    if (!is.null(eval.data$prc.all)) {
        aucspr.all = eval.data$auprc.all
        for (c in seq_len(length(eval.data$prc.all))) {
            pr = eval.data$prc.all[[c]]
            lines(pr$recall, pr$precision, col = alpha(colour, alpha=0.5))
            if (verbose > 2)
                message(paste(
                    "+++ AU-PRC (resampled run ",
                    c,
                    "): ",
                    format(aucspr.all[c], digits = 3)
                ))
        }
    }

    pr = eval.data$prc
    lines(pr$recall, pr$precision, col = colour, lwd = 2)
    auprc = eval.data$auprc


    if (!is.null(eval.data$roc.all)) {
        if (verbose > 1)
            message(
                paste(
                    "+ AU-PRC:\n+++ mean-prediction:",
                    format(auprc, digits = 3),
                    "\n+++ averaged       :",
                    format(mean(aucspr.all), digits = 3),
                    "\n+++ sd             :",
                    format(sd(aucspr.all), digits = 4)
                )
            )


    } else {
        if (verbose > 1)
            message("+ AU-PRC:", format(auprc, digits = 3), "\n")
    }
    return(auprc)
}

single.roc.plot <- function(siamcat, colour, verbose) {

    eval.data <- eval_data(siamcat)

    if (!is.null(eval.data$roc.all)){
        aucs = eval.data$auroc.all
        for (c in seq_along(eval.data$roc.all)) {
            roc.c = eval.data$roc.all[[c]]
            lines(1 - roc.c$specificities, roc.c$sensitivities,
                col = alpha(colour, alpha=0.5))
            if (verbose > 2) {
                message(paste('+++ AU-ROC (resampled run ',
                                c, "): ", format(aucs[c], digits=3)))
            }
        }
    }

    roc.summ = eval.data$roc
    lines(1 - roc.summ$specificities, roc.summ$sensitivities,
        col = colour, lwd = 2)
    auroc = eval.data$auroc

    # plot CI
    x = as.numeric(rownames(roc.summ$ci))
    yl = roc.summ$ci[, 1]
    yu = roc.summ$ci[, 3]
    polygon(1 - c(x, rev(x)), c(yl, rev(yu)),
        col = alpha(colour, alpha=0.1),
        border = NA)

    if (!is.null(eval.data$roc.all)){
        if (verbose > 1)
            message(
                paste(
                    "+ AU-ROC:\n+++ mean-prediction:",
                    format(auroc, digits = 3),
                    "\n+++ averaged       :",
                    format(mean(aucs), digits = 3),
                    "\n+++ sd             :",
                    format(sd(aucs), digits = 4)
                )
            )
    } else {
        if (verbose > 1)
            message(paste("+ AU-ROC:", format(auroc, digits = 3)))
    }

    return(as.numeric(auroc))
}
