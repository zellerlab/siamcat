#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Visualize associations between features and classes as volcano plot
#'
#' @description This function creates a volcano plot to vizualize the 
#' association between features and the label
#'
#' @usage volcano.plot(siamcat, fn.plot=NULL, color.scheme="RdYlBu", 
#' annotate=5)
#'
#' @param siamcat object of class \link{siamcat-class}
#'
#' @param fn.plot string, filename for the pdf-plot. If \code{fn.plot} is
#' \code{NULL}, the plot will be produced in the active graphics device.
#'
#' @param color.scheme valid R color scheme or vector of valid R colors (must 
#' be of the same length as the number of classes), defaults to \code{'RdYlBu'}
#' 
#' @param annotate integer, number of features to annotate with the name
#'
#' @return Does not return anything, but produces a volcano plot based on
#' association measures
#'
#' @keywords SIAMCAT volcano.plot
#'
#' @details bla bla bla
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
#' volcano.plot(siamcat_example, fn.plot='./volcano.pdf')
volcano.plot <- function(siamcat, fn.plot=NULL, 
    color.scheme="RdYlBu", annotate=5){

    associations <- associations(siamcat, verbose = 0)
    if (is.null(associations)){
        stop("SIAMCAT object does not contain association testing results! ",
            "Exiting...")
    }
    assoc.param <- assoc_param(siamcat)
    if (assoc.param$feature.type == 'original'){
        feat <- get.orig_feat.matrix(siamcat)
        feat <- log10(feat + assoc.param$log.n0)
    } else if (assoc.param$feature.type == 'filtered'){
        feat <- get.filt_feat.matrix(siamcat)
        feat <- log10(feat + assoc.param$log.n0)
    } else if (assoc.param$feature.type == 'normalized'){
        feat <- get.norm_feat.matrix(siamcat)
    }

    mean.ab <- rowMeans(feat)
    associations$mean.ab <- mean.ab[rownames(associations)]
    if (assoc.param$feature.type != 'normalized'){
        associations$circle.size <- (associations$mean.ab -
            log10(assoc.param$log.n0)) *
            8/-log10(assoc.param$log.n0) + 1
    } else {
        associations$circle.size <- (associations$mean.ab +
            abs(min(associations$mean.ab))) *
            6/max(associations$mean.ab) + 1
    }
    # check color.scheme
    col <- check.color.scheme.volcano(color.scheme)

    # if fn.plot
    if (!is.null(fn.plot)){
        pdf(fn.plot, paper = 'special', height = 6, width = 6)
    }

    # plot
    plot(associations$fc, -log10(associations$p.adj),
        xlab = "Effect size",
        ylab='-log10(adjusted P value)', type = 'n')
    abline(h=-log10(assoc.param$alpha), lty=3)
    abline(v=0, lty=3)

    # alpha value colouring
    non.signif <- associations[associations$p.adj > assoc.param$alpha,]
    symbols(x=non.signif$fc, y=-log10(non.signif$p.adj),
        circles=non.signif$c, inches=1/9,
        bg=ifelse(non.signif$fc > 0, alpha(col[2], 0.5), alpha(col[1], 0.5)),
        fg=alpha('black', 0.4), add=TRUE)
    if (any(associations$p.adj <= assoc.param$alpha)){
        signif <- associations[associations$p.adj <= assoc.param$alpha,]
        symbols(x=signif$fc, y=-log10(signif$p.adj),
            circles=signif$c, inches=1/9,
            bg=ifelse(signif$fc > 0, alpha(col[2], 0.85), alpha(col[1], 0.85)),
            fg=alpha('black', 0.8), add=TRUE)


        # annotate top x features
        if (annotate > 0){
            signif <- signif[sort(signif$p.adj, index.return=TRUE)$ix,]
            if (nrow(signif) > annotate){
                signif.red <- signif[seq_len(annotate),]
            } else {
                message("Fewer significant features at alpha ", 
                    assoc.param$alpha,
                    " than desired features for annotation (", annotate, ")!")
                signif.red <- signif
            }
            for (i in seq_len(nrow(signif.red))){
                text(signif.red$fc[i], -log10(signif.red$p.adj[i]),
                    rownames(signif.red)[i],
                col=ifelse(signif.red$fc[i] > 0, alpha(col[2], 0.5), 
                    alpha(col[1], 0.5)),
                pos=ifelse(signif.red$fc[i] > 0, 2, 4))
            }
        }
    }

    if (!is.null(fn.plot)){
        tmp <- dev.off()
    }
}

check.color.scheme.volcano <- function(color.scheme) {

    if (length(color.scheme) == 1 &&
        is.character(color.scheme)) {
    
        # if color scheme and binary label, make colors as before
        if (!color.scheme %in% row.names(brewer.pal.info)) {
            warning("Not a valid RColorBrewer palette name! Defaulting to",
                "RdYlBu.\n See brewer.pal.info for more information about",
                " RColorBrewer palettes.")
            color.scheme <- 'RdYlBu'
        }
        colors <- rev(colorRampPalette(
            brewer.pal(brewer.pal.info[color.scheme, 'maxcolors'],
                color.scheme))(2))
    } else if (length(color.scheme < 2) && all(is.color(color.scheme))) {
        # if colors, check that all strings are real colors and check that
        # the same length as n classes
        # convert color names to hex representation
        colors <- vapply(color.scheme, FUN = function(x) {
            rgb(t(col2rgb(x)), maxColorValue = 255)}, FUN.VALUE = character(1),
            USE.NAMES = FALSE)
    } else {
        stop("Not enough colors or no valid colors supplied")
    }
    return(colors)
}
