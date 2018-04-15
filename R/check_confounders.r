#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Check for potential confounders in the metadata
#' @description This function checks for associations between class labels and
#'        potential confounders (e.g. age, sex, or BMI) that are present in the
#'        metadata. Statistical testing is performed with Fisher's exact test
#'        or Wilcoxon test, while associations are visualized either as barplot
#'        or Q-Q plot, depending on the type of metadata.
#' @param siamcat an object of class \link{siamcat-class}
#' @param fn.plot string, filename for the pdf-plot
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'        for only information about progress and success, \code{2} for normal
#'        level of information and \code{3} for full debug information,
#'        defaults to \code{1}
#' @keywords SIAMCAT check.confounders
#' @export
#' @return Does not return anything, but produces a single plot for each
#'        metadata category.
#' @examples
#'  # Example data
#'  data(siamcat_example)
#'  # since the whole pipeline has been run in the example data, exchange the
#'  # normalized features with the original features
#'  siamcat_example <- reset.features(siamcat_example)
#'
#'  # Simple working example
#'  check.confounders(siamcat_example, './conf_plot.pdf')
#'
#'  # Additional information with verbose
#'  \dontrun{check.confounders(siamcat_example, './conf_plot.pdf', verbose=2)}
check.confounders <- function(siamcat, fn.plot, verbose = 1) {
    if (verbose > 1)
        message("+ starting check.confounders")
    s.time <- proc.time()[3]
    # TODO: implement color.scheme selection as function parameter
    pdf(fn.plot, onefile = TRUE)
    label      <- label(siamcat)
    case.count <- sum(label$p.idx)
    ctrl.count <- sum(label$n.idx)

    if (verbose > 2)
        message("+++ checking group sizes")
    if (case.count > ctrl.count) {
        lgr <- label$p.idx
        smlr <- label$n.idx
        bp.labs <- c(label$p.lab, label$n.lab)
    } else {
        lgr <- label$n.idx
        smlr <- label$p.idx
        bp.labs <- c(label$n.lab, label$p.lab)
    }

    if (verbose > 2)
        message("+++ setting up color scheme")
    len.diff <- abs(case.count - ctrl.count)

    hmapcolors <- colorRampPalette(brewer.pal(10, "RdYlGn"))
    color.scheme <- rev(colorRampPalette(brewer.pal(
        brewer.pal.info["BrBG", "maxcolors"], "BrBG"))(100))

    for (m in seq_len(ncol(meta(siamcat)))) {
        mname <- gsub("[_.-]", " ", colnames(meta(siamcat))[m])
        mname <- paste(toupper(substring(mname, 1, 1)), substring(mname, 2),
         sep = "")
        if (verbose > 1)
            message(paste("+++ checking", mname, "as a potential confounder"))

        mvar <- as.numeric(unlist(meta(siamcat)[, m]))
        u.val <- unique(mvar)
        u.val <- u.val[!is.na(u.val)]
        colors <- brewer.pal(5, "Spectral")
        histcolors <- brewer.pal(9, "YlGnBu")
        if (length(u.val) == 1) {
            if (verbose > 1) {
                message("+++  skipped because all subjects have the same
                 value")
            }
        } else if (length(u.val) <= 5) {
            if (verbose > 1)
                message("++++ discreet variable, using a bar plot")

            ct <- vapply(u.val, FUN=function(x){
              return(c(sum(mvar[label$n.idx] == x, na.rm=TRUE),
                       sum(mvar[label$p.idx] == x, na.rm=TRUE)))},
              USE.NAMES = FALSE, FUN.VALUE = integer(2))

            freq <- t(ct/rowSums(ct))

            if (verbose > 2)
                message("++++ plotting barplot")
            # packages required for grid + base
            layout(matrix(c(1, 1, 2)))
            # barplot
            par(mar = c(4.1, 9.1, 4.1, 9.1))
            vps <- baseViewports()
            pushViewport(vps$figure)
            vp1 <- plotViewport()
            bar.plot <- barplot(freq, ylim = c(0, 1), main = mname, 
                names.arg = c(label$n.lab, label$p.lab),
                col = colors)
            legend(2.5, 1, legend = u.val, xpd = NA, lwd = 2, col = colors,
             inset = 0.5, bg = "grey96", cex = 0.8)
            p.val <- fisher.test(ct)$p.value
            mtext(paste("Fisher test p-value:", format(p.val, digits = 4)),
             cex = 0.6, side = 1, line = 2)
            popViewport()

            if (verbose > 2)
                message("++++ drawing contingency table")
            # contingency table
            plot.new()
            vps <- baseViewports()
            pushViewport(vps$figure)
            niceLabel <- rep(label$p.lab, length(label$label))
            names(niceLabel) <- names(label$label)
            niceLabel[label$n.idx] <- label$n.lab
            vp1 <- plotViewport()
            t <- addmargins(table(mvar, niceLabel, dnn = c(mname, "Label")))
            grid.table(t, theme = ttheme_minimal())
            popViewport()
            par(mfrow = c(1, 1), bty = "o")
        } else {
            if (verbose > 1)
                message("++++ continuous variable, using a Q-Q plot")

            layout(rbind(c(1, 2), c(3, 4)))

            if (verbose > 2)
                message("++++ panel 1/4: Q-Q plot")
            par(mar = c(4.5, 4.5, 2.5, 1.5), mgp = c(2.5, 1, 0))
            ax.int <- c(min(mvar, na.rm = TRUE), max(mvar, na.rm = TRUE))
            qqplot(mvar[label$n.idx], mvar[label$p.idx], xlim = ax.int,
             ylim = ax.int, pch = 16, cex = 0.6,
                xlab = label$n.lab, ylab = label$p.lab, main = 
                paste("Q-Q plot for", mname))
            abline(0, 1, lty = 3)
            p.val <- wilcox.test(mvar[label$n.idx], mvar[label$p.idx], 
                exact = FALSE)$p.value
            text(ax.int[1] + 0.9 * (ax.int[2] - ax.int[1]), 
                ax.int[1] + 0.1 * (ax.int[2] - ax.int[1]), cex = 0.8, 
                paste("MWW test p-value:",format(p.val, digits = 4)), pos = 2)

            if (verbose > 2)
                message("++++ panel 2/4: X histogram")
            par(mar = c(4, 2.5, 3.5, 1.5))
            hist(mvar[label$n.idx], main = label$n.lab, xlab = mname, 
                col = histcolors, breaks = seq(min(mvar,
                na.rm = TRUE), max(mvar, na.rm = TRUE), length.out = 10))
            mtext(paste("N =", length(mvar[label$n.idx])), cex = 0.6, side = 3,
             adj = 1, line = 1)

            if (verbose > 2)
                message("++++ panel 3/4: X boxplot")
            par(mar = c(2.5, 4.5, 2.5, 1.5))
            combine <- data.frame(mvar[lgr], c(mvar[smlr], rep(NA, len.diff)))


            boxplot(combine[, 1], na.omit(combine[, 2]), range=1.5, 
                use.cols = TRUE, names = bp.labs, ylab = mname, 
                main = paste("Boxplot for",
                mname), col = histcolors, outpch=NA)
            stripchart(combine, vertical = TRUE, add = TRUE,
             method = "jitter", pch = 20)

            if (verbose > 2)
                message("++++ panel 4/4: Y histogram")
            par(mar = c(4.5, 2.5, 3.5, 1.5))
            hist(mvar[label$p.idx], main = label$p.lab, xlab = mname,
             col = histcolors, breaks = seq(min(mvar,
                na.rm = TRUE), max(mvar, na.rm = TRUE), length.out = 10))
            mtext(paste("N =", length(mvar[label$p.idx])), cex = 0.6, side = 3,
             adj = 1, line = 1)
            par(mfrow = c(1, 1))
        }
    }

    e.time <- proc.time()[3]
    tmp <- dev.off()
    if (verbose > 1) {
        message(paste("+ finished check.confounders in", 
            formatC(e.time - s.time, digits=3), "s"))
    }
    if (verbose == 1) {
        message(paste("Checking metadata for confounders finished, results 
            plotted to:", fn.plot))
    }
}
