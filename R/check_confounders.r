#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between Microbial Communities And host phenoTypes EMBL Heidelberg
### 2012-2018 GNU GPL 3.0

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
        cat("+ starting check.confounders\n")
    s.time <- proc.time()[3]
    # TODO: implement color.scheme selection as function parameter
    pdf(fn.plot, onefile = TRUE)
    case.count <- length(siamcat@label@label[siamcat@label@p.idx])
    ctrl.count <- length(siamcat@label@label[siamcat@label@n.idx])
    
    if (verbose > 2) 
        cat("+++ checking group sizes\n")
    if (case.count > ctrl.count) {
        lgr <- siamcat@label@p.idx
        smlr <- siamcat@label@n.idx
        bp.labs <- c(siamcat@label@p.lab, siamcat@label@n.lab)
    } else {
        if (verbose > 1) 
            cat("++++ continuous variable, using a Q-Q plot\n")
        # discretize continuous variable; split at median for now
        dct <- matrix(NA, nrow = 2, ncol = 2)
        dct[1, ] <- c(sum(mvar[siamcat@label@n.idx] <= median(mvar, na.rm = TRUE), na.rm = TRUE), sum(mvar[siamcat@label@p.idx] <= 
            median(mvar, na.rm = TRUE), na.rm = TRUE))
        dct[2, ] <- c(sum(mvar[siamcat@label@n.idx] > median(mvar, na.rm = TRUE), na.rm = TRUE), sum(mvar[siamcat@label@p.idx] > 
            median(mvar, na.rm = TRUE), na.rm = TRUE))
        rownames(dct) <- c(paste(mname, "<= med"), paste(mname, "> med"))
        hmap <- rbind(hmap, dct)
        layout(rbind(c(1, 2), c(3, 4)))
        
        if (verbose > 2) 
            cat("++++ panel 1/4: Q-Q plot\n")
        par(mar = c(4.5, 4.5, 2.5, 1.5), mgp = c(2.5, 1, 0))
        ax.int <- c(min(mvar, na.rm = TRUE), max(mvar, na.rm = TRUE))
        qqplot(mvar[siamcat@label@n.idx], mvar[siamcat@label@p.idx], xlim = ax.int, ylim = ax.int, pch = 16, cex = 0.6, 
            xlab = siamcat@label@n.lab, ylab = siamcat@label@p.lab, main = paste("Q-Q plot for", mname))
        abline(0, 1, lty = 3)
        p.val <- wilcox.test(mvar[siamcat@label@n.idx], mvar[siamcat@label@p.idx], exact = FALSE)$p.value
        text(ax.int[1] + 0.9 * (ax.int[2] - ax.int[1]), ax.int[1] + 0.1 * (ax.int[2] - ax.int[1]), cex = 0.8, paste("MWW test p-value:", 
            format(p.val, digits = 4)), pos = 2)
        
        if (verbose > 2) 
            cat("++++ panel 2/4: X histogram\n")
        par(mar = c(4, 2.5, 3.5, 1.5))
        hist(mvar[siamcat@label@n.idx], main = siamcat@label@n.lab, xlab = mname, col = histcolors, breaks = seq(min(mvar, 
            na.rm = TRUE), max(mvar, na.rm = TRUE), length.out = 10))
        mtext(paste("N =", length(mvar[siamcat@label@n.idx])), cex = 0.6, side = 3, adj = 1, line = 1)
        
        if (verbose > 2) 
            cat("++++ panel 3/4: X boxplot\n")
        par(mar = c(2.5, 4.5, 2.5, 1.5))
        combine <- data.frame(mvar[lgr], c(mvar[smlr], rep(NA, len.diff)))
        
        
        boxplot(combine[, 1], na.omit(combine[, 2]), use.cols = TRUE, names = bp.labs, ylab = mname, main = paste("Boxplot for", 
            mname), col = histcolors)
        stripchart(combine, vertical = TRUE, add = TRUE, method = "jitter", pch = 20)
        
        if (verbose > 2) 
            cat("++++ panel 4/4: Y histogram\n")
        par(mar = c(4.5, 2.5, 3.5, 1.5))
        hist(mvar[siamcat@label@p.idx], main = siamcat@label@p.lab, xlab = mname, col = histcolors, breaks = seq(min(mvar, 
            na.rm = TRUE), max(mvar, na.rm = TRUE), length.out = 10))
        mtext(paste("N =", length(mvar[siamcat@label@p.idx])), cex = 0.6, side = 3, adj = 1, line = 1)
        par(mfrow = c(1, 1))
    }
    
    e.time <- proc.time()[3]
    if (verbose > 1) {
        cat("+ finished check.confounders in", e.time - s.time, "s\n")
    }
    if (verbose == 1) {
        cat("Checking metadata for confounders finished, results plotted to:", fn.plot, "\n")
    }
}
