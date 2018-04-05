### SIAMCAT - Statistical Inference of Associations between Microbial Communities And host phenoTypes R flavor EMBL
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
        lgr <- siamcat@label@n.idx
        smlr <- siamcat@label@p.idx
        bp.labs <- c(siamcat@label@n.lab, siamcat@label@p.lab)
    }
    
    if (verbose > 2) 
        cat("+++ setting up color scheme\n")
    len.diff <- abs(case.count - ctrl.count)
    hmap <- data.frame()
    hmapcolors <- colorRampPalette(brewer.pal(10, "RdYlGn"))
    color.scheme <- rev(colorRampPalette(brewer.pal(brewer.pal.info["BrBG", "maxcolors"], "BrBG"))(100))
    
    for (m in 1:ncol(siamcat@phyloseq@sam_data)) {
        mname <- gsub("[_.-]", " ", colnames(siamcat@phyloseq@sam_data)[m])
        mname <- paste(toupper(substring(mname, 1, 1)), substring(mname, 2), sep = "")
        if (verbose > 1) 
            cat("+++ checking", mname, "as a potential confounder\n")
        
        mvar <- as.numeric(unlist(siamcat@phyloseq@sam_data[, m]))
        u.val <- unique(mvar)
        u.val <- u.val[!is.na(u.val)]
        colors <- brewer.pal(5, "Spectral")
        histcolors <- brewer.pal(9, "YlGnBu")
        if (length(u.val) == 1) {
            if (verbose > 1) {
                cat("+++  skipped because all subjects have the same value\n")
            }
        } else if (length(u.val) <= 5) {
            if (verbose > 1) 
                cat("++++ discreet variable, using a bar plot\n")
            
            ct <- matrix(NA, nrow = 2, ncol = length(u.val))
            nms <- c()
            for (i in 1:length(u.val)) {
                ct[1, i] = sum(mvar[siamcat@label@n.idx] == u.val[i], na.rm = TRUE)  # ctr
                ct[2, i] = sum(mvar[siamcat@label@p.idx] == u.val[i], na.rm = TRUE)  # cases
                nms <- c(nms, paste(mname, u.val[i]))
            }
            freq <- t(ct)
            temp <- freq
            rownames(temp) <- nms
            hmap <- rbind(hmap, temp)
            for (i in 1:dim(freq)[2]) {
                freq[, i] <- freq[, i]/sum(freq[, i])
            }
            
            if (verbose > 2) 
                cat("++++ plotting barplot\n")
            # packages required for grid + base
            layout(matrix(c(1, 1, 2)))
            # barplot
            par(mar = c(4.1, 9.1, 4.1, 9.1))
            vps <- baseViewports()
            pushViewport(vps$figure)
            vp1 <- plotViewport()
            bar.plot <- barplot(freq, ylim = c(0, 1), main = mname, names.arg = c(siamcat@label@n.lab, siamcat@label@p.lab), 
                col = colors)
            legend(2.5, 1, legend = u.val, xpd = NA, lwd = 2, col = colors, inset = 0.5, bg = "grey96", cex = 0.8)
            p.val <- fisher.test(ct)$p.value
            mtext(paste("Fisher test p-value:", format(p.val, digits = 4)), cex = 0.6, side = 1, line = 2)
            popViewport()
            
            if (verbose > 2) 
                cat("++++ drawing contingency table\n")
            # contingency table
            plot.new()
            vps <- baseViewports()
            pushViewport(vps$figure)
            niceLabel <- rep(siamcat@label@p.lab, length(siamcat@label@label))
            names(niceLabel) <- names(siamcat@label@label)
            niceLabel[siamcat@label@n.idx] <- siamcat@label@n.lab
            vp1 <- plotViewport()
            t <- addmargins(table(mvar, niceLabel, dnn = c(mname, "Label")))
            grid.table(t, theme = ttheme_minimal())
            popViewport()
            par(mfrow = c(1, 1), bty = "o")
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
            # combine <- list(mvar[siamcat@label@n.idx], mvar[siamcat@label@p.idx])
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
    }
    
    # ## HEATMAP / SPECTRAL PLOT... TBD if(0){ # transform heatmap values colnames(hmap) <- c(siamcat@label@n.lab,
    # siamcat@label@p.lab) hmap.ratios <- as.matrix(log(hmap[,2]/hmap[,1])) -
    # log(sum(siamcat@label@p.idx)/sum(siamcat@label@p.idx)) rownames(hmap.ratios) <- as.matrix(rownames(hmap)) # plot
    # as labeled barplot...  par(mar=c(20, 5, 10, 5)) grad.len <- 1000 barplot(as.matrix(rep(1,grad.len)), col =
    # viridisLite::viridis(grad.len), horiz=TRUE, border=0, ylab='', axes=FALSE) # make sure to add this to the
    # Description file once it is needed bp.bound <- max(abs(range(hmap.ratios[is.finite(hmap.ratios)]))) bp.ticks <-
    # seq(-round(bp.bound), round(bp.bound), length.out=7) bp.label <- paste('Relative ratios for metadata,',
    # siamcat@label@p.lab,':',siamcat@label@n.lab) axis(1, at=seq(0,grad.len,length.out=7), labels=bp.ticks,
    # cex.axis=0.7) # draw lines to represent each metadata category for (val in hmap.ratios) {
    # abline(v=(val*grad.len/6) + (grad.len/2)) } } tmp <- dev.off()
    e.time <- proc.time()[3]
    if (verbose > 1) {
        cat("+ finished check.confounders in", e.time - s.time, "s\n")
    }
    if (verbose == 1) {
        cat("Checking metadata for confounders finished, results plotted to:", fn.plot, "\n")
    }
}
