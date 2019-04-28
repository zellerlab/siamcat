#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Check for potential confounders in the metadata
#' @description Checks potential confounders in the metadata and produces
#'     some visualizations
#' @usage check.confounders(siamcat, fn.plot, meta.in = NULL, verbose = 1)
#' @param siamcat an object of class \link{siamcat-class}
#' @param fn.plot string, filename for the pdf-plot
#' @param meta.in vector, specific metadata variable names to analyze,
#'     defaults to NULL (all metadata variables will be analyzed)
#' @param verbose integer, control output: \code{0} for no output at all,
#'     \code{1} for only information about progress and success, \code{2} for
#'     normal level of information and \code{3} for full debug information,
#'     defaults to \code{1}
#' @keywords SIAMCAT check.confounders
#' @details This function checks for associations between class labels and
#'     potential confounders (e.g. Age, Sex, or BMI) that are present in the
#'     metadata. Statistical testing is performed with Fisher's exact test or
#'     Wilcoxon test, while associations are visualized either as barplot or
#'     Q-Q plot, depending on the type of metadata.
#'
#'     Additionally, it evaluates associations among metadata variables using
#'     conditional entropy and associations with the label using generalized
#'     linear models, producing a correlation heatmap and appropriate
#'     quantitative barplots, respectively.
#' @export
#' @return Does not return anything, but outputs plots to specified pdf file
#' @examples
#' # Example data
#' data(siamcat_example)
#'
#' # Simple working example
#' check.confounders(siamcat_example, './conf_plot.pdf')

check.confounders <- function(siamcat, fn.plot, meta.in = NULL, verbose = 1) {

    pdf(fn.plot, paper = 'special', height = 8.27, width = 11.69)
    if (verbose > 1) message("+ starting check.confounders")
    s.time <- proc.time()[3]
    label <- label(siamcat)
    meta <- meta(siamcat)
    if (is.null(meta)) {
        stop('SIAMCAT object does not contain any metadata.\nExiting...')
    }
    meta <- factorize.metadata(meta) # creates data.frame

    # check validity of input metadata conditions
    if (!is.null(meta.in)){
    if (!all(meta.in %in% colnames(meta))){
        meta.in <- meta.in[which(meta.in %in% colnames(meta))]
        warning(paste0("Some specified metadata were not in metadata file.\n",
                        "Continuing with: ",  paste(c(meta.in), collapse=" ")))
    }
        meta <- meta[,meta.in]
    }
    if (ncol(meta) > 10){
        warning(paste0("The recommended number of metadata variables is 10.\n",
                    "Please be aware that some visualizations may not work."))
    }

    # remove nested variables
    indep <- vapply(colnames(meta), FUN=function(x) {
        return(condentropy(label$label, discretize(meta[,x])))},
        FUN.VALUE = numeric(1))
    if ((verbose >= 1) & (length(names(which(indep == 0))) > 0)){
        message("++ metadata variables:\n\t",
                paste(c(names(which(indep == 0))), collapse = " & "),
                "\n++ are nested inside the label and ",
                "have been removed from this analysis")
    }
    meta <- meta[,names(which(indep != 0))]
    # remove metavariables with less than 2 levels
    n.levels <- vapply(meta,
        FUN = function(x){length(unique(x))},
        FUN.VALUE = integer(1))
    if (any(n.levels < 2)){
        s.name <- names(which(n.levels < 2))
        if (verbose >= 1){
            message("++ remove metadata variables, since all ",
                "subjects have the same value\n\t", s.name)
        }
        meta <- meta[,which(n.levels > 1)]
    }

    # FIRST PLOT - conditional entropies for metadata variables
    if (verbose > 1)
        message("+++ plotting conditional entropies for metadata variables")
    confounders.corrplot(meta, label)

    # SECOND PLOT - glm regression coefficients + significance + AU-ROCs
    layout(matrix(c(1, 2, 3), nrow = 1, ncol = 3),
            widths = c(0.5, 0.3, 0.2), heights = c(1, 1, 1))
    if (verbose > 2)
        message("+++ building logistic regression classifiers for metadata")
    glm.data <- confounders.build.glms(meta, label)
    if (verbose > 2)
        message("+++ plotting regression coefficients")
    confounders.glm.reg.coef.plot(glm.data)
    if (verbose > 2)
        message("+++ plotting regression coefficient significance")
    confounders.glm.reg.pval.plot(glm.data)
    if (verbose > 2)
        message("+++ plotting au-roc values")
    confounders.glm.auroc.plot(glm.data)

    # THIRD PLOT(S) - original confounder check descriptive stat plots
    confounders.descriptive.plots(meta(siamcat)[,colnames(meta)],
        label, verbose)
    dev.off()

    e.time <- proc.time()[3]
    if (verbose > 1) {
        message(paste("+ finished check.confounders in",
                        formatC(e.time - s.time, digits = 3), "s"))
    }
    if (verbose == 1) {
        message(paste("Finished checking metadata for confounders,",
            "results plotted to:", fn.plot))
    }
}

#'@keywords internal
confounders.corrplot <- function(meta, label) {

    meta.temp <- data.frame(meta, Label = label$label)
    entropies <- matrix(NA, nrow=ncol(meta)+1, ncol=ncol(meta)+1,
                        dimnames=list(c(colnames(meta), 'Label'),
                                    c(colnames(meta), 'Label')))

    for (n in seq_along(colnames(meta.temp))) {
        entropies[,n] <- vapply(colnames(meta.temp), FUN=function(x) {
            idx <- intersect(which(!is.na(meta.temp[,n]) == TRUE),
                            which(!is.na(meta.temp[,x]) == TRUE))
            return(condentropy(discretize(meta.temp[,x][idx]),
                            discretize(meta.temp[idx,n])))},
            FUN.VALUE = numeric(1))}

    fix.names <- vapply(colnames(entropies), FUN=function(x) {
        x <- tolower(x)
        return(gsub('[_.-]', ' ', paste(toupper(substring(x, 1, 1)),
                                        substring(x, 2), sep = "")))},
        FUN.VALUE = character(1))
    colnames(entropies) <- fix.names
    rownames(entropies) <- fix.names

    col = c(rev(colorRampPalette(brewer.pal(9, 'Reds'))(100)),
                rev(colorRampPalette(brewer.pal(9, 'Blues'))(100)))

    corrplot(entropies,
            is.corr = FALSE, order = 'FPC', method = 'color',
            pch.cex = 0.9,
            cl.lim = c(min(entropies), ceiling(max(entropies))),
            tl.col = 'black', tl.srt = 45,
            col=col,
            number.cex = 0.7, number.digits = 2, addCoef.col = 'black',
            mar = c(3.1, 2.1, 5.1, 4.1))
    title(main = 'Conditional Entropies\n H(Row | Col)', cex = 0.7)

}

#'@keywords internal
confounders.build.glms <- function(meta, label) {

    reg.coef <- vector()
    reg.ci   <- list()
    reg.pval <- vector()
    rocs     <- list()
    aucs     <- vector()
    rm       <- vector()

    for (m in seq_along(colnames(meta))) {
        d <- data.frame(x = as.ordered(meta[,m]),
                        y = factor(label$label, levels = c(-1,1),
                                labels = c(0,1)))

        model <- glm(y ~ x, family = binomial("logit"), data = d)

        # check glms; conditional entropy should catch
        if (model$converged == TRUE) {
            reg.coef[m] <- model$coefficients[2]
            reg.ci[[m]] <- confint(profile(model))[2,]
            reg.pval[m] <- coef(summary(model))[2,4]
            rocs[[m]]   <- roc(y ~ x, data=d, direction='<', ci=TRUE, auc=TRUE)
            aucs[m]     <- as.numeric(rocs[[m]]$auc)}
        else {rm <- c(rm, colnames(meta)[m])}}

    idx <- !(colnames(meta) %in% rm)
    names(reg.ci)   <- colnames(meta)
    names(reg.coef) <- colnames(meta)
    names(reg.pval) <- colnames(meta)
    names(rocs)     <- colnames(meta)
    names(aucs)     <- colnames(meta)

    # order on plot(s) by regression significance
    order <- names(sort(reg.pval, decreasing=TRUE, na.last=FALSE))

    return(list("reg.ci" = reg.ci[which(lapply(reg.ci, is.null) == FALSE)],
                "reg.coef" = reg.coef[!(names(reg.coef) %in% rm)],
                "reg.pval" = reg.pval[!(names(reg.pval) %in% rm)],
                "rocs" = rocs[which(lapply(rocs, is.null) == FALSE)],
                "aucs" = aucs[!(names(aucs) %in% rm)],
                "colors" = rep('slategrey', length(colnames(meta))),
                "order" = order[!(order %in% rm) == TRUE]))
}

#'@keywords internal
confounders.glm.reg.coef.plot <- function(glm.data) {

    order <- glm.data$order
    x.min <- min(unlist(lapply(glm.data$reg.ci, min)), na.rm=TRUE, finite=TRUE)
    x.max <- max(unlist(lapply(glm.data$reg.ci, max)), na.rm=TRUE, finite=TRUE)
    margins <- c(floor(x.min), ceiling(x.max))
    x.ticks <- seq(margins[1], margins[2])
    if (length(x.ticks) < 10){
        x.ticks <- seq(margins[1], margins[2], by=0.5)
    }
    x.tick.labels <- formatC(x.ticks, digits=2)
    y.ticks <- seq(1, length(order))

    text.names <- vapply(order, FUN=function(x) {
        x <- tolower(x)
        return(gsub('[_.-]', ' ',
                        paste(toupper(substring(x, 1, 1)),
                        substring(x, 2), sep = "")))
        }, FUN.VALUE = character(1))

    par(mar = c(5.1, 10.1, 4.1, 1.1))
    plot(NULL, xlab = '', ylab = '', xaxs = 'i', yaxs = 'i', axes = FALSE,
            xlim = c(margins[1], margins[2]),
            ylim = c(0.5, length(order) + 0.5),
            type = 'n')
    abline(v = x.ticks, lty = 3, col = 'lightgrey')
    abline(v = 0, lty = 3, col = 'lightgrey')
    axis(side = 1, at = x.ticks, labels = x.tick.labels, cex.axis = 0.9)
    axis(side = 2, at = y.ticks, labels = text.names, las=2)
    title(main = 'Single Covariate Logistic Regression',
            xlab = 'Coefficients', cex=0.9)

    for (i in seq_along(order)) {
        if (!is.na(glm.data$reg.coef[order[i]])) {
            segments(x0 = glm.data$reg.ci[[order[i]]][1],
                    x1 = glm.data$reg.ci[[order[i]]][2],
                    y0 = i, col = 'darkgrey', lwd = 1.5)
            points(glm.data$reg.coef[order[i]], i, pch = 18,
                    col=glm.data$colors[i])
            points(glm.data$reg.coef[order[i]], i, pch = 5,
                    col=glm.data$colors[i])}}
}

#'@keywords internal
confounders.glm.reg.pval.plot <- function(glm.data) {

    order <- glm.data$order
    pvals.log <- -log10(glm.data$reg.pval[order])
    x.max <- max(ceiling(abs(range(pvals.log, na.rm=TRUE, finite=TRUE))))
    x.ticks <- seq(0, x.max)
    x.tick.labels <- parse(text = paste0('10^',
        signif(x.ticks*-1, digits = 2)))

    par(mar = c(5.1, 2.1, 4.1, 1.1))
    plot(NULL, xlab = '', ylab = '', xaxs = 'i', yaxs = 'i', axes = FALSE,
        xlim = c(0,x.max), ylim = c(0.5, length(order) + 0.5),
        type = 'n')
    abline(v = x.ticks, lty = 3, col = 'lightgrey')
    axis(side = 1, at = x.ticks, labels = x.tick.labels, cex.axis = 0.9)
    title(main = 'Coefficient Significance',
        xlab = 'P Value', cex=0.9)

    for (i in seq_along(pvals.log)) {
        if (!is.na(pvals.log[i])) {
            rect(0, i-0.1, pvals.log[i], i+0.1, col=glm.data$colors[i])}}
    abline(v = -log10(0.5), lty = 1, col = 'red')
}

#'@keywords internal
confounders.glm.auroc.plot <- function(glm.data) {

    x.ticks <- seq(0, 1, length.out = 5)
    x.tick.labels <- formatC(x.ticks, digits = 2)
    order <- glm.data$order

    par(mar = c(5.1, 2.1, 4.1, 4.1))
    plot(NULL, xlab = '', ylab = '', xaxs = 'i', yaxs = 'i', axes = FALSE,
        xlim = c(0,1), ylim = c(0.5, length(order) + 0.5), type = 'n')
    abline(v = x.ticks, lty = 3, col = 'lightgrey')
    abline(v = 0.5, lty = 1, col = 'darkgrey')
    axis(side = 1, at = x.ticks, labels = x.tick.labels, cex.axis = 0.9)
    title(main = 'ROC Analysis', xlab = 'AU-ROC', cex=0.9)

    for (i in seq_along(order)) {
        if (!is.na(glm.data$aucs[order[i]])) {
            segments(x0 = glm.data$rocs[[order[i]]]$ci[1],
                    x1 = glm.data$rocs[[order[i]]]$ci[3],
                    y0 = i, col = 'lightgrey', lwd = 1.5)
            points(glm.data$rocs[[order[i]]]$ci[2], i, pch = 18,
                    col=glm.data$colors[i])
            points(glm.data$rocs[[order[i]]]$ci[2], i, pch = 5,
                    col=glm.data$colors[i])}}
}

#'@keywords internal
confounders.descriptive.plots <- function(meta, label, verbose) {

    cases <- which(label$label == max(label$info))
    controls <- which(label$label == min(label$info))
    p.lab <- names(which(label$info == max(label$info)))
    n.lab <- names(which(label$info == min(label$info)))
    var.level.names <- get.names(meta) # for contingency plots & legends

    for (m in seq_along(meta)) {
        mname <- gsub("[_.-]", " ", colnames(meta)[m])
        mname <- paste(toupper(substring(mname, 1, 1)), substring(mname, 2),
            sep = "")
        if (verbose > 1)
            message(paste("+++ checking",mname,"as a potential confounder"))
        mvar <- meta[[m]]
        if (is.character(mvar)) mvar <- as.factor(mvar)
        mvar <- as.numeric(mvar)
        names(mvar) <- rownames(meta)
        u.val <- unique(mvar)[!is.na(unique(mvar))]
        colors <- brewer.pal(6, "Spectral")
        histcolors <- brewer.pal(9, "YlGnBu")

        if (length(u.val) == 1) {
            if (verbose > 1) {
                message("+++ skipped because all subjects have the",
                    "same value")}}
        else if (length(u.val) <= 6) {
            if (verbose > 1) message("++++ discrete variable, using a bar plot")

            # create contingency table
            ct <- vapply(u.val, FUN = function(x) {
                # deal with NAs...?
                return(c(length(intersect(which(mvar == x), controls)),
                        length(intersect(which(mvar == x), cases))))},
                USE.NAMES = FALSE,
                FUN.VALUE = integer(2))

            freq <- t(ct / rowSums(ct))
            mvar <- factor(mvar, levels = unique(na.omit(mvar)),
                            labels = var.level.names[[m]])

            if (verbose > 2)
                message("++++ plotting barplot")

            layout(matrix(c(1, 1, 2)))

            # barplot
            par(mar = c(4.1, 9.1, 4.1, 9.1))
            vps <- baseViewports()
            pushViewport(vps$figure)
            vp1 <- plotViewport()
            bar.plot <- barplot(freq, ylim = c(0, 1), main = mname,
                                names.arg = names(label$info), col = colors)
            legend(2.5, 1, legend = var.level.names[[m]],
                    xpd = NA, lwd = 2, col = colors,
                    inset = 0.5, bg = "grey96", cex = 0.8)
            ifelse(length(u.val) > 4,
                    p.val <- fisher.test(ct, simulate.p.value = TRUE,
                                        B = 2000)$p.value,
                    p.val <- fisher.test(ct)$p.value)
            mtext(
                paste("Fisher Test P Value:", format(p.val, digits = 4)),
                cex = 0.8, side = 1, line = 2)
            popViewport()

            if (verbose > 2)
                message("++++ drawing contingency table")

            # contingency table
            plot.new()
            vps <- baseViewports()
            pushViewport(vps$figure)
            niceLabel <- factor(label$label, levels = label$info,
                                labels = names(label$info))
            vp1 <- plotViewport()
            t <- addmargins(table(mvar, niceLabel, dnn = c(mname, "Label")))
            grid.table(t, theme = ttheme_minimal())
            popViewport()
            par(mfrow = c(1, 1), bty = "o")}
        else {
            if (verbose > 1)
                message("++++ continuous variable, using a Q-Q plot")

            layout(rbind(c(1, 2), c(3, 4)))

            if (verbose > 2)
                message("++++ panel 1/4: Q-Q plot")
            par(mar = c(4.5, 4.5, 2.5, 1.5), mgp = c(2.5, 1, 0))
            ax.int <- c(min(mvar, na.rm = TRUE), max(mvar, na.rm = TRUE))
            qqplot(mvar[controls], mvar[cases], xlim = ax.int,
                    ylim = ax.int, pch = 16, cex = 0.6, xlab = n.lab,
                    ylab = p.lab, main = paste("Q-Q plot for", mname))
            abline(0, 1, lty = 3)
            p.val <- wilcox.test(mvar[controls], mvar[cases],
                                exact = FALSE)$p.value
            text(ax.int[1] + 0.9 * (ax.int[2] - ax.int[1]),
                ax.int[1] + 0.1 * (ax.int[2] - ax.int[1]),
                cex = 0.8,paste("MWW Test P Value:",
                                format(p.val, digits = 4)),
                pos = 2)

            if (verbose > 2)
                message("++++ panel 2/4: X histogram")
            par(mar = c(4, 2.5, 3.5, 1.5))
            hist(mvar[controls], main = n.lab, xlab = mname,
                col = histcolors, breaks = seq(min(mvar, na.rm = TRUE),
                                                max(mvar, na.rm = TRUE),
                                                length.out = 10))
            mtext(paste("N =", length(mvar[controls])), cex = 0.6,
                    side = 3, adj = 1, line = 1)

            if (verbose > 2)
                message("++++ panel 3/4: X boxplot")
            par(mar = c(2.5, 4.5, 2.5, 1.5))
            boxplot(mvar ~ label$label, range = 1.5,
                    use.cols = TRUE, names = names(label$info),
                    ylab = mname, main = paste("Boxplot for", mname),
                    col = histcolors, outpch = NA)
            stripchart(mvar ~ label$label, vertical = TRUE, add = TRUE,
                    method = "jitter", pch = 20)

            if (verbose > 2)
                message("++++ panel 4/4: Y histogram")
            par(mar = c(4.5, 2.5, 3.5, 1.5))
            hist(mvar[cases], main = p.lab, xlab = mname,
                    col = histcolors, breaks = seq(min(mvar, na.rm = TRUE),
                                                max(mvar, na.rm = TRUE),
                                                length.out = 10))
            mtext(paste("N =", length(mvar[cases])), cex = 0.6,
                    side = 3, adj = 1, line = 1)
            par(mfrow = c(1, 1))
        }
            }
}

#'@keywords internal
get.names <- function(meta) {
    temp <- lapply(meta, FUN=function(x) {
        if (is.numeric(x) | is.character(x)) {
            if (length(unique(x)) <= 6) {return(levels(as.factor(x)))}
            else {return(NULL)}}
        else if (is.factor(x)) {return(levels(x))}
        else {return(levels(x))}})
}

#'@keywords internal
factorize.metadata <- function(meta) {

    if ('BMI' %in% toupper(colnames(meta))) {
        idx <- match('BMI', toupper(colnames(meta)))
        meta[,idx] <- factorize.bmi(meta[,idx])}

    factorized <- as.data.frame(lapply(meta, FUN=function(x) {
        if (is.numeric(x) & (length(unique(x)) > 5)) {
            quart <- quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE)
            temp <- cut(x, unique(quart), include.lowest = TRUE)
            return(factor(temp, labels = seq_along(levels(temp))))}
        else (return(as.factor(x)))}))
    rownames(factorized) <- rownames(meta)
    return(factorized)
}

#'@keywords internal
factorize.bmi <- function(bmi) {
    # ranges taken from CDC
    # https://www.cdc.gov/healthyweight/assessing/bmi/adult_bmi/index.html

    if (!is.matrix(bmi)) bmi <- as.matrix(bmi)
    temp <- vapply(bmi, FUN=function(x) {
        if (is.na(x)) {return(as.character(NA))}
        else if (x < 18.5) {return("Underweight")}
        else if ((x >= 18.5) & (x <= 24.9)) {return("Healthy")}
        else if ((x > 24.9) & (x <= 29.9)) {return("Overweight")}
        else if (x > 29.9) {return("Obese")}},
        FUN.VALUE = character(1), USE.NAMES = TRUE)
    #names(temp) <- rownames(bmi)
    return(as.factor(temp))
}
