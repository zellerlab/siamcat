#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Model Interpretation Plot
#'
#' @description This function produces a plot for model interpretation
#'
#' @usage model.interpretation.plot(siamcat, fn.plot = NULL,
#' color.scheme = "BrBG", consens.thres = 0.5, heatmap.type = "zscore",
#' limits = c(-3, 3), log.n0 = 1e-06, max.show = 50, prompt=TRUE,
#' verbose = 1)
#'
#' @param siamcat object of class \link{siamcat-class}
#'
#' @param fn.plot string, filename for the pdf-plot
#'
#' @param color.scheme color scheme for the heatmap, defaults to \code{'BrBG'}
#'
#' @param consens.thres float, minimal ratio of models incorporating a feature
#' in order to include it into the heatmap, defaults to \code{0.5}
#' \strong{Note that for \code{'randomForest'} models, this cutoff specifies the
#' minimum median Gini coefficient for a feature to be included and
#' should therefore be much lower, e.g. \code{0.01}}
#'
#' @param heatmap.type string, type of the heatmap, can be either \code{'fc'}
#' or \code{'zscore'}, defaults to \code{'zscore'}
#'
#' @param limits vector, cutoff for extreme values in the heatmap,
#' defaults to \code{c(-3, 3)}
#'
#' @param log.n0 float, pseudocount to be added before log-transformation
#' of features, defaults to \code{1e-06}
#'
#' @param max.show integer, maximum number of features to be shown in the model
#' interpretation plot, defaults to 50
#'
#' @param prompt boolean, turn on/off prompting user input when not plotting
#' into a pdf-file, defaults to TRUE
#'
#' @param verbose control output: \code{0} for no output at all, \code{1}
#' for only information about progress and success, \code{2} for normal
#' level of information and \code{3} for full debug information,
#' defaults to \code{1}
#'
#' @keywords SIAMCAT model.interpretation.plot
#'
#' @return Does not return anything, but produces the model interpretation plot.
#'
#' @export
#'
#' @encoding UTF-8
#'
#' @details Produces a plot consisting of \itemize{
#' \item a barplot showing the feature weights and their robustness (i.e. in
#' what proportion of models have they been incorporated)
#' \item a heatmap showing the z-scores of the metagenomic features across
#' samples
#' \item another heatmap displaying the metadata categories (if applicable)
#' \item a boxplot displaying the poportion of weight per model that is
#' actually shown for the features that are incorporated into more than
#' \code{consens.thres} percent of the models.
#' }
#'
#' @examples
#' data(siamcat_example)
#'
#' # simple working example
#' siamcat_example <- train.model(siamcat_example, method='lasso')
#' model.interpretation.plot(siamcat_example, fn.plot='./interpretion.pdf',
#'     heatmap.type='zscore')
model.interpretation.plot <-
    function(siamcat,
        fn.plot=NULL,
        color.scheme = "BrBG",
        consens.thres = 0.5,
        heatmap.type = "zscore",
        limits = c(-3, 3),
        log.n0 = 1e-06,
        max.show = 50,
        prompt=TRUE,
        verbose = 1) {

        if (verbose > 1)
            message("+ starting model.interpretation.plot")

        s.time <- proc.time()[3]

        # ######################################################################
        # general checks
        stopifnot(length(heatmap.type) == 1)
        stopifnot(heatmap.type %in% c('zscore', 'fc', 'log'))
            
        if (is.null(model_list(siamcat, verbose=0))){
            stop("SIAMCAT object does not contain any models. Exiting...")
        }
        if (is.null(pred_matrix(siamcat, verbose=0))){
            stop("SIAMCAT object does not contain any predictions. Exiting...")
        }

        if (heatmap.type %in% c('fc', 'log')){
            if (any(orig_feat(siamcat) < 0) |
                any(colSums(orig_feat(siamcat)) > 1.01)){
                    stop("Original data should be compositional for ",
                        "heatmap.type `fc` or `log`. Exiting...")
            }
        }
        label <- label(siamcat)
        if (heatmap.type == 'fc' & label$type == "CONTINUOUS"){
            stop("Regression-type model cannot be combined with fold-change ",
                "heatmap! Exiting...")
        }


        model.type <- model_type(siamcat)
        if (model.type=='SVM'){
            stop("Interpretation heatmap not possible for kernel SVMs!")
        }
        feature.type <- feature_type(siamcat)
        feature.weights <- feature_weights(siamcat)
        weight.matrix <- weight_matrix(siamcat)

        # ######################################################################
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

        # ######################################################################
        # some color pre-processing
        if (verbose > 2)
            message("+++ preprocessing color scheme")
        if (!color.scheme %in% row.names(brewer.pal.info)) {
            warning(
                "Not a valid RColorBrewer palette name, defaulting to BrBG...\n
                See brewer.pal.info for more information about
                RColorBrewer palettes."
            )
            color.scheme <- "BrBG"
        }
        if (brewer.pal.info[color.scheme, 'category'] == 'seq'){
            color.scheme <- colorRampPalette(
                    brewer.pal(brewer.pal.info[color.scheme, "maxcolors"],
                                color.scheme))(100)
        } else {
            color.scheme <-
                rev(colorRampPalette(
                    brewer.pal(brewer.pal.info[color.scheme, "maxcolors"],
                        color.scheme))(100))
        }
        #browser()
        # ######################################################################
        # preprocess models
        if (verbose > 2)
            message("+++ preprocessing models")
        sel.idx <-
            model.interpretation.select.features(
                feature.weights = feature.weights,
                model.type = model.type,
                consens.thres = consens.thres,
                label = label,
                max.show = max.show,
                verbose = verbose
            )
        num.sel.f <- length(sel.idx)

        # ######################################################################
        # aggregate predictions and sort
        # patients by score aggregate predictions of several models if more than
        # one is given
        mean.agg.pred <- rowMeans(pred_matrix(siamcat))
        # idx to sort samples according to their class membership and prediction
        # score
        if (label$type=='BINARY'){
            srt.idx <-
                sort(label$label + mean.agg.pred, index.return = TRUE)$ix
        } else if (label$type=='CONTINUOUS'){
            srt.idx <- sort(label$label, index.return=TRUE)$ix
        }

        # ######################################################################
        # prepare heatmap
        if (verbose > 2)
            message("+++ preparing heatmap")
        if (heatmap.type == "zscore") {
            if (feature.type == 'original'){
                feat <- get.orig_feat.matrix(siamcat)
            } else if (feature.type == 'filtered') {
                feat <- get.filt_feat.matrix(siamcat)
            } else if (feature.type == 'normalized') {
                feat <- get.norm_feat.matrix(siamcat)
            }
            img.data <- model.interpretation.prepare.heatmap.zscore(
                heatmap.data = feat[sel.idx, srt.idx],
                limits = limits,
                verbose = verbose)
        } else if (heatmap.type == "fc") {
            feat <- get.orig_feat.matrix(siamcat)
            if (is.null(log.n0)) {
                warning(
                    "WARNING: Pseudo-count before log-transformation
                    not supplied! Estimating it as 5% percentile..."
                )
                log.n0 <- quantile(feat[feat != 0], 0.05)
            }
            img.data <- model.interpretation.prepare.heatmap.fc(
                heatmap.data = feat[, srt.idx],
                sel.feat = names(sel.idx),
                limits = limits,
                meta = meta(siamcat),
                label = label,
                log.n0 = log.n0,
                verbose = verbose)
        } else if (heatmap.type == 'log') {
            feat <- get.orig_feat.matrix(siamcat)
            if (is.null(log.n0)) {
                warning(
                    "WARNING: Pseudo-count before log-transformation
                    not supplied! Estimating it as 5% percentile..."
                )
                log.n0 <- quantile(feat[feat != 0], 0.05)
            }

            img.data <- model.interpretation.prepare.heatmap.log(
                heatmap.data = feat[names(sel.idx), srt.idx],
                log.n0 = log.n0,
                verbose = verbose)
            limits[1] <- log10(log.n0)
            limits[2] <- 0
        } else {
            stop("! unknown heatmap.type: ", heatmap.type)
        }

        # ######################################################################
        # start plotting model properties
        if (verbose > 2)
            message("+++ plotting model properties")
        if (!is.null(fn.plot)) {
            pdf(fn.plot, paper = "special", height = 8.27, width = 11.69,
                onefile = TRUE)
        }

        ### plot layout
        sel.f.cex <- max(0.3, 0.8 - 0.01 * num.sel.f)
        lmat <- rbind(c(1, 2, 3, 4), c(5, 6, 0, 7), c(0, 8, 0, 0))
        h_t <- 0.1
        h_m <- ifelse(is.null(meta(siamcat)), 0.8,
            max(0.5, 0.7 - 0.01 * ncol(meta(siamcat))))
        h_b <- 1 - h_t - h_m
        if (verbose >2)
            message(paste0("+++ Layout height values: ", h_t,
                ", ", h_m, ", ", h_b))
        layout(lmat, widths = c(0.14, 0.58, 0.1, 0.14),
            heights = c(h_t, h_m, h_b))
        par(oma = c(3, 4, 3, 4))

        ### header row
        ########################################################################
        # Title of Feature Weights
        if (verbose > 2)
            message("+++ plotting titles")

        par(mar = c(0, 1.1, 3.1, 1.1))
        plot(NULL, type = "n", xlim = c(-0.1, 0.1), xaxt = "n", xlab = "",
            ylim = c(-0.1, 0.1), yaxt = "n", ylab = "", bty = "n"
        )
        mtext("Feature Weights", side = 3, line = 2, at = 0.04,
            cex = 1, adj = 0.5)

        # ######################################################################
        # Title of heatmap and brackets for classes
        par(mar = c(0, 4.1, 3.1, 5.1))
        hm.label <- label$label[srt.idx]
        plot(NULL, type = "n", xlim = c(0, length(hm.label)), xaxs = "i",
            xaxt = "n", ylim = c(-0.5, 0.5), yaxs = "i", yaxt = "n",
            xlab = "", ylab = "", bty = "n")
        if (label$type=='BINARY'){
            ul <- unique(hm.label)

            for (l in seq_along(ul)) {
                idx <- which(ul[l] == hm.label)
                lines(c(idx[1] - 0.8, idx[length(idx)] - 0.2), c(0, 0))
                lines(c(idx[1] - 0.8, idx[1] - 0.8), c(-0.2, 0))
                lines(c(idx[length(idx)] - 0.2, idx[length(idx)] - 0.2),
                    c(-0.2, 0))
                h <- (idx[1] + idx[length(idx)]) / 2
                t <- gsub("_", " ",
                    names(label$info)[label$info == ul[l]])
                t <- paste(t, " (n=", length(idx), ")", sep = "")
                mtext(t, side = 3, line = -0.5, at = h, cex = 0.7, adj = 0.5)
            }
        }
        mtext("Metagenomic Features", side = 3, line = 2,
            at = length(hm.label) / 2, cex = 1, adj = 0.5)

        # ######################################################################
        # Heatmap legend
        if (verbose > 2)
            message("+++ plotting legend")
        par(mar = c(3.1, 1.1, 1.1, 1.1))
        barplot(
            as.matrix(rep(1, 100)),
            col = color.scheme,
            horiz = TRUE,
            border = NA,
            ylab = "",
            axes = FALSE
        )
        if (heatmap.type == "fc") {
            key.ticks <- seq(round(min(img.data, na.rm = TRUE), digits = 1),
                round(max(img.data, na.rm = TRUE), digits = 1),
                length.out = 7)
            key.label <- "Feature fold change over controls"
        } else if (heatmap.type == "zscore") {
            key.ticks <- seq(limits[1], limits[2], length.out = 7)
            key.label <- "Feature z-score"
        } else if (heatmap.type == 'log'){
            key.ticks <- seq(limits[1], limits[2], length.out = 7)
            key.label <- "Log10 relative abundance"
        }
        axis(side = 1, at = seq(0, 100, length.out = 7), labels = key.ticks)
        mtext(key.label, side = 3, line = 0.5, at = 50, cex = 0.7, adj = 0.5)

        # ######################################################################
        # Model header (model sensitive)
        par(mar = c(0, 6.1, 3.1, 1.1))
        plot(NULL, type = "n", xlim = c(-0.1, 0.1), xaxt = "n", xlab = "",
            ylim = c(-0.1, 0.1), yaxt = "n", ylab = "", bty = "n")
        mtext(paste0(model.type, " model"), side = 3, line = 2, at = 0.04,
            cex = 0.7, adj = 0.5)
        mtext(paste("(|W| = ", num.sel.f, ")", sep = ""), side = 3, line = 1,
            at = 0.04, cex = 0.7, adj = 0.5)

        # ######################################################################
        # Feature weights ( model sensitive)
        if (verbose > 2)
            message("+++ plotting feature weights")
        if (any(colSums(abs(weight.matrix)) == 0)){
            stop("Some models do not contain any features!")
        }
        rel.weights <- t(t(weight.matrix)/
            colSums(abs(weight.matrix), na.rm=TRUE))
        model.interpretation.feature.weights.plot(
            rel.weights = rel.weights,
            sel.idx = sel.idx,
            mod.type = model.type,
            label = label, verbose=verbose)

        # ######################################################################
        # Heatmap
        if (verbose > 2)
            message("+++ plotting heatmap")
        if (model.type != "randomForest") {
            model.interpretation.heatmap.plot(
                image.data = img.data,
                limits = limits,
                color.scheme = color.scheme,
                effect.size = rowMedians(rel.weights[sel.idx,]),
                verbose = verbose)
        } else {
            auroc.effect <- apply(img.data, 2,
                FUN = function(f) {
                    roc(predictor = f, response = label$label,
                        direction = "<", levels = label$info)$auc
                })
            bin.auroc.effect <- ifelse(auroc.effect >= 0.5, 1, 0)
            model.interpretation.heatmap.plot(
                image.data = img.data,
                limits = limits,
                color.scheme = color.scheme,
                effect.size = bin.auroc.effect,
                verbose = verbose
            )
        }

        # ######################################################################
        # Proportion of weights shown
        if (verbose > 2)
            message("+++ plotting proportion of weights shown")
        model.interpretation.proportion.of.weights.plot(
            s.idx = sel.idx,
            weights = weight.matrix,
            verbose = verbose
        )

        # ######################################################################
        # Metadata and prediction
        if (verbose > 2)
            message("+++ plotting metadata and predictions")
        model.interpretation.pred.and.meta.plot(
            prediction = mean.agg.pred[srt.idx],
            label = label,
            meta = meta(siamcat)[srt.idx,],
            verbose = verbose
        )

        if(!is.null(fn.plot)) {
            tmp <- dev.off()
        } else {
            par(par.old)
        }
        e.time <- proc.time()[3]
        if (verbose > 1)
            message(paste(
                "+ finished model.interpretation.plot in",
                formatC(e.time - s.time, digits = 3),
                "s"
            ))
        if (verbose == 1 & !is.null(fn.plot))
            message(paste(
                "Successfully plotted model interpretation plot to:",
                fn.plot
            ))
        }

# function to plot the feature weights
#' @keywords internal
model.interpretation.feature.weights.plot <-
    function(rel.weights,
        sel.idx,
        mod.type,
        label,
        verbose = 0) {
        if (verbose > 2)
            message("+ model.interpretation.feature.weights.plot")
        if (mod.type == 'randomForest'){
            relative.weights <- -rel.weights[sel.idx, ]
        } else {
            relative.weights <- rel.weights[sel.idx, ]
        }
        med <- rowMedians(relative.weights)
        low.qt <- rowQuantiles(relative.weights, probs = 0.25)
        upp.qt <- rowQuantiles(relative.weights, probs = 0.75)

        if (mod.type != "randomForest") {
            par(mar = c(0.1, 1.1, 0, 1.1))
            mi <- min(-med - (abs(low.qt - upp.qt)))
            mx <- max(-med + (abs(low.qt - upp.qt)))

            barplot(-med, horiz = TRUE, width = 1, space = 0, yaxs = "i",
                col = "gray30",
                xlim = c(-max(abs(mi), abs(mx)),
                    max(abs(mi), max(mx))),
                ylim = c(0, dim(relative.weights)[1]),
                xlab = "", ylab = "", yaxt = "n")
            x <- par("usr")
            rect(x[1], x[3], x[2], x[4], col = "lightgray")
            # grid lines
            grid(NULL, NA, lty = 2, col = "gray99")
            # plot the barplot again due to the idiotic way of how to change the
            # background color of a plot in r
            barplot(-med, horiz = TRUE, width = 1, space = 0, yaxs = "i",
                col = "gray30",
                xlim = c(-max(abs(mi), abs(mx)),
                    max(abs(mi), max(mx))),
                ylim = c(0, dim(relative.weights)[1]),
                xlab = "", ylab = "", yaxt = "n", xaxt = "n", add = TRUE)

            # error bars
            arrows(y0 = c(seq_along(med)) - 0.5, x0 = -upp.qt, x1 = -low.qt,
                angle = 90, code = 3, length = 0.04)
            mtext("median relative feat. weight", side = 1, line = 2,
                at = 0, cex = 0.7, adj = 0.5)
            # robustness indicated as percentage of models including a given
            # feature (to the right of the barplot)
            for (f in seq_along(sel.idx)) {
                t <- paste(format(
                    100 * sum(rel.weights[sel.idx[f], ] != 0) /
                        dim(rel.weights)[2],
                    digits = 1, scientific = FALSE), "%", sep = "")
                mtext(t, side = 4, line = 2.5, at = (f - 0.5),
                    cex = max(0.3, 0.8 - 0.01 * length(sel.idx)),
                    las = 2, adj = 1)
            }
            if (label$type=='BINARY'){
                # label positive/negative class
                mtext(gsub("_", " ",
                    names(which(label$info == min(label$info)))),
                    side = 2, at = floor(length(sel.idx) / 2), line = -2)
                mtext(gsub("_", " ",
                    names(which(label$info == max(label$info)))),
                    side = 4, at = floor(length(sel.idx) / 2), line = -2)
            }

            mtext("effect size", side = 3, line = 1, at = (mx / 2),
                cex = 0.7, adj = 1)
            mtext("robustness", side = 3, line = 1, at = mx, cex = 0.7,
                adj = 0)
        } else {
            par(mar = c(0.1, 1.1, 0, 1.1))
            mx <- max(-med + (abs(low.qt - upp.qt)))

            barplot(-med, horiz = TRUE, width = 1, space = 0, yaxs = "i",
                col = "gray30", xlim = c(0, mx),
                ylim = c(0, dim(relative.weights)[1]),
                xlab = "", ylab = "", yaxt = "n")
            x <- par("usr")
            rect(x[1], x[3], x[2], x[4], col = "lightgray")
            # grid lines
            grid(NULL, NA, lty = 2, col = "gray99")
            # plot the barplot again due to the idiotic way of how to change the
            # background color of a plot in r
            barplot(-med, horiz = TRUE, width = 1, space = 0, yaxs = "i",
                col = "gray30", xlim = c(0, mx),
                ylim = c(0, dim(relative.weights)[1]),
                xlab = "", ylab = "", yaxt = "n", xaxt='n', add = TRUE)

            # error bars
            arrows(y0 = c(seq_along(med)) - 0.5, x0 = -upp.qt, x1 = -low.qt,
                angle = 90, code = 3, length = 0.04)
            # labels
            mtext("median relative Gini coefficient", side = 1, line = 2,
                at = mx / 2, cex = 0.7, adj = 0.5)
            mtext("effect size", side = 3, line = 1,
                at = (mx / 2), cex = 0.7, adj = 0.5)
        }

        box(lwd = 1)
        if (verbose > 2)
            message("+ finished model.interpretation.feature.weights.plot")
    }

# function to plot the predictions and metadata
#' @keywords internal
model.interpretation.pred.and.meta.plot <-
    function(prediction, label, meta = NULL, verbose = 0) {

        if (verbose > 2)
            message("+ model.interpretation.pred.and.meta.plot")
        par(mar = c(1.1, 4.1, 0.3, 5.1))
        if (label$type == 'BINARY'){
            img.data <- as.matrix(prediction)
            colnames(img.data) <- "Classification result"
        } else {
            img.data <- cbind(as.matrix(prediction),
                            as.matrix(label$label[names(prediction)]))
            colnames(img.data) <- c('Predicted value', 'True value')
        }
        img.data.processed <- NULL
        if (!is.null(meta)) {
            img.data <- cbind(meta[, ncol(meta):1], img.data)
            ### transform any categorial column into a numeric one
            for (m in seq_len(ncol(img.data))) {
                cur.processed.data <- NULL
                temp.metadata <- img.data[,m]
                if (!all(is.numeric(temp.metadata))) {

                    temp.metadata <- factor(temp.metadata)
                    temp.metadata <- as.numeric(temp.metadata)
                }
                cur.processed.data <- (temp.metadata - min(temp.metadata,
                        na.rm = TRUE))
                if (max(temp.metadata, na.rm = TRUE) != 0) {
                    cur.processed.data <- cur.processed.data /
                        max(cur.processed.data, na.rm = TRUE)
                }
                img.data.processed <-
                    cbind(img.data.processed, cur.processed.data)
            }
        } else {
            img.data.processed <- img.data
        }

        grays <- rev(gray(seq(0, 1, length.out = 100)))
        image(
            as.matrix(img.data.processed),
            col = grays,
            xaxt = "n",
            yaxt = "n",
            xlab = "",
            ylab = "",
            bty = "n"
        )
        box(lwd = 1)
        if (label$type == 'BINARY'){
            # add deliminator between the different classes
            abline(v = length(which(label$label == min(label$info))) /
                    length(label$label),
                col = "red")
        }

        meta.cex <- max(0.3, 0.7 - 0.01 * ncol(img.data))
        for (m in seq_len(ncol(img.data))) {
            t <- colnames(img.data)[m]
            t <- gsub("\\.|_", " ", t)
            mtext(
                t,
                side = 4,
                line = 1,
                at = (m - 1) / (dim(img.data)[2] - 1),
                cex = meta.cex,
                las = 2
            )
        }
        # mark missing values
        x.increment <- 1/((nrow(img.data) - 1) * 2)
        y.increment <- 1/((ncol(img.data) - 1) * 2)
        for (m in seq_len(ncol(img.data))) {
            idx <- which(is.na(img.data[, m]))
            for (i in idx) {
                x.left <- ((i-1) * x.increment * 2) - x.increment
                x.right <- (i * x.increment * 2) - x.increment
                y.bottom <- ((m-1) * y.increment * 2) - y.increment
                y.top <-  (m * y.increment * 2) - y.increment
                rect(x.left, y.bottom, x.right, y.top, col = "yellow",
                    border=NA)
            }
        }
        if (verbose > 2)
            message("+ finished model.interpretation.pred.and.meta.plot")
    }

# function to plot the proportion of weights
#' @keywords internal
model.interpretation.proportion.of.weights.plot <-
    function(s.idx, weights, verbose = 0) {
        if (verbose > 2)
            message("+ model.interpretation.proportion.of.weights.plot")
        par(mar = c(0.1, 6.1, 0, 1.1))
        boxplot(colSums(abs(weights[s.idx,])) / colSums(abs(weights)),
            ylim = c(0, 1))
        mtext("proportion of", side = 1, line = 1, at = 1, adj = 0.5, cex = 0.7)
        mtext("weight shown", side = 1, line = 2, at = 1, adj = 0.5, cex = 0.7)
        if (verbose > 2)
            message(
                "+ finished model.interpretation.proportion.of.weights.plot")
    }

# function to plot the percentage of features
#' @keywords internal
plot.percentage.of.features.plot <-
    function(selected.weights, all.weights,
        verbose = 0) {
        if (verbose > 2)
            message("+ plot.percentage.of.features")
        par(mar = c(0.1, 6.1, 0, 1.1))
        boxplot(dim(selected.weights)[1] / colSums(all.weights != 0),
            ylim = c(0, 1))
        mtext(
            "percentage of",
            side = 1,
            line = 1,
            at = 1,
            adj = 0.5,
            cex = 0.7
        )
        mtext(
            "features shown",
            side = 1,
            line = 2,
            at = 1,
            adj = 0.5,
            cex = 0.7
        )
        if (verbose > 2)
            message("+ finished plot.percentage.of.features")
    }

# function to plot the heatmap
#' @keywords internal
model.interpretation.heatmap.plot <-
    function(image.data,
        limits,
        color.scheme,
        effect.size,
        verbose = 0) {
        if (verbose > 2)
            message("+ model.interpretation.heatmap.plot")
        par(mar = c(0.1, 4.1, 0, 5.1))
        image(image.data, zlim = limits, col = color.scheme, xaxt = "n",
            yaxt = "n", xlab = "", ylab = "", bty = "n")
        if (!is.null(effect.size)) {
            for (f in seq_len(ncol(image.data))) {
                mtext(colnames(image.data)[f], side = 4, line = 1,
                    at = (f - 1) / (dim(image.data)[2] - 1), las = 2,
                    cex = max(0.3, 0.8 - 0.01 * dim(image.data)[2]),
                    col = ifelse(effect.size[f] > 0,
                        color.scheme[1 + 4],
                        color.scheme[length(color.scheme) - 4]))
            }
        } else {
            for (f in seq_len(ncol(image.data))) {
                mtext(colnames(image.data)[f], side = 4, line = 1,
                    at = (f - 1) /
                        (dim(image.data)[2] - 1), las = 2,
                    cex = max(0.3, 0.8 - 0.01 * dim(image.data)[2]),
                    col = "black")
            }
        }
        box(lwd = 1)
        if (verbose > 2)
            message("+ finished model.interpretation.heatmap.plot")
    }

# function to prepare the data to plot fold change heatmap
#' @keywords internal
model.interpretation.prepare.heatmap.fc <-
    function(heatmap.data,
        limits,
        sel.feat,
        meta = NULL,
        label,
        log.n0,
        verbose = 0) {
        if (verbose > 2)
            message("+ model.interpretation.prepare.heatmap.fc")
        n.label <- min(label$info)
        n.idx <- which(label$label == n.label)
        if (!any(grepl("META", sel.feat))) {
            feat.log <- log10(heatmap.data[sel.feat,] + log.n0)
            img.data <- t(feat.log -
                    log10(rowMedians(heatmap.data[sel.feat, n.idx]) +
                            log.n0))
        } else {
            img.data <- matrix(NA,
                nrow = length(sel.feat),
                ncol = ncol(heatmap.data))
            row.names(img.data) <- sel.feat
            if (verbose > 2)
                message("+ Selected features:")
            for (f in sel.feat) {
                if (verbose > 2)
                    message(paste("+++", f))
                if (!grepl("META", f)) {
                    median.ctr <-
                        suppressWarnings(median(as.numeric(
                            heatmap.data[f, n.idx])))
                    img.data[f, ] <-
                        log10(heatmap.data[f, ] + log.n0) -
                        log10(median.ctr + log.n0)
                } else {
                    meta.data <- meta[, grep(
                        strsplit(f, "_")[[1]][2],
                        colnames(meta),
                        ignore.case = TRUE,
                        value = TRUE
                    )]
                    meta.data <- data.frame(meta.data)[,1]
                    # transform metadata to zscores
                    meta.data <- as.numeric(meta.data)
                    meta.data <-
                        (meta.data - mean(meta.data, na.rm = TRUE)) /
                        sd(meta.data, na.rm = TRUE)
                    names(meta.data) <- rownames(meta)
                    img.data[f, ] <-
                        meta.data[colnames(heatmap.data)]
                }
            }
            img.data <- t(img.data)
        }
        img.data[img.data < limits[1]] <- limits[1]
        img.data[img.data > limits[2]] <- limits[2]

        if (verbose > 2)
            message("+ finished model.interpretation.heatmap.plot")
        return(img.data)
    }

# function to prepare the data to plot zscore heatmap
#' @keywords internal
model.interpretation.prepare.heatmap.zscore <-
    function(heatmap.data, limits,
        verbose = 0) {
        if (verbose > 2)
            message("+ prepare.heatmap.zscore")
        # data is transposed and transformed to feature z-scores for display
        img.data <-
            (heatmap.data - rowMeans(heatmap.data)) / rowSds(heatmap.data)
        img.data[img.data < limits[1]] <- limits[1]
        img.data[img.data > limits[2]] <- limits[2]
        if (verbose > 2)
            message("+ finished prepare.heatmap.zscore")
        return(t(img.data))
    }
#' @keywords internal
model.interpretation.prepare.heatmap.log <-
    function(heatmap.data, log.n0, verbose = 0) {
        if (verbose > 2)
            message("+ prepare.heatmap.log")
        # data is transposed and transformed to feature z-scores for display
        img.data <- log10(heatmap.data + log.n0)
        if (verbose > 2)
            message("+ finished prepare.heatmap.log")
        return(t(img.data))
    }

# function to select the features to plot on the heatmap
#' @keywords internal
model.interpretation.select.features <-
    function(feature.weights,
        model.type,
        consens.thres,
        label,
        max.show,
        verbose = 0) {

        if (verbose > 2) message("+ model.interpretation.select.features")
        # for linear models, select those that have been selected more than
        # consens.thres percent of the models
        if (model.type != "randomForest") {
            sel.idx <- which(feature.weights$percentage > consens.thres)
            names(sel.idx) <- rownames(feature.weights)[sel.idx]
            # normalize by model size and order features by
            #   relative model weight
            median.sorted.features <-
                sort(feature.weights$median.rel.weight[sel.idx],
                    decreasing = TRUE,
                    index.return = TRUE)
            # restrict to plot at maximum fifty features
            if (length(sel.idx) > max.show) {
                warning(paste0("WARNING: restricting amount of features",
                    " to be plotted to ", max.show))
                median.sorted.features.abs <- sort(
                    abs(feature.weights$median.rel.weight),
                    decreasing = TRUE,
                    index.return = TRUE)
                idx <- head(median.sorted.features.abs$ix, n = max.show)
                median.sorted.features <- sort(
                    feature.weights$mean.rel.weight[idx],
                    decreasing = TRUE,
                    index.return = TRUE)
                sel.idx <- idx[median.sorted.features$ix]
            } else {
                sel.idx <- sel.idx[median.sorted.features$ix]
            }
        } else {
        # for Random Forest, caluclate relative median feature weights and sort
        # by auroc as effect measure
            median.sorted.features <-
                sort(feature.weights$median.rel.weight,
                    decreasing = FALSE,
                    index.return = TRUE)
            # take the feature with median higher than consens.threshold
            sel.idx <-
                median.sorted.features$ix[which(median.sorted.features$x >=
                        consens.thres)]
            names(sel.idx) <- rownames(feature.weights)[sel.idx]

            if (length(sel.idx) > max.show) {
                sel.idx <- tail(sel.idx, n = max.show)
            }
        }

        if (verbose > 2)
            message(paste("+++ generating plot for a model with",
                length(sel.idx), "selected features"))
        if (length(sel.idx)==0) stop("No features were selected for plotting!")
        if (length(sel.idx)==1)
            stop("Not enough features were selected for plotting!")
        if (verbose > 2)
            message("+ finished model.interpretation.select.features")
        return(sel.idx)
    }
