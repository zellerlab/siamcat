#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Make predictions on a test set
#'
#' @description This function takes a \link{siamcat-class}-object containing
#'     a model trained by \link{train.model} and performs predictions on
#'     a given test-set.
#'
#' @usage make.predictions(siamcat, siamcat.holdout = NULL,
#'     normalize.holdout = TRUE, verbose = 1)
#'
#' @param siamcat object of class \link{siamcat-class}
#'
#' @param siamcat.holdout optional, object of class \link{siamcat-class} on
#'     which to make predictions, defaults to \code{NULL}
#'
#' @param normalize.holdout boolean, should the holdout features be normalized
#'     with a frozen normalization (see \link{normalize.features}) using the
#'     normalization parameters in \code{siamcat}?, defaults to \code{TRUE}
#'
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'
#'     for only information about progress and success, \code{2} for normal
#'     level of information and \code{3} for full debug information,
#'     defaults to \code{1}
#'
#' @export
#'
#' @keywords SIAMCAT make.predictions
#'
#' @return object of class \link{siamcat-class} with the slot \code{pred_matrix}
#'     filled or a matrix containing the predictions for the holdout set
#'
#' @details This functions uses the model in the \code{model_list}-slot of the
#'     \code{siamcat} object to make predictions on a given test set. The
#'     test set can either consist of the test instances in the cross-
#'     validation, saved in the \code{data_split}-slot of the same
#'     \code{siamcat} object, or a completely external feature set, given in
#'     the form of another \code{siamcat} object (\code{siamcat.holdout}).
#'
#' @examples
#'     data(siamcat_example)
#'     # Simple example
#'     siamcat.pred <- make.predictions(siamcat_example)
#'
#'     # Predictions on a holdout-set
#'     \dontrun{pred.mat <- make.predictions(siamcat.trained, siamcat.holdout,
#'     normalize.holdout=TRUE)}
#'
make.predictions <- function(siamcat,
    siamcat.holdout = NULL,
    normalize.holdout = TRUE,
    verbose = 1) {
    s.time <- proc.time()[3]

    # if holdout is NULL, make predictions on data in siamcat
    if (is.null(siamcat.holdout)) {
        if (verbose > 1)
            message("+ starting make.predictions on siamcat object")

        feat <- t(features(siamcat))
        label <- label(siamcat)
        data.split <- data_split(siamcat)
        models <- models(siamcat)

        label.fac <-
            factor(label$label,
                levels = c(label$negative.lab,
                    label$positive.lab))

        # assert that there is a split
        stopifnot(!is.null(data_split(siamcat)))

        num.folds <- data.split$num.folds
        num.resample <- data.split$num.resample

        pred <-
            matrix(
                NA,
                ncol = num.resample,
                nrow = length(label.fac),
                dimnames = list(names(label.fac), paste0("CV_rep",
                    seq_len(
                        num.resample
                    )))
            )
        i = 1
        if (verbose == 1 || verbose == 2)
            pb <-
            txtProgressBar(max = num.folds * num.resample, style = 3)
        for (f in seq_len(num.folds)) {
            for (r in seq_len(num.resample)) {
                test.label <- label.fac[data.split$test.folds[[r]][[f]]]
                data <-
                    as.data.frame(feat[data.split$test.folds[[r]][[f]], ])

                # assert stuff
                stopifnot(nrow(data) == length(test.label))
                stopifnot(all(rownames(data) == names(test.label)))

                model <- models[[i]]

                stopifnot(!any(
                    rownames(model$task$env$data) %in%
                        rownames(data)
                ))

                # subselect features for each model
                # needs to be added due to feature selection
                data <- data[,model$features]

                data$label <- test.label

                if (verbose > 2)
                    message(
                        paste0(
                            "Applying ",
                            model_type(siamcat),
                            " on cv_fold",
                            f,
                            "_rep",
                            r,
                            " (",
                            i,
                            " of ",
                            num.resample * num.folds,
                            ")..."
                        )
                    )

                task <-
                    makeClassifTask(data = data, target = "label")
                pdata <- predict(model, task = task)

                # rescale posterior probabilities between -1 and 1 (this works
                # only for binary data!!!!)
                # TODO: Will need adjustment and generalization in the future)
                p <- label$negative.lab + abs(label$positive.lab -
                        label$negative.lab) * pdata$data[, 4]
                names(p) <- rownames(pdata$data)
                pred[names(p), r] <- p
                i <- i + 1
                if (verbose == 1 || verbose == 2)
                    setTxtProgressBar(pb, i)
            }
        }
        stopifnot(!any(is.na(pred)))
        pred <- pred_matrix(pred)
        pred_matrix(siamcat) <- pred
        return.object <- siamcat
    } else {
        if (verbose > 1)
            message("+ starting make.predictions on external dataset")

        if (normalize.holdout) {
            if (verbose > 1)
                message("+ Performing frozen normalization on holdout set")
            siamcat.holdout <- normalize.features(siamcat.holdout,
                norm.param = norm_param(siamcat),
                verbose = verbose)
        } else {
            message("WARNING: holdout set is not being normalized!")
        }
        feat.test <- t(features(siamcat.holdout))
        feat.ref <- t(features(siamcat))
        label <- label(siamcat.holdout)
        data.split <- data_split(siamcat)
        models <- models(siamcat)
        # data sanity checks
        stopifnot(all(colnames(feat.ref) %in% colnames(feat.test)))

        # prediction
        num.models <- data.split$num.folds * data.split$num.resample

        pred <-
            matrix(
                NA,
                ncol = num.models,
                nrow = nrow(feat.test),
                dimnames = list(rownames(feat.test), paste0("Model_",
                    seq_len(num.models)))
            )
        if (verbose == 1 || verbose == 2)
            pb <-
            txtProgressBar(max = data.split$num.folds * data.split$num.resample,
            style = 3)
        for (i in seq_len(num.models)) {
            data <- as.data.frame(feat.test)
            model <- models[[i]]

            data <- data[, model$features]
            data$label <- as.factor(label$label)

            if (verbose > 2)
                message(
                    paste0(
                        "Applying ",
                        model_type(siamcat),
                        " on complete external dataset",
                        " (",
                        i,
                        " of ",
                        num.models,
                        ")..."
                    )
                )

            task <- makeClassifTask(data = data, target = "label")
            pdata <- predict(model, task = task)

            p <- label$negative.lab + abs(label$positive.lab -
                    label$negative.lab) * pdata$data[, 4]
            names(p) <- rownames(pdata$data)
            pred[names(p), i] <- p
            pred <- pred_matrix(pred)
            if (verbose == 1 || verbose == 2)
                setTxtProgressBar(pb, i)
        }
        return.object <- pred
    }

    # print correlation matrix
    if (verbose > 1)
        message(paste("Total number of predictions made:", length(pred)))
    correlation <- cor(pred, method = "spearman")
    if (verbose > 1)
        message("Correlation between predictions from repeated CV:")
    if (verbose > 1)
        message(paste(
            "Min: ",
            min(correlation),
            ", Median: ",
            median(correlation),
            ", Mean: ",
            mean(correlation)
        ))

    # print out time
    e.time <- proc.time()[3]
    if (verbose > 1)
        message(paste(
            "+ finished make.predictions in",
            formatC(e.time - s.time, digits = 3),
            "s"
        ))
    if (verbose == 1)
        message("Made predictions successfully.")

    return(return.object)
}
