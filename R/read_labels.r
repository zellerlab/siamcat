#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Read labels file
#'
#' @description This file reads in the tsv file with labels and converts it
#' into a label object.
#'
#' First row is expected to be \code{#BINARY:1=[label for cases];
#' -1=[label for controls]}.
#' Second row should contain the sample identifiers as tab-separated list
#' (consistent with feature and metadata).
#'
#' Third row is expected to contain the actual class labels (tab-separated):
#' \code{1} for each case and \code{-1} for each control.
#'
#' Note: Labels can take other numeric values (but not characters or strings);
#' importantly, the label for cases has to be greater than the one for controls
#'
#' @param fn.in.label name of the tsv file containing labels
#'
#' @export
#'
#' @return label object containing several entries:\itemize{
#' \item \code{$label} named vector containing the numerical labels from the
#' file;
#' \item \code{$info} information about the classes in the label;
#' \item \code{$type} information about the label type (e.g. \code{BINARY});
#'}
#'
#' @examples
#'     # run with example data
#' fn.label <- system.file('extdata', 'label_crc_zeller_msb_mocat_specI.tsv',
#'     package = 'SIAMCAT')
#'
#' labels <- read.labels(fn.label)
read.labels <- function(fn.in.label) {
    if (is.null(fn.in.label))
        stop("Filename for labels file not provided!\n")
    label <-
        read.table(
            file = fn.in.label,
            sep = "\t",
            header = TRUE,
            row.names = NULL,
            stringsAsFactors = FALSE,
            check.names = FALSE,
            quote = "",
            comment.char = "#",
            encoding = "latin1"
        )
    label <- as.matrix(label)
    if (dim(label)[1] > dim(label)[2]) {
        temp <- names(label)
        names(label) <- NULL
        label <- rbind(temp, label)
        rownames(label) <- label[, 1]
        label[, 1] <- NULL
        label <- t(label)
    }
    namesL <- colnames(label)
    label <- as.numeric(label)
    names(label) <- namesL

    # Check general suitablity of supplied dataset
    classes <- unique(label)
    for (i in classes) {
        if (sum(label == i) <= 5)
            stop(
                "Data set has only",
                sum(label == i),
                "training examples of class",
                i,
                " This is not enough for
                SIAMCAT to proceed"
            )
        if (sum(label == i) < 10) {
            message(
                paste(
                    "Data set has only",
                    sum(label == i),
                    "training examples of class",
                    i,
                    " . Note that a dataset this
                    small/skewed is not necessarily suitable for analysis in
                    this pipeline."
                )
                )
        }
        }

    # Check label header!
    con <- file(fn.in.label, "rt")
    header <- readLines(con, 1)
    if (substring(header, 1, 1) != "#") {
        stop("Label header seems to be missing or broken.")
    }
    close(con)
    label <- list(label = label)
    label.info <- parse.label.header(header)
    label$info <- label.info$class.descr
    label$type <- label.info$type
    stopifnot(label$type == "BINARY")
    labelRes <- label(label)
    invisible(labelRes)
    }

##### auxiliary function to trim whitespace from string returns string without
##### leading or trailing whitespace
#' @keywords internal
trim <- function(x) {
    gsub("^\\s+|\\s+$", "", x)
}

#' @title Parse label header
#'
#' @description This function parses the header of a label file
#'
#' @param  label.header - string in the format: #<TYPE>:<L1>=<class1>;
#' <L2>=<class2>[;<L3>=<class3>] where <TYPE> is a string specifying the type
#' of label variable such as BINARY (for binary classification), CATEGORICAL
#' (for multi-class classification), or CONTINUOUS (for regression)
#' <L1> is a short numeric label for the first class with description <class1>
#' (similarly for the other classes)
#'
#' @return a list with tow items \itemize{
#' \item \code{$type} type of the label: BINARY CONTINUOUS or CATEGORICAL
#' \item \code{$class.descr} lables and information on what do they mean
#'}
#'
#' @keywords internal
parse.label.header <- function(label.header) {
    s <- strsplit(label.header, ":")[[1]]
    type <- trim(s[1])
    if (substr(type, 1, 1) == "#")
        type <- trim(substr(type, 2, nchar(type)))
    class.descr <-
        unlist(strsplit(strsplit(trim(s[2]), ";")[[1]], "="))
    l <- class.descr[seq(2, length(class.descr), 2)]
    class.descr <-
        as.numeric(class.descr[seq(1, length(class.descr) - 1, 2)])
    names(class.descr) <- l

    label.info <- list()
    label.info$type <- type
    label.info$class.descr <- sort(class.descr)
    return(label.info)
}
