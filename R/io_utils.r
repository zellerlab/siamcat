#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Read feature file
#'
#' @description This file reads in the tsv file with features and
#' converts it into a matrix.
#'
#' The file should be oragnized as follows:
#' features (in rows) x samples (in columns).
#'
#' First row should contain sample labels (consistent with label data), while
#' the first column should contain feature labels (e.g. taxonomic identifiers).
#' The remaining entries are expected to be real values \code{>= 0} that
#' quantify the abundance of each feature in each sample.
#'
#' @param fn.in.feat name of the tsv file containing features
#'
#' @param transpose should the features table be transposed?
#'
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'     for information about progress and time, defaults to \code{0}
#'
#' @export
#'
#' @return \code{otu_table} containing features from the file
#'
#' @examples
#'     # run with example data
#'     fn.feat <- system.file('extdata', 'feat_crc_zeller_msb_mocat_specI.tsv',
#'     package = 'SIAMCAT')
#'
#'     features <- read.features(fn.feat)
read.features <- function(fn.in.feat, transpose = FALSE, verbose = 0) {
    if (verbose > 1)
        message("+ starting read.features")
    s.time <- proc.time()[3]
    if (is.null(fn.in.feat))
        stop("Filename for features file not provided!\n")

    feat <- read.table(
        file = fn.in.feat,
        sep = "\t",
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        quote = "",
        encoding = "latin1"
    )
    feat <- as.matrix(feat)
    if(transpose) feat <- t(feat)

    featNames <- make.names(rownames(feat))

    if (any(rownames(feat) != featNames)) {
        message(
            "The provided feature names were not semantically correct for
            use in R, they were updated."
        )
        rownames(feat) <- featNames
    }
    e.time <- proc.time()[3]
    if (verbose > 0)
        message(paste(
            "+ finished read.features in",
            formatC(e.time - s.time, digits = 3),
            "s"
        ))
    invisible(otu_table(feat, taxa_are_rows = TRUE))
}

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
#' @param fn.in.label name of the tsv file containing labels
#'
#' @export
#'
#' @return label object containing several entries:\itemize{
#' \item \code{$label} named vector containing the numerical labels from the
#' file;
#' \item \code{$header} first row of the label file;
#' \item \code{$info} information about the type of label (e.g. \code{BINARY});
#' \item \code{$positive.lab} numerical label for controls, e.g. \code{-1};
#' \item \code{$negative.lab} numerical label for cases, e.g. \code{1};
#' \item \code{$n.idx} logical vector of labels (\code{TRUE} for controls,
#' \code{FALSE} otherwise);
#' \item \code{$n.lab} label for controls, e.g. \code{healthy};
#' \item \code{$p.idx} logical vector of labels (\code{TRUE} for cases,
#' \code{FALSE} otherwise);
#' \item \code{$p.lab} label for cases, e.g. \code{cancer}
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
    label <- list(label = label, header = header)
    label$info <- parse.label.header(label$header)
    stopifnot(label$info$type == "BINARY")
    label$positive.lab <- max(label$info$class.descr)
    label$negative.lab <- min(label$info$class.descr)

    label$n.idx <- label$label == label$negative.lab
    label$n.lab <- gsub("[_.-]", " ", names(label$info$class.descr)
        [label$info$class.descr == label$negative.lab])

    label$p.idx <- label$label == label$positive.lab
    label$p.lab <- gsub("[_.-]", " ", names(label$info$class.descr)
        [label$info$class.descr == label$positive.lab])

    labelRes <-
        label(
            list(
                label = label$label,
                header = label$header,
                info = label$info,
                positive.lab = label$positive.lab,
                negative.lab = label$negative.lab,
                n.idx = label$n.idx,
                p.idx = label$p.idx,
                n.lab = label$n.lab,
                p.lab = label$p.lab
            )
        )
    invisible(labelRes)
        }

#' @title Read metadata file
#'
#' @description This file reads in the tsv file with numerical metadata and
#' converts it into a matrix.
#'
#' The file should be organized as follows:
#' samples (in rows) x metadata (in columns). Metadata needs to be converted to
#' numerical values by the user.
#'
#' Metadata may be optional for the SIAMCAT workflow, but are necessary for
#' heatmap displays, see \link{model.interpretation.plot}
#'
#' @param fn.in.meta name of the tsv file containing metadata
#'
#' @export
#'
#' @return \code{sample_data} object
#'
#' @examples
#'     # run with example data
#' fn.meta  <- system.file('extdata',
#' 'num_metadata_crc_zeller_msb_mocat_specI.tsv',
#' package = 'SIAMCAT')
#'
#' meta_data <- read.meta(fn.meta)

read.meta <- function(fn.in.meta) {
    if (is.null(fn.in.meta) || toupper(fn.in.meta) == "NULL" ||
            toupper(fn.in.meta) == "NONE" || toupper(fn.in.meta) ==
            "UNKNOWN") {
        warning("Filename for metadata file not provided, continuing
            without it.\n")
    } else {
        meta <-
            read.table(
                file = fn.in.meta,
                sep = "\t",
                header = TRUE,
                row.names = 1,
                check.names = FALSE,
                quote = "",
                encoding = "latin1"
            )
    }
    invisible(sample_data(meta))
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
    label.info$class.descr <- class.descr
    return(label.info)
}