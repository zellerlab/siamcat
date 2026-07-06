#' @title Test for differences in beta diversity (PERMANOVA)
#'
#' @description Tests whether groups of samples differ in their
#' beta-diversity (distance-based) profile, using \code{\link[vegan]{adonis2}}
#' (PERMANOVA) on the distance matrix stored in the \link{siamcat-class}
#' object's ordination slot. Optionally also tests homogeneity of group
#' dispersions via \code{\link[vegan]{betadisper}}.
#'
#' @usage test.ordination(siamcat, group.by, permutations=999,
#' test.homogeneity=TRUE, verbose=1)
#'
#' @param siamcat object of class \link{siamcat-class}, must already contain
#' an ordination, see \link{make.ordination}
#'
#' @param group.by string, name of a single column in the sample metadata
#' to test group differences for (e.g. \code{"Group"}).
#'
#' @param permutations integer, number of permutations for
#' \code{\link[vegan]{adonis2}} and \code{\link[vegan]{permutest}},
#' defaults to \code{999}
#'
#' @param test.homogeneity logical, should
#' \code{\link[vegan]{betadisper}}/\code{\link[vegan]{permutest}} also be run
#' to test homogeneity of group dispersions? A significant PERMANOVA result
#' can be driven by dispersion differences rather than centroid differences,
#' so this is recommended as a companion test, defaults to \code{TRUE}
#'
#' @param verbose integer, control output: \code{0} for no output at all,
#' \code{1} for only information about progress and success, \code{2} for
#' normal level of information and \code{3} for full debug information,
#' defaults to \code{1}
#'
#' @export
#'
#' @return object of class \link{siamcat-class} with the test results stored
#' in \code{ordination(siamcat)$test}, a list with elements \code{adonis}
#' (the \code{adonis2} result) and, if \code{test.homogeneity=TRUE},
#' \code{betadisper}, itself a list with elements \code{disp} (the
#' \code{betadisper} object) and \code{permtest} (the dispersion
#' permutation test result)

test.ordination <- function(siamcat, group.by, permutations=999,
    test.homogeneity=TRUE, verbose=1){

    if (verbose > 1) message("+++ Starting test.ordination")
    if (is.null(ordination(siamcat))) {
        stop("This siamcat object does not contain an ordination. Generate one with make.ordination.")
    }
    dist_matrix <- ordination(siamcat)$distmat
    if (is.null(dist_matrix)) {
        stop("This siamcat object's ordination does not contain a distance matrix. Please regenerate it with make.ordination.")
    }

    if (verbose > 1) message("+++ Extracting metadata")
    meta <- sample_data(siamcat@phyloseq)
    # samples may differ because of dropping zero-abundance samples
    # see make.ordination
    meta <- meta[rownames(meta) %in% labels(dist_matrix)]
    meta <- as(meta, "data.frame")
    # make sure metadata is in the same order as the distance matrix
    meta <- meta[labels(dist_matrix), , drop=FALSE]
    
    if (length(group.by) != 1 || !is.character(group.by)) {
        stop("group.by must be a non-empty string specifying a column in the sample metadata.")
    } else if (!(group.by %in% colnames(meta))) {
        stop("group.by column(s) not found in sample data.")
    }

    if (verbose > 1) message("+++ Running PERMANOVA (adonis2)")
    form <- as.formula(paste("dist_matrix ~", group.by))
    adonis.res <- vegan::adonis2(form, data = meta, permutations = permutations)

    if (verbose > 0) {
        message("+++ PERMANOVA result:")
        print(adonis.res)
    }

    result <- list(adonis = adonis.res)

    if (test.homogeneity) {
        if (verbose > 1) message("+++ Testing homogeneity of dispersions (betadisper)")
        disp <- vegan::betadisper(dist_matrix, meta[[group.by]])
        permtest <- vegan::permutest(disp, permutations = permutations)
        if (verbose > 0) {
            message("+++ Dispersion homogeneity result:")
            print(permtest)
        }
        result$betadisper <- list(disp=disp, permutest=permtest)
    }

    if (verbose > 1) message("+++ Storing test results in ordination slot")
    ord <- ordination(siamcat)
    ord$test <- result
    ordination(siamcat) <- ord

    return(siamcat)
}