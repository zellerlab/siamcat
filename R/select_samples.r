#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Select samples based on metadata
#' 
#' @description This functions selects labels and metadata based on
#'         a specific column in the metadata. Provided with a column-name in the
#'         metadata and a range or a set of allowed values, the function will
#'         filter the \link{siamcat-class} object accordingly.
#'         
#' @usage select.samples(siamcat, filter, allowed.set = NULL, 
#' allowed.range = NULL, verbose = 1) 
#'         
#' @param siamcat an object of class \link{siamcat-class}
#' 
#' @param filter string, name of the meta variable on which the selection
#'         should be done
#'         
#' @param allowed.set a vector of allowed values
#' 
#' @param allowed.range a range of allowed values
#' 
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'         for only information about progress and success, \code{2} for normal
#'         level of information and \code{3} for full debug information,
#'         defaults to \code{1}
#'         
#' @keywords SIAMCAT select.samples
#' 
#' @export
#' 
#' @return an object of class \link{siamcat-class} with labels and metadata
#'         filtered in order to contain only allowed values
#'         
#' @examples
#'     data(siamcat_example)
#'     # Select all samples that fall into an Age-range between 20 and 80 years
#'     siamcat_selected <- select.samples(siamcat_example, 'age',
#'     allowed.range=c(20, 80))
#'
#'     # Select all samples for which information about the gender is given
#'     # Provide additional information with verbose
#'     \dontrun{siamcat_selected <- select.samples(siamcat_example, 'gender',
#'     allowed.set=c(1, 2), verbose=2)}
#'     
select.samples <- function(siamcat,
    filter,
    allowed.set = NULL,
    allowed.range = NULL,
    verbose = 1) {
    if (verbose > 1)
        message("+ starting select.samples")
    s.time <- proc.time()[3]
    
    if (verbose > 2)
        message("+++ checking allowed values")
    if (!is.null(allowed.range)) {
        allowed.range <- gsub("\\[|\\]", "", allowed.range)
        allowed.range <-
            as.numeric(unlist(strsplit(allowed.range, ",")))
        stopifnot(length(allowed.range) == 2)
        stopifnot(allowed.range[1] <= allowed.range[2])
    }
    
    if (!is.null(allowed.set)) {
        # parse set
        allowed.set <- gsub("\\{|\\}", "", allowed.set)
        allowed.set <-
            as.numeric(unlist(strsplit(allowed.set, ",")))
        allowed.set <- sort(unique(allowed.set))
    }
    
    if (!xor(is.null(allowed.range), is.null(allowed.set))) {
        stop(
            "Neither allowed.range nor allowed.set (or both at the same time)
            have been provided, exiting!\n"
        )
    } else {
        if (!is.null(allowed.range)) {
            if (verbose > 2)
                message(paste0(
                    "+++ allowed.range  = [",
                    paste(allowed.range,
                        collapse = ","),
                    "]"
                ))
        } else {
            if (verbose > 2)
                message(paste0(
                    "+++ allowed.set = {",
                    paste(allowed.set,
                        collapse = ","),
                    "}"
                ))
        }
    }
    
    if (!filter %in% colnames(meta(siamcat)))
        stop("! The filter name is not present in colnames of the sample data.
            Stopping.\n")
    filter.var <- meta(siamcat)[, filter]
    
    if (!is.null(allowed.range)) {
        s.idx <- !is.na(filter.var) & filter.var >= allowed.range[1] &
            filter.var <= allowed.range[2]
        if (verbose > 1) {
            message(
                paste0(
                    "+++ removed ",
                    sum(!s.idx),
                    " samples with ",
                    filter,
                    " not in [",
                    paste(allowed.range, collapse = ", "),
                    "] (retaining ",
                    sum(s.idx),
                    ")"
                )
            )
        }
        
    } else {
        s.idx <- !is.na(filter.var) & filter.var %in% allowed.set
        if (verbose > 1) {
            message(
                paste0(
                    "+++ removed ",
                    sum(!s.idx),
                    " samples with ",
                    filter,
                    " not in {",
                    paste(allowed.set, collapse = ", "),
                    "} (retaining ",
                    sum(s.idx),
                    ")"
                )
            )
        }
    }
    
    s.names <- rownames(meta(siamcat))[s.idx]
    physeq(siamcat) <-
        prune_samples(x = physeq(siamcat), samples = s.names)
    siamcat <- filter.label(siamcat, verbose = verbose)
    e.time <- proc.time()[3]
    if (verbose > 1)
        message(paste(
            "+ finished select.samples in",
            formatC(e.time - s.time,
                digits = 3),
            "s"
        ))
    if (verbose == 1)
        message("Selecting samples finished")
    return(siamcat)
}
