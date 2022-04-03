#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Summarize features
#'
#' @description This function summarize features on a specific taxonomic level
#'
#' @usage summarize.features(siamcat, level = 'g__', 
#' feature.type='original', verbose=1)
#'
#' @param siamcat object of class \link{siamcat-class}
#'
#' @param level string, at which level to summarize (e.g. \code{g__} = genus)
#'
#' @param feature.type string, on which type of features should the function
#' work? Can be either \code{"original"}, \code{"filtered"}, or
#' \code{"normalized"}. Please only change this paramter if you know what
#' you are doing!
#'
#' @param verbose integer, control output: \code{0} for no output at all,
#' \code{1} for only information about progress and success, \code{2} for
#' normal level of information and \code{3} for full debug information,
#' defaults to \code{1}
#'
#' @keywords internal
#'
#' @export summarize.features
#'
#' @encoding UTF-8
#'
#' @details This function will summarize features at different taxonomic
#' levels, e.g. transform species-level relative abundance into genus-level
#' taxonomic profiles.
#'
#' The function expects a SIAMCAT object that either contains an entry in the
#' \link[phyloseq]{tax_table} slot of its \link{phyloseq} object, \strong{OR}
#' a set of feature names which encode taxonomic information, e.g.
#'
#' \code{k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Acidimicrobiales;..}
#'
#' Then, for a given taxonomic level (e.g. \code{g__}), the function will
#' sum up all the relative abundances of features belonging to the same group
#' at that specific taxonomic level.
#'
#' \strong{Please note that this function is currently maturing and not
#' necessarily reliable!!!}
#'
#' @return object of class \link{siamcat-class} with a summarized feature table
#'
#' @examples
#' ## load the phyloseq example data
#' data("GlobalPatterns")
#' ## create an example label
#' label <- create.label(meta=sample_data(GlobalPatterns),
#'     label = "SampleType",
#'     case = c("Freshwater", "Freshwater (creek)", "Ocean"))
#' # run the constructor function
#' siamcat <- siamcat(phyloseq=GlobalPatterns, label=label, verbose=1)
#' siamcat <- summarize.features(siamcat, level='Genus', verbose=3)
summarize.features <- function(siamcat, level = "g__",
    feature.type='original', verbose=1){

    if (verbose > 1) message("+ starting summarize.features")

    s.time <- proc.time()[3]

    if (verbose > 2) message("+++ summarizing on level: ",level)

    if (!feature.type %in% c('original', 'filtered', 'normalized')){
        stop("Unrecognised feature type, exiting...\n")
    }
    # get the right features
    if (feature.type == 'original'){
        feat <- get.orig_feat.matrix(siamcat)
    } else if (feature.type == 'filtered'){
        if (is.null(filt_feat(siamcat, verbose=0))){
            stop('Features have not yet been filtered, exiting...\n')
        }
        feat <- get.filt_feat.matrix(siamcat)
    } else if (feature.type == 'normalized'){
        if (is.null(norm_feat(siamcat, verbose=0))){
            stop('Features have not yet been normalized, exiting...\n')
        }
        feat <- get.norm_feat.matrix(siamcat)
    }

    if (is.null(tax_table(physeq(siamcat), errorIfNULL = FALSE))){
        warning(paste0("Tax_table slot is empty! Will try to infer the tax",
            " table from the feature names"))
        # make sure that seperating characters are dots
        taxa.names <- make.names(rownames(feat))
        # TODO
        # checks and balances
        # check that all feature names are on the same taxonomic Level

        # check that desired level is in the names
        tax.table <- str_split(taxa.names,
            pattern='(\\.[a-z]__)|(^[a-z]__)',
            simplify = TRUE)[,-1]
        colnames(tax.table) <- str_extract_all(taxa.names[1],
            '[a-z]__')[[1]]
        rownames(tax.table) <- rownames(feat)
        tax_table(physeq(siamcat)) <- tax.table
    } else {
        tax.table <- tax_table(physeq(siamcat))
    }

    tax.table <- tax.table[rownames(feat),]
    tax.table <- tax.table@.Data
    # search for bins
    idx <- match(level, colnames(tax.table))
    if (is.na(idx)){
        stop('Level ', level, ' not found in the feature names. Exiting...')
    }
    # rename unclassified features
    tax.table[tax.table[,idx] == '', idx] <- 'unclassified'
    tax.table[is.na(tax.table[,idx]), idx] <- 'unclassified'
    bins.unique <- unique(tax.table[,idx])

    # make empty matrix
    summarized.feat <- matrix(NA, nrow=length(bins.unique), ncol=ncol(feat),
        dimnames=list(bins.unique, colnames(feat)))

    if (verbose > 1) pb <- progress_bar$new(total=length(bins.unique))
    for (x in bins.unique){
        summarized.feat[x,] <- colSums(feat[which(tax.table[,idx] == x),,
            drop=FALSE])
        if (verbose > 1) pb$tick()
    }

    if (verbose > 2) message("+++ summarized features table contains: ",
        nrow(summarized.feat)," features")

    if (feature.type == 'original'){
        # adjust tax table and entire phyloseq object for reduced feat table
        tax.table.red <- tax.table[,seq_len(idx), drop=FALSE]
        tax.table.red <- tax.table.red[tax.table.red[,idx]!='unclassified',]
        tax.table.red <- tax.table.red[!duplicated(tax.table.red),]
        if (nrow(tax.table.red) > (nrow(summarized.feat)-1)){
            warning(paste0("Tax table does not seem to be consistent in ",
                "all cases...\nWill be collapsed at level ", level))
            for (g in unique(tax.table.red[,idx])){
                temp <- tax.table.red[tax.table.red[,idx]==g,,drop=FALSE]
                if (nrow(temp)!=1){
                    tax.table.red <- tax.table.red[-which(
                        tax.table.red[,idx]==g),]
                    temp <- vapply(colnames(temp), FUN=function(x){
                        y=temp[,x]
                        if(length(unique(y))==1){
                            return(unique(y))
                        } else {
                            return(paste(y, collapse='|'))
                        }
                    }, FUN.VALUE=character(1))
                    tax.table.red <- rbind(tax.table.red, temp)
                }
            }
        }
        # add unclassified again
        tax.table.red <- rbind(tax.table.red, NA)
        tax.table.red[nrow(tax.table.red), idx] <- 'unclassified'
        rownames(tax.table.red) <- tax.table.red[,idx]
        tax.table.red <- tax.table.red[rownames(summarized.feat),]
        temp.phyloseq <- phyloseq(tax_table=tax_table(tax.table.red),
            otu_table=otu_table(summarized.feat, taxa_are_rows=TRUE))
        if (!is.null(sample_data(physeq(siamcat), errorIfNULL=FALSE))){
            sample_data(temp.phyloseq) <- sample_data(physeq(siamcat))
        }
        if (!is.null(phy_tree(physeq(siamcat), errorIfNULL=FALSE))){
            warning("Phylogenetic tree in original SIAMCAT object had ",
                "to be deleted...")
        }
        physeq(siamcat) <- temp.phyloseq
    } else if (feature.type == 'filtered'){
        filt_feat(siamcat) <- list(
            filt.feat=otu_table(summarized.feat,taxa_are_rows = TRUE),
            filt.param=filt_params(siamcat))
    } else if (feature.type == 'normalized'){
        norm_feat(siamcat) <- list(
            norm.feat=otu_table(summarized.feat,taxa_are_rows = TRUE),
            norm.param=norm_params(siamcat))
    }

    e.time <- proc.time()[3]
    if (verbose > 1)
    message(paste("+ finished summarize.features in",
        formatC(e.time - s.time, digits = 3), "s" ))

    if (verbose == 1)
        message("Summarized features successfully.")

    return(siamcat)
}
