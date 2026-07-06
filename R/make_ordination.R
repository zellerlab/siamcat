
#' @title Generate ordination of a SIAMCAT object
#'
#' @description A thin wrapper around 
#' \code{\link[phyloseq]{ordinate}}.
#' Computes and stores an ordination of the
#' phyloseq object contained within a \link{siamcat-class} object
#' in the ordination(siamcat) slot.
#' 
#' @usage make.ordination(siamcat, distance="bray", method="Pcoa")
#'
#' @param siamcat object of class \link{siamcat-class}
#' 
#' @param method string, ordination method passed to
#' \code{\link[phyloseq]{ordinate}}. Supported methods include
#' \code{c("PCoA", "MDS", "NMDS", "DPCoA", "CAP", "RDA", "CCA", "DCA")},
#' defaults to \code{"PCoA"}
#'
#' @param distance string, distance metric passed to
#' \code{\link[phyloseq]{ordinate}}, defaults to \code{"bray"}
#' (Bray-Curtis dissimilarity)
#' 
#' @param verbose integer, control output: \code{0} for no output at all,
#' \code{1} for only information about progress and success, \code{2} for
#' normal level of information and \code{3} for full debug information,
#' defaults to \code{1}
#' 
#' @param feature.type string, on which type of features should the function
#' work? Can be either \code{"original"}, \code{"filtered"}.
#' Please only change this paramter if you know what
#' you are doing!

make.ordination <- function(siamcat, distance="bray", method="PCoA", feature.type="filtered", verbose=1){
    # get the right features
    if (feature.type == 'original'){
        feat <- get.orig_feat.matrix(siamcat)
        if (verbose > 1) message('+ using original features')
    } else if (feature.type == 'filtered'){
        if (is.null(filt_feat(siamcat, verbose=0))){
            stop('Features have not yet been filtered, exiting...\n')
        }
        feat <- get.filt_feat.matrix(siamcat)
    } else if (feature.type == 'normalized'){
        stop("Normalised features are not allowed for ordination.")
    }
    # drop samples with complete zeros to avoid mathematical errors
    # these can arise because of filtering
    feat <- feat[,colSums(feat) > 1e-4]
    temp_phyloseq <- phyloseq(otu_table=otu_table(feat, taxa_are_rows=TRUE))
    dmat <- phyloseq::distance(temp_phyloseq, method = distance)
    ordination(siamcat) <- list(
        ord = phyloseq::ordinate(temp_phyloseq, method = method, distance = dmat),
        distmat = dmat,
        distance = distance,
        method = method
    )
    return(siamcat)
}