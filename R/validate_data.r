###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# RScript flavor
#
# written by Georg Zeller
# with additions by Nicolai Karcher and Konrad Zych
# EMBL Heidelberg 2012-2017
#
# version 0.2.0
# file last updated: 26.06.2017
# GNU GPL 3.0
###

#' @title Validate samples in metadata, labels, and features
#' @description This function checks if metadata is available for all samples in labels and vice versa. If yes, order metadata according to the order found in labels.
#' @param feat features object
#' @param label labels object
#' @param meta metadata object
#' @keywords SIAMCAT validate.data
#' @export
#' @return a list containing values after filtering: \itemize{
#'  \item \code{$feat} = features;
#'  \item \code{$label} = labels;
#'  \item \code{$meta} = metadata
#' }
validate.data <- function(feat, label, meta = NULL){
  # TODO attempt multiple reading attempts!
  # TODO 2: should this function also validate the features? (at the moment, nothing happens to the features...)
  # TODO 3: this functions takes label$label at the moment, not an label object. also does not return an label obejct. Should be changed, i guess?

  if (!is.null(meta)) {
  stopifnot(all(names(label) %in% rownames(meta)) && all(rownames(meta) %in% names(label)))
  m    <- match(names(label), rownames(meta))
  meta <- meta[m,]
  stopifnot(all(names(label) == rownames(meta)))
  return(list("meta"=meta,
              "label"=label,
              "feat"=feat))
  }
  return(list("label"=label,
              "feat"=feat))
}
  ### run data checks / validation / format conversion
  # TODO !!!
