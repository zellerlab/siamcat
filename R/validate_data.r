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
#' @description This function checks if labels are available for all samples in features. Additionally validates metadata, if given
#' @param feat features object
#' @param label labels object
#' @param meta metadata object
#' @keywords SIAMCAT validate.data
#' @export
#' @return a list containing values after validation: \itemize{
#'  \item \code{$feat} = features;
#'  \item \code{$label} = labels;
#'  \item \code{$meta} = metadata
#' }
validate.data <- function(feat, label, meta = NULL){
  # TODO attempt multiple reading attempts!

  # Check if labels are available for all samples in features
  if (length(label$label) == dim(feat)[2]){
    stopifnot(all(names(label$label) %in% colnames(feat)) && all(colnames(feat) %in% names(label$label)))
    # if of the same length, everything should match and be in the same order
    m <- match(names(label$label), colnames(feat))
    feat <- feat[,m]
    stopifnot(all(names(label$label) == colnames(feat)))

  } else if (length(label$label) >= dim(feat)[2]) {
    # if there are more labels than samples in features, remove them in labels
    stopifnot(all(colnames(feat) %in% names(label$label)))

    m <- match(colnames(feat), names(label$label))
    # create new labels object with reduced labels
    stopifnot(label$info$type == 'BINARY')
    labels_new      <- list("label" = label$label[m],
                            "header" = label$header,
                            "info" = label$info,
                            "positive.lab" = label$positive.lab,
                            "negative.lab" = label$negative.lab,
                            "n.lab" = label$n.lab,
                            "p.lab" = label$p.lab)
    labels_new$n.idx <- labels_new$label == labels_new$negative.lab
    labels_new$p.idx <- labels_new$label == labels_new$positive.lab

    cat('Warning: Removing labels of', length(label$label) - dim(feat)[2],'sample(s) for which no features are available.\n')

    label <- labels_new

  } else if (length(label$label) <= dim(feat)[2]) {
    # if there are more samples in features, remove them and keep only the ones for which labels are available
    stopifnot(all(names(label$label) %in% colnames(feat)))

    cat('Warning: Removing', dim(feat)[2]-length(label$label),'sample(s) for which no labels are available.\n')

    m <- match(names(label$label), colnames(feat))
    feat <- feat[,m]

  }
  # Check for sample number in the different classes
  for (i in label$info$class.descr){
    if(sum(label$label==i) <= 5) stop("Data set has only",sum(label$label==i), "training examples of class",i," This is not enough for SIAMCAT to proceed")
    if (sum(label$label==i) < 10){
      cat("Data set has only",sum(label$label==i), "training examples of class",i," . Note that a dataset this small/skewed is not necessarily suitable for analysis in this pipe line." )
    }
  }

  # if metadata is available, check for overlap in labels
  if (!is.null(meta)) {
    if (length(label$label) == dim(meta)[1]){
      stopifnot(all(names(label$label) %in% rownames(meta)) && all(rownames(meta) %in% names(label$label)))
      m    <- match(names(label$label), rownames(meta))
      meta <- meta[m,]
      stopifnot(all(names(label$label) == rownames(meta)))
    } else if (length(label$label) <= dim(meta)[1]) {
      stopifnot(all(names(label$label) %in% rownames(meta)))
      m <- match(names(label$label), rownames(meta))
      cat('Warning: Removing metadata information for', dim(meta)[1]-length(label$label),'superfluous samples.\n')
      meta <- meta[m,]
    } else {
      stop('Metadata is not available for all samples!')
    }
  return(list("meta"=meta,
              "label"=label,
              "feat"=feat))
  }
  return(list("label"=label,
              "feat"=feat))
}
  ### run data checks / validation / format conversion
  # TODO !!!
