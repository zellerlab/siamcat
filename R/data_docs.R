#' @title SIAMCAT example
#'
#' @description Reduced version of the CRC dataset from Zeller et al. MSB 2014
#' (see \url{http://msb.embopress.org/content/10/11/766}), containing 100
#' features (15 associated features at 5\% FDR in the original dataset and
#' 85 random other features) and 141 samples, saved after the complete SIAMCAT
#' pipelinehas been run.
#'
#' Thus, the example dataset contains entries in every slot of the SIAMCAT
#' object (see \link{siamcat-class}), e.g, \code{eval_data} or
#' \code{data_split}.
#'
#' Mainly used for running the examples in the function documentation.
#'
#' @name siamcat_example
#'
#' @source \url{http://msb.embopress.org/content/10/11/766}
#'
#' @docType data
#'
#' @keywords data
NULL

#' @title Example feature matrix
#'
#' @description Feature matrix (as data.frame) of the CRC dataset from Zeller
#' et al. MSB 2014 (see \url{http://msb.embopress.org/content/10/11/766}),
#' containing 141 samples and 1754 bacterial species (features).
#'
#' @name feat.crc.zeller
#'
#' @source \url{http://msb.embopress.org/content/10/11/766}
#'
#' @docType data
#'
#' @keywords data
NULL

#' @title Example metadata matrix
#'
#' @description Metadata (as data.frame) of the CRC dataset from Zeller et
#' al. MSB 2014 (see \url{http://msb.embopress.org/content/10/11/766}),
#' containing 6 metadata variables variables (e.g. Age or BMI) for 141 samples.
#'
#' @name meta.crc.zeller
#'
#' @source \url{http://msb.embopress.org/content/10/11/766}
#'
#' @docType data
#'
#' @keywords data
NULL
