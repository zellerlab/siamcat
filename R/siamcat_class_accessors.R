#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

################################################################################
#' Universal slot accessor function for siamcat-class.
#'
#' This function is used internally by many accessors.
#'
#' @usage accessSlot(siamcat, slot, errorIfNULL=FALSE)
#'
#' @param siamcat an object of \link{siamcat-class}.
#'
#' @param slot A character string indicating the slot (not data class)
#'  of the component data type that is desired.
#'
#' @return Returns the component object specified by the argument \code{slot}.
#'  Returns NULL if slot does not exist.
#'
#' @export
#' @examples #
#' data(siamcat_example)
#' access(siamcat_example, "label")
#' access(siamcat_example, "model_list")
accessSlot <- function(siamcat, slot){
  if(!slot %in% slotNames(siamcat) ){
    # If slot is invalid, set out to NULL. Test later if this is an error.
    out = NULL
  } else {
    # By elimination, must be valid. Access slot
    out = eval(parse(text=paste("siamcat@", slot, sep="")))
  }
  if(is.null(out) ){
    # Only error regarding a NULL return value if errorIfNULL is TRUE.
    message(slot, " slot is empty.")
  }
  return(out)
}

################################################################################
#' Retrieve \link{model_list-class} from object.
#'
#'
#' @usage model_list(siamcat, errorIfNULL=TRUE)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'  that contains a model_list or instance of \link{model_list-class}.
#'
#'
#' @return The \link{model_list-class} object or NULL.
#'
#' @export
#' @rdname model_list-methods
#' @docType methods
#'
#' @examples
#'  data(siamcat_example)
#'  model_list(siamcat_example)
setGeneric("model_list", function(siamcat) standardGeneric("model_list"))
#' @rdname model_list-methods
#' @aliases model_list,ANY-method
setMethod("model_list", "ANY", function(siamcat){
  accessSlot(siamcat, "model_list")
})
# Return as-is if already a model_list object
#' @rdname model_list-methods
setMethod("model_list", "model_list", function(siamcat){ return(siamcat) })

################################################################################
#' Retrieve list of models from object.
#'
#'
#' @usage models(siamcat, errorIfNULL=TRUE)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'  that contains a model_list or instance of \link{model_list-class}.
#'
#'
#' @return The list of models or NULL.
#'
#' @export
#' @rdname models-methods
#' @docType methods
#'
#' @examples
#'  data(siamcat_example)
#'  models(siamcat_example)
setGeneric("models", function(siamcat) standardGeneric("models"))
#' @rdname models-methods
#' @aliases models,ANY-method
setMethod("models", "ANY", function(siamcat){
  eval(model_list(siamcat)@models)
})
# Return list of models if a model_list object
#' @rdname model_list-methods
setMethod("models", "model_list", function(siamcat){ return(siamcat@models) })

################################################################################
#' Retrieve model_type from object.
#'
#'
#' @usage model_type(siamcat)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'  that contains a model_list or instance of \link{model_list-class}.
#'
#'
#' @return The string describing type of model used or NULL.
#'
#' @export
#' @rdname model_type-methods
#' @docType methods
#'
#' @examples
#'  data(siamcat_example)
#'  model_type(siamcat_example)
setGeneric("model_type", function(siamcat) standardGeneric("model_type"))
#' @rdname model_type-methods
#' @aliases model_type,ANY-method
setMethod("model_type", "ANY", function(siamcat){
  eval(model_list(siamcat)@model_type)
})
# Return model type if a model_list object
#' @rdname model_list-methods
setMethod("model_type", "model_list", function(siamcat){ return(siamcat@model_type) })

################################################################################
#' Retrieve eval_data from object.
#'
#'
#' @usage eval_data(siamcat, errorIfNULL=TRUE)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'  that contains a eval_data..
#'
#' @return The eval_data list or NULL.
#'
#' @export
#' @rdname eval_data-methods
#' @docType methods
#'
#' @examples
#'  data(siamcat_example)
#'  eval_data(siamcat_example)
setGeneric("eval_data", function(siamcat) standardGeneric("eval_data"))
#' @rdname eval_data-methods
#' @aliases eval_data,ANY-method
setMethod("eval_data", "ANY", function(siamcat){
  accessSlot(siamcat, "eval_data")
})

################################################################################
#' Retrieve pred_matrix from object.
#'
#'
#' @usage pred_matrix(siamcat, errorIfNULL=TRUE)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'  that contains a pred_matrix
#'
#' @return The pred_matrix matrix or NULL.
#'
#' @export
#' @rdname pred_matrix-methods
#' @docType methods
#'
#' @examples
#'  data(siamcat_example)
#'  pred_matrix(siamcat_example)
setGeneric("pred_matrix", function(siamcat) standardGeneric("pred_matrix"))
#' @rdname pred_matrix-methods
#' @aliases pred_matrix,ANY-method
setMethod("pred_matrix", "ANY", function(siamcat){
  accessSlot(siamcat, "pred_matrix")
})

################################################################################
#' Retrieve norm_param from object.
#'
#'
#' @usage norm_param(siamcat, errorIfNULL=TRUE)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'  that contains a norm_param
#'
#'
#' @return The norm_param list or NULL.
#'
#' @export
#' @rdname norm_param-methods
#' @docType methods
#'
#' @examples
#'  data(siamcat_example)
#'  norm_param(siamcat_example)
setGeneric("norm_param", function(siamcat) standardGeneric("norm_param"))
#' @rdname norm_param-methods
#' @aliases norm_param,ANY-method
setMethod("norm_param", "ANY", function(siamcat){
  accessSlot(siamcat, "norm_param")
})

################################################################################
#' Retrieve a \link{label-class} object from object.
#'
#'
#' @usage label(siamcat, errorIfNULL=TRUE)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'  that contains a label or instance of \link{label-class}.
#'
#'
#' @return The \link{label-class} object or NULL.
#'
#' @export
#' @rdname label-methods
#' @docType methods
#'
#' @examples
#'  data(siamcat_example)
#'  label(siamcat_example)
setGeneric("label", function(siamcat) standardGeneric("label"))
#' @rdname label-methods
#' @aliases label,ANY-method
setMethod("label", "ANY", function(siamcat){
  accessSlot(siamcat, "label")
})
# Return as-is if already a label object
#' @rdname label-methods
setMethod("label", "label", function(siamcat){ return(siamcat) })

################################################################################
#' Retrieve a \link{data_split-class} object from object.
#'
#'
#' @usage data_split(siamcat, errorIfNULL=TRUE)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'  that contains a label or instance of \link{data_split-class}.
#'
#'
#' @return The \link{data_split-class} object or NULL.
#'
#' @export
#' @rdname data_split-methods
#' @docType methods
#'
#' @examples
#'  data(siamcat_example)
#'  data_split(siamcat_example)
setGeneric("data_split", function(siamcat) standardGeneric("data_split"))
#' @rdname data_split-methods
#' @aliases data_split,ANY-method
setMethod("data_split", "ANY", function(siamcat){
  accessSlot(siamcat, "data_split")
})
# Return as-is if already a data_split object
#' @rdname data_split-methods
setMethod("data_split", "data_split", function(siamcat){ return(siamcat) })

################################################################################
#' Retrieve a \link[phyloseq]{otu_table-class} object from orig_feat slot.
#'
#'
#' @usage orig_feat(siamcat, errorIfNULL=TRUE)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'  that contains a label or instance of \link[phyloseq]{otu_table-class}.
#'
#'
#' @return The \link[phyloseq]{otu_table-class} object or NULL.
#'
#' @export
#' @rdname orig_feat-methods
#' @docType methods
#'
#' @examples
#'  data(siamcat_example)
#'  data_split(siamcat_example)
setGeneric("orig_feat", function(siamcat) standardGeneric("orig_feat"))
#' @rdname orig_feat-methods
#' @aliases orig_feat,ANY-method
setMethod("orig_feat", "ANY", function(siamcat){
  accessSlot(siamcat, "orig_feat")
})
# Return as-is if already a otu_table object
#' @rdname orig_feat-methods
setMethod("orig_feat", "otu_table", function(siamcat){ return(siamcat) })

################################################################################
#' Retrieve a \link[phyloseq]{phyloseq-class} object from object.
#'
#'
#' @usage physeq(siamcat, errorIfNULL=TRUE)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'  that contains a label or instance of \link[phyloseq]{phyloseq-class}.
#'
#'
#' @return The \link[phyloseq]{phyloseq-class} object or NULL.
#'
#' @export
#' @rdname physeq-methods
#' @docType methods
#'
#' @examples
#'  data(siamcat_example)
#'  physeq(siamcat_example)
setGeneric("physeq", function(siamcat) standardGeneric("physeq"))
#' @rdname physeq-methods
#' @aliases physeq,ANY-method
setMethod("physeq", "ANY", function(siamcat){
  accessSlot(siamcat, "phyloseq")
})
# Return as-is if already a phyloseq object
#' @rdname physeq-methods
setMethod("physeq", "phyloseq", function(siamcat){ return(siamcat) })

################################################################################
#' Retrieve a \link[phyloseq]{otu_table-class} object from object.
#'
#'
#' @usage features(siamcat, errorIfNULL=TRUE)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'  that contains a label or instance of \link[phyloseq]{otu_table-class} .
#'
#'
#' @return The \link[phyloseq]{otu_table-class} object or NULL.
#'
#' @export
#' @rdname features-methods
#' @docType methods
#'
#' @examples
#'  data(siamcat_example)
#'  features(siamcat_example)
setGeneric("features", function(siamcat) standardGeneric("features"))
#' @rdname features-methods
#' @aliases features,ANY-method
setMethod("features", "ANY", function(siamcat){
  otu_table(physeq(siamcat))
})
# Return as-is if already a otu_table object
#' @rdname features-methods
setMethod("features", "otu_table", function(siamcat){ return(siamcat) })


################################################################################
#' Retrieve a \link[phyloseq]{sample_data-class} object from object.
#'
#'
#' @usage meta(siamcat, errorIfNULL=TRUE)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'  that contains a label or instance of \link[phyloseq]{sample_data-class}.
#'
#'
#' @return The \link[phyloseq]{sample_data-class} object or NULL.
#'
#' @export
#' @rdname meta-methods
#' @docType methods
#'
#' @examples
#'  data(siamcat_example)
#'  meta(siamcat_example)
setGeneric("meta", function(siamcat) standardGeneric("meta"))
#' @rdname meta-methods
#' @aliases meta,ANY-method
setMethod("meta", "ANY", function(siamcat){
  sample_data(physeq(siamcat))
})
# Return as-is if already a sample_data object
#' @rdname meta-methods
setMethod("meta", "sample_data", function(siamcat){ return(siamcat) })

#' Access features in siamcat@phylose@otu_table
#' @title get.features
#' @name get.features
#' @description Function to access features in siamcat@phylose@otu_table
#' @param siamcat an object of class \link{siamcat-class}
#' @return Features as \link[phyloseq]{otu_table-class} object
#' @export
#' @examples
#'  data(siamcat_example)
#'  feat <- get.features(siamcat_example)
get.features <- function(siamcat) {
    return(siamcat@phyloseq@otu_table)
}

#' Access phyloseq object in siamcat@phyloseq
#' @title get.phyloseq
#' @name get.phyloseq
#' @description Function to access phyloseq object in siamcat@phylose
#' @param siamcat an object of class \link{siamcat-class}t
#' @return Object of class \link[phyloseq]{phyloseq-class}
#' @export
#' @examples
#'  data(siamcat_example)
#'  phyloseq <- get.phyloseq(siamcat_example)
get.phyloseq <- function(siamcat) {
    return(siamcat@phyloseq)
}

#' Access features in siamcat@phylose@otu_table as matrix
#' @title get.features.matrix
#' @name get.features.matrix
#' @description Function to access features in siamcat@phylose@otu_table
#' @param siamcat an object of class \link{siamcat-class}t
#' @return Features as a matrix
#' @export
#' @examples
#'  data(siamcat_example)
#'  feat <- get.features.matrix(siamcat_example)
<<<<<<< HEAD
features.as.matrix <- function(otu_table) {
    return(matrix(otu_table, nrow = nrow(otu_table), ncol = ncol(otu_table),
        dimnames = list(rownames(otu_table), colnames(otu_table))))
=======
get.features.matrix <- function(siamcat) {
    return(siamcat@phyloseq@otu_table@.Data)
>>>>>>> 806eaaef6b4b101a426d496696b9b7ee526b60e8
}

#' Access original features in siamcat@orig_feat as matrix
#' @title get.orig_feat.matrix
#' @name get.orig_feat.matrix
#' @description Function to access original features in siamcat@orig_feat
#' @param siamcat an object of class \link{siamcat-class}t
#' @return Original features as a matrix
#' @export
#' @examples
#'  data(siamcat_example)
#'  orig_feat <- get.orig_feat.matrix(siamcat_example)
get.orig_feat.matrix <- function(siamcat) {
    return(siamcat@orig_feat@.Data)
}

#' Access features in siamcat@phylose@otu_table
#' @title get.model_list
#' @name get.model_list
#' @description Function to access model_list in siamcat@get.model_list
#' @param siamcat an object of class \link{siamcat-class}t
#' @return List of models in a \link{model_list-class} object
#' @export
#' @examples
#'  data(siamcat_example)
#'  model_list <- get.model_list(siamcat_example)
get.model_list <- function(siamcat) {
    return(siamcat@model_list)
}
#' Access model type in siamcat@model_list@model.type
#' @title get.model.type
#' @name get.model.type
#' @description Function to access features in siamcat@model_list@model.type
#' @param siamcat an object of class \link{siamcat-class}
#' @return Character string specifing the name of the method used to construct the model.
#' @export
#' @examples
#'  data(siamcat_example)
#'  get.model.type(siamcat_example)
get.model.type <- function(siamcat) {
    return(model.type = siamcat@model_list@model.type)
}

#' Access models in siamcat@model_list@models
#' @title get.models
#' @name get.models
#' @description Function to access models in siamcat@model_list@models
#' @param siamcat an object of class \link{siamcat-class}
#' @return List of models
#' @export
#' @examples
#'  data(siamcat_example)
#'  models <- get.models(siamcat_example)
get.models <- function(siamcat) {
    return(siamcat@model_list@models)
}
#' Access evaluation data in siamcat@eval_data
#' @title get.eval_data
#' @name get.eval_data
#' @description Function to access evaluation data in siamcat@eval_data
#' @param siamcat an object of class \link{siamcat-class}t
#' @return Evaluation data
#' @export
#' @examples
#'  data(siamcat_example)
#'  eval_data <- get.eval_data(siamcat_example)
get.eval_data <- function(siamcat) {
    return(siamcat@eval_data)
}
#' Access Prediction matrix  in siamcat@pred_matrix
#' @title get.pred_matrix
#' @name get.pred_matrix
#' @description Function to access prediction matrix  in siamcat@pred_matrix
#' siamcat@orig_feat in an object of class \link{siamcat-class}
#' @param siamcat an object of class \link{siamcat-class}t
#' @return Prediction matrix
#' @export
#' @examples
#'  data(siamcat_example)
#'  pred_matrix <- get.pred_matrix(siamcat_example)
get.pred_matrix <- function(siamcat) {
    return(siamcat@pred_matrix)
}

#' Access labels in siamcat@label@label
#' @title get.label.label
#' @name get.label.label
#' @description Function to access labels in siamcat@label@label
#' @param siamcat an object of class \link{siamcat-class}
#' @return Labels as a vector
#' @export
#' @examples
#'  data(siamcat_example)
#'  label <- get.label.label(siamcat_example)
get.label.label <- function(siamcat) {
    return(siamcat@label@label)
}

#' Access label object in siamcat@label
#' @title get.label
#' @name get.label
#' @description Function to access labels in siamcat@label
#' @param siamcat an object of class \link{siamcat-class}
#' @return an object of class \link{label-class}
#' @export
#' @examples
#'  data(siamcat_example)
#'  label <- get.label(siamcat_example)
get.label <- function(siamcat) {
    return(siamcat@label)
}

#' Access labels in siamcat@label
#' @title get.label.list
#' @name get.label.list
#' @description Function to access labels in siamcat@label
#' @param siamcat an object of class \link{siamcat-class}
#' @return Label object converted to a list
#' @export
#' @examples
#'  data(siamcat_example)
#'  label <- get.label.list(siamcat_example)
get.label.list <- function(siamcat) {
    return(list(label = siamcat@label@label, header = siamcat@label@header, info = siamcat@label@info, positive.lab = siamcat@label@positive.lab,
        negative.lab = siamcat@label@negative.lab, n.idx = siamcat@label@n.idx, p.idx = siamcat@label@p.idx, n.lab = siamcat@label@n.lab,
        p.lab = siamcat@label@p.lab))
}

#' Access label info in siamcat@label@info
#' @title get.label.info
#' @name get.label.info
#' @description Function to access label info in siamcat@label@info
#' @param siamcat an object of class \link{siamcat-class}
#' @return List of label informations
#' @export
#' @examples
#'  data(siamcat_example)
#'  label_info <- get.label.info(siamcat_example)
get.label.info <- function(siamcat) {
    return(siamcat@label@info)
}

#' Access the data split in siamcat@data_split
#' @title get.data.split
#' @name get.data.split
#' @description Function to access the data split info in siamcat@data_split
#' @param siamcat an object of class \link{siamcat-class}
#' @return List of the data split
#' @export
#' @examples
#'  data(siamcat_example)
#'  data_split <- get.data.split(siamcat_example)
get.data.split <- function(siamcat) {
    return(list(num.folds=siamcat@data_split@num.folds, num.resample=siamcat@data_split@num.resample, test.folds=siamcat@data_split@test.folds, training.folds=siamcat@data_split@training.folds))
}
