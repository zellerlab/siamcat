#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

################################################################################
#' Universal slot accessor function for siamcat-class.
#'
#' This function is used internally by many accessors.
#'
#' @usage accessSlot(siamcat, slot)
#'
#' @param siamcat an object of \link{siamcat-class}.
#'
#' @param slot A character string indicating the slot (not data class)
#'     of the component data type that is desired.
#'
#' @return Returns the component object specified by the argument \code{slot}.
#'     Returns NULL if slot does not exist.
#'
#' @export
#' @examples #
#' data(siamcat_example)
#' accessSlot(siamcat_example, "label")
#' accessSlot(siamcat_example, "model_list")
accessSlot <- function(siamcat, slot) {
    if (!slot %in% slotNames(siamcat)) {
        # If slot is invalid, set out to NULL. Test later if this is an error.
        out = NULL
    } else {
        # By elimination, must be valid. Access slot
        out = eval(parse(text = paste("siamcat@", slot, sep = "")))
    }
    if (is.null(out)) {
        # Only error regarding a NULL return value if errorIfNULL is TRUE.
        message(slot, " slot is empty.")
    }
    return(out)
}

###############################################################################
#' Retrieve \link{model_list-class} from object.
#'
#'
#' @usage model_list(siamcat)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a model_list or instance of \link{model_list-class}.
#'
#'
#' @return The \link{model_list-class} object or NULL.
#'
#' @export
#' @rdname model_list-methods
#' @docType methods
#'
#' @examples
#'     data(siamcat_example)
#'     model_list(siamcat_example)
setGeneric("model_list", function(siamcat)
    standardGeneric("model_list"))
#' @rdname model_list-methods
#' @aliases model_list,ANY-method
setMethod("model_list", "ANY", function(siamcat) {
    accessSlot(siamcat, "model_list")
})
# Return as-is if already a model_list object
#' @rdname model_list-methods
setMethod("model_list", "model_list", function(siamcat) {
    return(siamcat)
})

###############################################################################
#' Retrieve list of models from object.
#'
#'
#' @usage models(siamcat)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a model_list or instance of \link{model_list-class}.
#'
#'
#' @return The list of models or NULL.
#'
#' @export
#' @rdname models-methods
#' @docType methods
#'
#' @examples
#'     data(siamcat_example)
#'     models(siamcat_example)
setGeneric("models", function(siamcat)
    standardGeneric("models"))
#' @rdname models-methods
#' @aliases models,ANY-method
setMethod("models", "ANY", function(siamcat) {
    eval(model_list(siamcat)@models)
})
# Return list of models if a model_list object
#' @rdname model_list-methods
setMethod("models", "model_list", function(siamcat) {
    return(siamcat@models)
})

###############################################################################
#' Retrieve model_type from object.
#'
#'
#' @usage model_type(siamcat)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a model_list or instance of \link{model_list-class}.
#'
#'
#' @return The string describing type of model used or NULL.
#'
#' @export
#' @rdname model_type-methods
#' @docType methods
#'
#' @examples
#'     data(siamcat_example)
#'     model_type(siamcat_example)
setGeneric("model_type", function(siamcat)
    standardGeneric("model_type"))
#' @rdname model_type-methods
#' @aliases model_type,ANY-method
setMethod("model_type", "ANY", function(siamcat) {
    eval(model_list(siamcat)@model.type)
})
# Return model type if a model_list object
#' @rdname model_list-methods
setMethod("model_type", "model_list", function(siamcat) {
    return(siamcat@model_type)
})

###############################################################################
#' Retrieve eval_data from object.
#'
#'
#' @usage eval_data(siamcat)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a eval_data..
#'
#' @return The eval_data list or NULL.
#'
#' @export
#' @rdname eval_data-methods
#' @docType methods
#'
#' @examples
#'     data(siamcat_example)
#'     eval_data(siamcat_example)
setGeneric("eval_data", function(siamcat)
    standardGeneric("eval_data"))
#' @rdname eval_data-methods
#' @aliases eval_data,ANY-method
setMethod("eval_data", "ANY", function(siamcat) {
    accessSlot(siamcat, "eval_data")
})
# constructor; for creating eval_data from a list
#' @rdname eval_data-methods
#' @aliases eval_data
setMethod("eval_data", "list", function(siamcat) {
    return(new("eval_data",
        siamcat))
})

###############################################################################
#' Retrieve pred_matrix from object.
#'
#'
#' @usage pred_matrix(siamcat)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a pred_matrix
#'
#' @return The pred_matrix matrix or NULL.
#'
#' @export
#' @rdname pred_matrix-methods
#' @docType methods
#'
#' @examples
#'     data(siamcat_example)
#'     pred_matrix(siamcat_example)
setGeneric("pred_matrix", function(siamcat)
    standardGeneric("pred_matrix"))
#' @rdname pred_matrix-methods
#' @aliases pred_matrix,ANY-method
setMethod("pred_matrix", "ANY", function(siamcat) {
    accessSlot(siamcat, "pred_matrix")
})
# constructor; for creating pred_matrix from a matrix
#' @rdname pred_matrix-methods
#' @aliases pred_matrix
setMethod("pred_matrix", "matrix", function(siamcat) {
    return(new("pred_matrix",
        siamcat))
})


###############################################################################
#' Retrieve norm_param from object.
#'
#'
#' @usage norm_param(siamcat)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a norm_param
#'
#'
#' @return The norm_param list or NULL.
#'
#' @export
#' @rdname norm_param-methods
#' @docType methods
#'
#' @examples
#'     data(siamcat_example)
#'     norm_param(siamcat_example)
setGeneric("norm_param", function(siamcat)
    standardGeneric("norm_param"))
#' @rdname norm_param-methods
#' @aliases norm_param,ANY-method
setMethod("norm_param", "ANY", function(siamcat) {
    accessSlot(siamcat, "norm_param")
})
# constructor; for creating norm_param from a list
#' @rdname norm_param-methods
#' @aliases norm_param
setMethod("norm_param", "list", function(siamcat) {
    return(new("norm_param",
        siamcat))
})

################################################################################
#' Retrieve a \link{label-class} object from object.
#'
#'
#' @usage label(siamcat)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a label or instance of \link{label-class} or a list.
#'
#'
#' @return The \link{label-class} object or NULL.
#'
#' @export
#' @rdname label-methods
#' @docType methods
#'
#' @examples
#'     data(siamcat_example)
#'     label(siamcat_example)
setGeneric("label", function(siamcat)
    standardGeneric("label"))
#' @rdname label-methods
#' @aliases label,ANY-method
setMethod("label", "ANY", function(siamcat) {
    accessSlot(siamcat, "label")
})
# Return as-is if already a label object
#' @rdname label-methods
#' @aliases label
setMethod("label", "label", function(siamcat) {
    return(siamcat)
})
# constructor; for creating label from a list
#' @rdname label-methods
#' @aliases label
setMethod("label", "list", function(siamcat) {
    return(new("label", siamcat))
})

################################################################################
#' Retrieve a \link{data_split-class} object from object.
#'
#'
#' @usage data_split(siamcat)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a label or instance of \link{data_split-class} or a list.
#'
#'
#' @return The \link{data_split-class} object or NULL.
#'
#' @export
#' @rdname data_split-methods
#' @docType methods
#'
#' @examples
#'     data(siamcat_example)
#'     data_split(siamcat_example)
setGeneric("data_split", function(siamcat)
    standardGeneric("data_split"))
#' @rdname data_split-methods
#' @aliases data_split,ANY-method
setMethod("data_split", "ANY", function(siamcat) {
    accessSlot(siamcat, "data_split")
})
# Return as-is if already a data_split object
#' @rdname data_split-methods
setMethod("data_split", "data_split", function(siamcat) {
    return(siamcat)
})
# Return as-is if already a data_split object
#' @rdname data_split-methods
setMethod("data_split", "list", function(siamcat) {
    return(new("data_split", siamcat))
})

################################################################################
#' Retrieve a \link[phyloseq]{otu_table-class} object from orig_feat slot.
#'
#'
#' @usage orig_feat(siamcat)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a label or instance of \link[phyloseq]{otu_table-class}.
#'
#'
#' @return The \link[phyloseq]{otu_table-class} object or NULL.
#'
#' @export
#' @rdname orig_feat-methods
#' @docType methods
#'
#' @examples
#'     data(siamcat_example)
#'     data_split(siamcat_example)
setGeneric("orig_feat", function(siamcat)
    standardGeneric("orig_feat"))
#' @rdname orig_feat-methods
#' @aliases orig_feat,ANY-method
setMethod("orig_feat", "ANY", function(siamcat) {
    accessSlot(siamcat, "orig_feat")
})
# Return as-is if already a orig_feat object
#' @rdname orig_feat-methods
setMethod("orig_feat", "orig_feat", function(siamcat) {
    return(siamcat)
})
# constructor; for creating label from a list
#' @rdname label-methods
#' @aliases label
setMethod("orig_feat", "otu_table", function(siamcat) {
    return(new("orig_feat", siamcat))
})

################################################################################
#' Retrieve a \link[phyloseq]{phyloseq-class} object from object.
#'
#'
#' @usage physeq(siamcat)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a label or instance of \link[phyloseq]{phyloseq-class}.
#'
#'
#' @return The \link[phyloseq]{phyloseq-class} object or NULL.
#'
#' @export
#' @rdname physeq-methods
#' @docType methods
#'
#' @examples
#'     data(siamcat_example)
#'     physeq(siamcat_example)
setGeneric("physeq", function(siamcat)
    standardGeneric("physeq"))
#' @rdname physeq-methods
#' @aliases physeq,ANY-method
setMethod("physeq", "ANY", function(siamcat) {
    accessSlot(siamcat, "phyloseq")
})
# Return as-is if already a phyloseq object
#' @rdname physeq-methods
setMethod("physeq", "phyloseq", function(siamcat) {
    return(siamcat)
})

################################################################################
#' Retrieve a \link[phyloseq]{otu_table-class} object from object.
#'
#'
#' @usage features(siamcat)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a label or instance of \link[phyloseq]{otu_table-class} .
#'
#'
#' @return The \link[phyloseq]{otu_table-class} object or NULL.
#'
#' @export
#' @rdname features-methods
#' @docType methods
#'
#' @examples
#'     data(siamcat_example)
#'     features(siamcat_example)
setGeneric("features", function(siamcat)
    standardGeneric("features"))
#' @rdname features-methods
#' @aliases features,ANY-method
setMethod("features", "ANY", function(siamcat) {
    otu_table(physeq(siamcat))
})
# Return as-is if already a otu_table object
#' @rdname features-methods
setMethod("features", "otu_table", function(siamcat) {
    return(siamcat)
})


################################################################################
#' Retrieve a \link[phyloseq]{sample_data-class} object from object.
#'
#'
#' @usage meta(siamcat)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a label or instance of \link[phyloseq]{sample_data-class}.
#'
#'
#' @return The \link[phyloseq]{sample_data-class} object or NULL.
#'
#' @export
#' @rdname meta-methods
#' @docType methods
#'
#' @examples
#'     data(siamcat_example)
#'     meta(siamcat_example)
setGeneric("meta", function(siamcat)
    standardGeneric("meta"))
#' @rdname meta-methods
#' @aliases meta,ANY-method
setMethod("meta", "ANY", function(siamcat) {
    sample_data(physeq(siamcat), errorIfNULL = FALSE)
})
# Return as-is if already a sample_data object
#' @rdname meta-methods
setMethod("meta", "sample_data", function(siamcat) {
    return(siamcat)
})


#' Access features in siamcat@phylose@otu_table as matrix
#' @title get.features.matrix
#' @name get.features.matrix
#' @description Function to access features in siamcat@phylose@otu_table
#' @param siamcat an object of class \link{siamcat-class}t
#' @return Features as a matrix
#' @export
#' @examples
#'     data(siamcat_example)
#'     feat <- get.features.matrix(siamcat_example)
get.features.matrix <- function(siamcat) {
    return(siamcat@phyloseq@otu_table@.Data)
}

#' Access original features in siamcat@orig_feat as matrix
#' @title get.orig_feat.matrix
#' @name get.orig_feat.matrix
#' @description Function to access original features in siamcat@orig_feat
#' @param siamcat an object of class \link{siamcat-class}t
#' @return Original features as a matrix
#' @export
#' @examples
#'     data(siamcat_example)
#'     orig_feat <- get.orig_feat.matrix(siamcat_example)
get.orig_feat.matrix <- function(siamcat) {
    return(siamcat@orig_feat@.Data)
}
