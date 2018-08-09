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
accessSlot <- function(siamcat, slot, verbose=1) {
    if (!slot %in% slotNames(siamcat)) {
        # If slot is invalid, set out to NULL. Test later if this is an error.
        out = NULL
    } else if (slot == 'associations') {
        out <- siamcat@associations@assoc.results
    } else if (slot == 'filt_feat') {
        out <- siamcat@filt_feat@filt.feat
    } else if (slot == 'norm_feat'){
        out <- siamcat@norm_feat@norm.feat
    } else {
        # By elimination, must be valid. Access slot
        out = eval(parse(text = paste("siamcat@", slot, sep = "")))
    }

    # make sure that NULL is returned when the object is empty
    if (slot %in% c('eval_data', 'data_split')){
        if (length(out) == 0) {
            out <- NULL
        }
    }
    if (slot %in% c('pred_matrix', 'associations')){
        if (nrow(out) == 0){
            out <- NULL
        }
    }
    if (slot == 'model_list'){
        if (length(out@model.type) == 0){
            out <- NULL
        }
    }

    if (is.null(out) & verbose > 0) {
        # Only error regarding a NULL return value if errorIfNULL is TRUE.
        message(slot, " slot is empty.")
    }
    return(out)
}


# phyloseq = "phyloseq",

################################################################################
#' Retrieve a \link[phyloseq]{phyloseq-class} object from object.
#'
#' @usage physeq(siamcat)
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a label or instance of \link[phyloseq]{phyloseq-class}.
#' @return The \link[phyloseq]{phyloseq-class} object or NULL.
#' @export
#' @rdname physeq-methods
#' @docType methods
#' @examples
#'     data(siamcat_example)
#'     physeq(siamcat_example)
setGeneric("physeq", function(siamcat, verbose=1)
    standardGeneric("physeq"))
#' @rdname physeq-methods
#' @aliases physeq,ANY-method
setMethod("physeq", "ANY", function(siamcat, verbose=1) {
    accessSlot(siamcat, "phyloseq", verbose)
})
# Return as-is if already a phyloseq object
#' @rdname physeq-methods
setMethod("physeq", "phyloseq", function(siamcat) {
    return(siamcat)
})

################################################################################
#' Retrieve a \link[phyloseq]{otu_table-class} object from otu_table slot in
#'  the phyloseq slot in a siamcat object
#'
#' @usage orig_feat(siamcat)
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a label or instance of \link[phyloseq]{otu_table-class}.
#' @return The \link[phyloseq]{otu_table-class} object or NULL.
#' @export
#' @rdname orig_feat-methods
#' @docType methods
#' @examples
#'     data(siamcat_example)
#'     orig_feat(siamcat_example)
setGeneric("orig_feat", function(siamcat)
    standardGeneric("orig_feat"))
#' @rdname orig_feat-methods
#' @aliases orig_feat,ANY-method
setMethod("orig_feat", "ANY", function(siamcat) {
    otu_table(physeq(siamcat))
})
# Return as-is if already a otu_table object
#' @rdname orig_feat-methods
setMethod("orig_feat", "otu_table", function(siamcat) {
    return(siamcat)
})

#' Access original features in siamcat@phyloseq as matrix
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
    return(siamcat@phyloseq@otu_table@.Data)
}

################################################################################
#' Retrieve the active \link[phyloseq]{otu_table-class} out of SIAMCAT
#'
#' @usage features(siamcat)
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a label or instance of \link[phyloseq]{otu_table-class} .
#' @return The \link[phyloseq]{otu_table-class} object or NULL.
#' @export
#' @rdname features-methods
#' @docType methods
#' @examples
#'     data(siamcat_example)
#'     features(siamcat_example)
setGeneric("features", function(siamcat)
    standardGeneric("features"))
#' @rdname features-methods
#' @aliases features,ANY-method
setMethod("features", "ANY", function(siamcat) {
    if (!is.null(norm_feat(siamcat, verbose=0))){
        return(norm_feat(siamcat))
    } else if (!is.null(filt_feat(siamcat, verbose=0))){
        return(filt_feat(siamcat))
    } else {
        return(orig_feat(siamcat))
    }
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
    return(features(siamcat)@.Data)
}

################################################################################
#' Retrieve a \link[phyloseq]{sample_data-class} object from object.
#'
#' @usage meta(siamcat)
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a label or instance of \link[phyloseq]{sample_data-class}.
#' @return The \link[phyloseq]{sample_data-class} object or NULL.
#' @export
#' @rdname meta-methods
#' @docType methods
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

################################################################################
#' Retrieve a \link{label-class} object from object.
#'
#' @usage label(siamcat)
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a label or instance of \link{label-class} or a list.
#' @return The \link{label-class} object or NULL.
#' @export
#' @rdname label-methods
#' @docType methods
#' @examples
#'     data(siamcat_example)
#'     label(siamcat_example)
setGeneric("label", function(siamcat, verbose=1)
    standardGeneric("label"))
#' @rdname label-methods
#' @aliases label,ANY-method
setMethod("label", "ANY", function(siamcat, verbose=1) {
    accessSlot(siamcat, "label", verbose)
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

###############################################################################
#' Retrieve filtered features form object
#'
#' @usage filt_feat(siamcat)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains an object in the filt_feat slot
#'
#' @return The filtered feature matrix or NULL.
#'
#' @export
#' @rdname filt_feat-methods
#' @docType methods
#'
#' @examples
#'     data(siamcat_example)
#'     filt_feat(siamcat_example)
setGeneric("filt_feat", function(siamcat, verbose=1)
    standardGeneric("filt_feat"))
#' @rdname filt_feat-methods
#' @aliases filt_feat,ANY-method
setMethod("filt_feat", "ANY", function(siamcat, verbose=1) {
    temp <- siamcat@filt_feat@filt.feat
    if (is.null(temp) & verbose > 0) {
        message("Filtered features slot is empty.")
    }
    return(temp)
})

#' Access features in siamcat@filt_feat@filt.feat as matrix
#' @title get.filt_features.matrix
#' @name get.filt_features.matrix
#' @description Function to access features in siamcat@filt_feat@filt.feat
#' @param siamcat an object of class \link{siamcat-class}
#' @return Filtered features as a matrix
#' @export
#' @examples
#'     data(siamcat_example)
#'     feat <- get.filt_features.matrix(siamcat_example)
get.filt_feat.matrix <- function(siamcat) {
    return(filt_feat(siamcat)@.Data)
}


###############################################################################
#' Retrieve the list of filtering parameters from object.
#'
#' @usage filt_params(siamcat)
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a filt_feat object
#' @return The list of filtering parameters or NULL.
#' @export
#' @rdname filt_params-methods
#' @docType methods
#'
#' @examples
#'     data(siamcat_example)
#'     filt_params(siamcat_example)
setGeneric("filt_params", function(siamcat, verbose=1)
    standardGeneric("filt_params"))
#' @rdname filt_param-methods
#' @aliases filt_param,ANY-method
setMethod("filt_params", "ANY", function(siamcat, verbose=1) {
    temp <- siamcat@filt_feat@filt.param
    if (length(temp) == 0) temp <- NULL
    if (is.null(temp) & verbose > 0) {
        message("Filtering paramters slot is empty.")
    }
    return(temp)
})


###############################################################################
#' Retrieve associations from object.
#'
#' @usage associations(siamcat)
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains an object in the associations slot
#' @return The results of the association testing or NULL.
#' @export
#' @rdname associations-methods
#' @docType methods
#' @examples
#'     data(siamcat_example)
#'     associations(siamcat_example)
setGeneric("associations", function(siamcat, verbose=1)
    standardGeneric("associations"))
#' @rdname associations-methods
#' @aliases associations,ANY-method
setMethod("associations", "ANY", function(siamcat, verbose=1) {
    accessSlot(siamcat, "associations", verbose)
})
# constructor; for creating an associations object from a list
#' @rdname associations-methods
#' @aliases associations
setMethod("associations", "associations", function(siamcat) {
    return(new("associations",
        siamcat))
})

###############################################################################
#' Retrieve parameters of association testing from object.
#'
#' @usage assoc_param(siamcat)
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains an object in the associations slot
#' @return The parameters of the assocation testing or NULL
#' @export
#' @rdname associations-methods
#' @docType methods
#' @examples
#'     data(siamcat_example)
#'     assoc_param(siamcat_example)
setGeneric("assoc_param", function(siamcat, verbose=1)
    standardGeneric("assoc_param"))
#' @rdname assoc-methods
#' @aliases assoc_param,ANY-method
setMethod("assoc_param", "ANY", function(siamcat, verbose=1) {
    temp <- siamcat@associations@assoc.param
    if (length(temp) == 0) temp <- NULL
    if (is.null(temp) & verbose > 0){
        message("Parameters for assocations testing are empty")
    }
    return(temp)
})

###############################################################################
#' Retrieve normalized features form object
#'
#' @usage norm_feat(siamcat)
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains an object in the norm_feat slot
#' @return The normalized feature matrix or NULL.
#' @export
#' @rdname norm_feat-methods
#' @docType methods
#' @examples
#'     data(siamcat_example)
#'     norm_feat(siamcat_example)
setGeneric("norm_feat", function(siamcat, verbose=1)
    standardGeneric("norm_feat"))
#' @rdname norm_feat-methods
#' @aliases norm_feat,ANY-method
setMethod("norm_feat", "ANY", function(siamcat, verbose=1) {
    temp <- siamcat@norm_feat@norm.feat
    if (length(temp) == 0) temp <- NULL
    if (is.null(temp) & verbose > 0){
        message("Normalized features are empty")
    }
    return(temp)
})

###############################################################################
#' Retrieve the list of normalization parameters from object.
#'
#' @usage norm_params(siamcat)
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a norm_feat
#' @return The list of normalization parameters or NULL.
#' @export
#' @rdname norm_params-methods
#' @docType methods
#'
#' @examples
#'     data(siamcat_example)
#'     norm_params(siamcat_example)
setGeneric("norm_params", function(siamcat, verbose=1)
    standardGeneric("norm_params"))
#' @rdname norm_params-methods
#' @aliases norm_params,ANY-method
setMethod("norm_params", "ANY", function(siamcat, verbose=1) {
    temp <-  siamcat@norm_feat@norm.param
    if (length(temp) == 0) temp <- NULL
    if (is.null(temp) & verbose > 0){
        message("Parameters for normalization are empty")
    }
    return(temp)
})

################################################################################
#' Retrieve a \link{data_split-class} object from object.
#'
#' @usage data_split(siamcat)
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a label or instance of \link{data_split-class} or a list.
#' @return The \link{data_split-class} object or NULL.
#' @export
#' @rdname data_split-methods
#' @docType methods
#' @examples
#'     data(siamcat_example)
#'     data_split(siamcat_example)
setGeneric("data_split", function(siamcat, verbose=1)
    standardGeneric("data_split"))
#' @rdname data_split-methods
#' @aliases data_split,ANY-method
setMethod("data_split", "ANY", function(siamcat, verbose=1) {
    accessSlot(siamcat, "data_split", verbose)
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

###############################################################################
#' Retrieve \link{model_list-class} from object.
#'
#' @usage model_list(siamcat)
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a model_list or instance of \link{model_list-class}.
#' @return The \link{model_list-class} object or NULL.
#' @export
#' @rdname model_list-methods
#' @docType methods
#' @examples
#'     data(siamcat_example)
#'     model_list(siamcat_example)
setGeneric("model_list", function(siamcat, verbose=1)
    standardGeneric("model_list"))
#' @rdname model_list-methods
#' @aliases model_list,ANY-method
setMethod("model_list", "ANY", function(siamcat, verbose=1) {
    accessSlot(siamcat, "model_list", verbose)
})
# Return as-is if already a model_list object
#' @rdname model_list-methods
setMethod("model_list", "model_list", function(siamcat) {
    return(siamcat)
})

###############################################################################
#' Retrieve list of models from object.
#'
#' @usage models(siamcat)
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a model_list or instance of \link{model_list-class}.
#' @return The list of models or NULL.
#' @export
#' @rdname models-methods
#' @docType methods
#' @examples
#'     data(siamcat_example)
#'     models(siamcat_example)
setGeneric("models", function(siamcat, verbose=1)
    standardGeneric("models"))
#' @rdname models-methods
#' @aliases models,ANY-method
setMethod("models", "ANY", function(siamcat, verbose=1) {
    temp <- siamcat@model_list@models
    if (length(temp) == 0) temp <- NULL
    if (is.null(temp) & verbose > 0){
        message("Model slot are empty")
    }
    return(temp)
})
# Return list of models if a model_list object
#' @rdname model_list-methods
setMethod("models", "model_list", function(siamcat) {
    temp <- siamcat@models
    if (length(temp) == 0) temp <- NULL
    if (is.null(temp) & verbose > 0){
        message("Model slot are empty")
    }
    return(temp)
})

###############################################################################
#' Retrieve model_type from object.
#'
#' @usage model_type(siamcat)
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a model_list or instance of \link{model_list-class}.
#' @return The string describing type of model used or NULL.
#' @export
#' @rdname model_type-methods
#' @docType methods
#' @examples
#'     data(siamcat_example)
#'     model_type(siamcat_example)
setGeneric("model_type", function(siamcat, verbose=1)
    standardGeneric("model_type"))
#' @rdname model_type-methods
#' @aliases model_type,ANY-method
setMethod("model_type", "ANY", function(siamcat, verbose=1) {
    temp <- siamcat@model_list@model.type
    if (length(temp) == 0) temp <- NULL
    if (is.null(temp) & verbose > 0){
        message("Model type are empty")
    }
    return(temp)
})
# Return model type if a model_list object
#' @rdname model_list-methods
setMethod("model_type", "model_list", function(siamcat, verbose=1) {
    temp <- siamcat@model_type
    if (length(temp) == 0) temp <- NULL
    if (is.null(temp) & verbose > 0){
        message("Model type are empty")
    }
    return(temp)
})

###############################################################################
#' Retrieve pred_matrix from object.
#'
#' @usage pred_matrix(siamcat)
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a pred_matrix
#' @return The pred_matrix matrix or NULL.
#' @export
#' @rdname pred_matrix-methods
#' @docType methods
#' @examples
#'     data(siamcat_example)
#'     pred_matrix(siamcat_example)
setGeneric("pred_matrix", function(siamcat, verbose=1)
    standardGeneric("pred_matrix"))
#' @rdname pred_matrix-methods
#' @aliases pred_matrix,ANY-method
setMethod("pred_matrix", "ANY", function(siamcat, verbose=1) {
    accessSlot(siamcat, "pred_matrix", verbose)
})
# constructor; for creating pred_matrix from a matrix
#' @rdname pred_matrix-methods
#' @aliases pred_matrix
setMethod("pred_matrix", "matrix", function(siamcat) {
    return(new("pred_matrix",
        siamcat))
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
setGeneric("eval_data", function(siamcat, verbose=1)
    standardGeneric("eval_data"))
#' @rdname eval_data-methods
#' @aliases eval_data,ANY-method
setMethod("eval_data", "ANY", function(siamcat, verbose=1) {
    accessSlot(siamcat, "eval_data", verbose)
})
# constructor; for creating eval_data from a list
#' @rdname eval_data-methods
#' @aliases eval_data
setMethod("eval_data", "list", function(siamcat) {
    return(new("eval_data",
        siamcat))
})
