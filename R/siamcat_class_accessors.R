#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

################################################################################
#' Universal slot accessor function for siamcat-class.
#'
#' This function is used internally by many accessors.
#'
#' @usage accessSlot(siamcat, slot, verbose=1)
#'
#' @param siamcat an object of \link{siamcat-class}.
#'
#' @param slot A character string indicating the slot (not data class)
#'     of the component data type that is desired.
#'
#' @param verbose If the slot is empty, should a message be printed? values
#'      can be either 0 (no output) or 1 (print message)
#'
#' @return Returns the component object specified by the argument \code{slot}.
#'     Returns NULL if slot does not exist.
#'
#' @export
#' @keywords internal
#' @examples #
#' data(siamcat_example)
#' accessSlot(siamcat_example, "label")
#' accessSlot(siamcat_example, "model_list")
accessSlot <- function(siamcat, slot, verbose=1) {
    if (!slot %in% slotNames(siamcat)) {
        stop('SIAMCAT object does not contain a slot: ', slot)
    } else if (slot == 'associations') {
        out <- siamcat@associations$assoc.results
    } else {
        # By elimination, must be valid. Access slot
        out <- eval(parse(text = paste("siamcat@", slot, sep = "")))
    }

    # make sure that NULL is returned when the object is empty
    if (slot %in% c(
        'eval_data', 'data_split', 'model_list', 'norm_feat', 'filt_feat')){
        if (length(out) == 0) {
            out <- NULL
        }
    }
    if (slot %in% c('pred_matrix')){
        if (nrow(out) == 0){
            out <- NULL
        }
    }

    if (is.null(out) & verbose > 0) {
        # Only error regarding a NULL return value if errorIfNULL is TRUE.
        message(slot, " slot is empty.")
    }
    return(out)
}

################################################################################
#' Retrieve a \link[phyloseq]{phyloseq-class} object from object.
#'
#' @usage physeq(siamcat, verbose=1)
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'     that contains a label or instance of \link[phyloseq]{phyloseq-class}.
#'
#' @param verbose If the slot is empty, should a message be printed? values
#'      can be either 0 (no output) or 1 (print message)
#'
#' @return The \link[phyloseq]{phyloseq-class} object or NULL.
#' @export
#' @keywords internal
#' @rdname physeq-methods
#' @docType methods
#' @examples
#' data(siamcat_example)
#' physeq(siamcat_example)
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
#' the phyloseq slot in a siamcat object
#'
#' @usage orig_feat(siamcat)
#' @param siamcat (Required). An instance of \link{siamcat-class}
#' that contains a label or instance of \link[phyloseq]{otu_table-class}.
#' @return The \link[phyloseq]{otu_table-class} object or NULL.
#' @export
#' @keywords internal
#' @rdname orig_feat-methods
#' @docType methods
#' @examples
#' data(siamcat_example)
#' temp <- orig_feat(siamcat_example)
setGeneric("orig_feat", function(siamcat)
    standardGeneric("orig_feat"))
#' @rdname orig_feat-methods
#' @aliases orig_feat
setMethod("orig_feat", "siamcat", function(siamcat) {
    otu_table(physeq(siamcat))
})
#' @rdname orig_feat-methods
#' @aliases orig_feat
setMethod("orig_feat", "otu_table", function(siamcat) {
    return(siamcat)
})

################################################################################
#' @title Retrieve the original features from a SIAMCAT object
#'
#' @description Function to retrieve the original features from a SIAMCAT
#' object
#'
#' @usage get.orig_feat.matrix(siamcat)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#'
#' @details The function returns the original features as matrix.
#'
#' @return A matrix containing the original features
#'
#' @export
#'
#' @examples
#' data(siamcat_example)
#' feat.original <- get.orig_feat.matrix(siamcat_example)
#' feat.original[1:3, 1:3]
get.orig_feat.matrix <- function(siamcat) {
    return(siamcat@phyloseq@otu_table@.Data)
}

################################################################################
#' @title Retrieve the metadata from a SIAMCAT object
#'
#' @description Retrieve the metadata from a SIAMCAT object
#'
#' @details This function will retrieve the metadata from a
#' SIAMCAT object. The metadata is a object of the
#' \link[phyloseq]{sample_data-class}.
#'
#' @usage meta(siamcat)
#'
#' @param siamcat (Required). A \link{siamcat-class} object
#'
#' @return The metadata as \link[phyloseq]{sample_data-class} object
#'
#' @export
#'
#' @rdname meta-methods
#'
#' @docType methods
#'
#' @examples
#' data(siamcat_example)
#' temp <- meta(siamcat_example)
#' head(temp)
setGeneric("meta", function(siamcat)
    standardGeneric("meta"))
#' @rdname meta-methods
#' @aliases meta
setMethod("meta", "siamcat", function(siamcat) {
    sample_data(physeq(siamcat), errorIfNULL = FALSE)
})
#' @rdname meta-methods
#' @aliases meta
setMethod("meta", "sample_data", function(siamcat) {
    return(siamcat)
})

################################################################################
#' @title Retrieve the label from a SIAMCAT object
#'
#' @description Retrieve the label from a SIAMCAT object
#'
#' @details This function will retrieve the label information from a
#' SIAMCAT object. The label will contain three entries: \itemize{
#' \item \code{label}: The label as named vector, in which the classes are
#' encoded numerically
#' \item \code{info}: Information about the different classes
#' \item \code{type}: What kind of label is it?
#' }
#'
#' @usage label(siamcat, verbose=1)
#'
#' @param siamcat (Required). A \link{siamcat-class} object
#'
#' @param verbose integer, if the slot is empty, should a message be printed?
#'      values can be either \code{0} (no output) or \code{1} (print message)
#'
#' @return The label or NULL.
#'
#' @export
#'
#' @rdname label-methods
#'
#' @docType methods
#'
#' @examples
#' data(siamcat_example)
#' temp <- label(siamcat_example)
#' head(temp$label)
#' temp$info
#' temp$type
setGeneric("label", function(siamcat, verbose=1)
    standardGeneric("label"))
#' @rdname label-methods
#' @aliases label
setMethod("label", "siamcat", function(siamcat, verbose=1) {
    accessSlot(siamcat, "label", verbose)
})

###############################################################################
#' @title Retrieve the information stored in the \code{filt_feat} slot within
#' a SIAMCAT object
#'
#' @description Function to retrieve the information stored in the
#' \code{filt_feat} slot within a SIAMCAT object
#'
#' @usage filt_feat(siamcat, verbose=1)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#' that contains filtered features
#'
#' @param verbose integer, if the slot is empty, should a message be printed?
#'      values can be either \code{0} (no output) or \code{1} (print message)
#'
#' @details The function will return a list containing the information stored
#' in the \code{filt_feat} slot of a SIAMCAT object. This list contains:
#' \itemize{
#' \item \code{filt.feat} - filtered features as matrix, see
#' \link{get.filt_feat.matrix}
#' \item \code{filt.param} - parameters used for feature filtering, see
#' \link{get.filt_feat.matrix}
#' }
#'
#' @return The list stored in the \code{filt_feat} slot of the SIAMCAT object
#' or \code{NULL}
#'
#' @export
#' @keywords internal
#'
#' @rdname filt_feat-methods
#'
#' @docType methods
#'
#' @examples
#' data(siamcat_example)
#' temp <- filt_feat(siamcat_example)
#' names(temp)
setGeneric("filt_feat", function(siamcat, verbose=1)
    standardGeneric("filt_feat"))
#' @rdname filt_feat-methods
#' @aliases filt_feat
setMethod("filt_feat", "siamcat", function(siamcat, verbose=1) {
    accessSlot(siamcat, "filt_feat", verbose)
})


################################################################################
#' @title Retrieve the filtered features from a SIAMCAT object
#'
#' @description Function to retrieve the filtered features from a SIAMCAT
#' object
#'
#' @usage get.filt_feat.matrix(siamcat)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#' containing filtered features
#'
#' @details The function returns the filtered features as matrix.
#' See \link{filter.features} for more details.
#'
#' @return A matrix containing the filtered features
#'
#' @export
#'
#' @examples
#' data(siamcat_example)
#' feat.filt <- get.filt_feat.matrix(siamcat_example)
#' feat.filt[1:3, 1:3]
get.filt_feat.matrix <- function(siamcat) {
    return(filt_feat(siamcat)$filt.feat)
}


################################################################################
#' @title Retrieve the list of parameters for feature filtering from a
#' SIAMCAT object
#'
#' @description Function to retrieve the list of parameters for feature
#' filtering
#'
#' @usage filt_params(siamcat, verbose=1)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#' containing filtered features
#'
#' @param verbose integer, if the slot is empty, should a message be printed?
#'      values can be either \code{0} (no output) or \code{1} (print message)
#'
#' @details The function returns the list of feature filtering parameters.
#' See \link{filter.features} for more details.
#'
#' @return A list of feature filtering parameters or \code{NULL}
#'
#' @export
#'
#' @rdname filt_params-methods
#'
#' @docType methods
#'
#' @examples
#' data(siamcat_example)
#' temp <- filt_params(siamcat_example)
#' names(temp)
setGeneric("filt_params", function(siamcat, verbose=1)
    standardGeneric("filt_params"))
#' @rdname filt_params-methods
#' @aliases filt_params
setMethod("filt_params", "siamcat", function(siamcat, verbose=1) {
    temp <- siamcat@filt_feat$filt.param
    if (length(temp) == 0) temp <- NULL
    if (is.null(temp) & verbose > 0) {
        message("Filtering paramters slot is empty.")
    }
    return(temp)
})

################################################################################
#' @title Retrieve the results of association testing from a SIAMCAT object
#'
#' @description Function to retrieve the results of association testing
#'
#' @usage associations(siamcat, verbose=1)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#' containing the results of association testing
#'
#' @param verbose integer, if the slot is empty, should a message be printed?
#'      values can be either \code{0} (no output) or \code{1} (print message)
#'
#' @details The function returns the results of the association testing
#' procedure as dataframe. See \link{check.associations} for more details.
#'
#' @return A \code{data.frame} of association testing results or \code{NULL}
#'
#' @export
#'
#' @rdname associations-methods
#'
#' @docType methods
#'
#' @examples
#' data(siamcat_example)
#' temp <- associations(siamcat_example)
#' head(temp)
setGeneric("associations", function(siamcat, verbose=1)
    standardGeneric("associations"))
#' @rdname associations-methods
#' @aliases associations
setMethod("associations", "siamcat", function(siamcat, verbose=1) {
    accessSlot(siamcat, "associations", verbose)
})

################################################################################
#' @title Retrieve the list of parameters for association testing from
#' a SIAMCAT object
#'
#' @description Function to retrieve the list of parameters for
#' association testing
#'
#' @usage assoc_param(siamcat, verbose=1)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#' containing the results from association testing
#'
#' @param verbose integer, if the slot is empty, should a message be printed?
#'      values can be either \code{0} (no output) or \code{1} (print message)
#'
#' @details The function returns the list of parameters used in association
#' testing. See \link{check.associations} for more details.
#'
#' @return A list of parameters for association testing or \code{NULL}
#'
#' @export
#'
#' @rdname assoc_param-methods
#'
#' @docType methods
#'
#' @examples
#' data(siamcat_example)
#' temp <- assoc_param(siamcat_example)
#' names(temp)
setGeneric("assoc_param", function(siamcat, verbose=1)
    standardGeneric("assoc_param"))
#' @rdname assoc_param-methods
#' @aliases assoc_param_param
setMethod("assoc_param", "siamcat", function(siamcat, verbose=1) {
    temp <- siamcat@associations$assoc.param
    if (length(temp) == 0) temp <- NULL
    if (is.null(temp) & verbose > 0){
        message("Parameters for associations testing are empty")
    }
    return(temp)
})

###############################################################################
#' @title Retrieve the information stored in the \code{norm_feat} slot within
#' a SIAMCAT object
#'
#' @description Function to retrieve the information stored in the
#' \code{norm_feat} slot within a SIAMCAT object
#'
#' @usage norm_feat(siamcat, verbose=1)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#' that contains normalized features
#'
#' @param verbose integer, if the slot is empty, should a message be printed?
#'      values can be either \code{0} (no output) or \code{1} (print message)
#'
#' @details The function will return a list containing the information stored
#' in the \code{norm_feat} slot of a SIAMCAT object. This list contains:
#' \itemize{
#' \item \code{norm.feat} - normalized features as matrix, see
#' \link{get.norm_feat.matrix}
#' \item \code{norm.param} - parameters used for normalization, see
#' \link{normalize.features}
#' }
#'
#' @return The list stored in the \code{norm_feat} slot of the SIAMCAT object
#' or \code{NULL}
#'
#' @export
#' @keywords internal
#'
#' @rdname norm_feat-methods
#'
#' @docType methods
#'
#' @examples
#' data(siamcat_example)
#' temp <- norm_feat(siamcat_example)
#' names(temp)
setGeneric("norm_feat", function(siamcat, verbose=1)
    standardGeneric("norm_feat"))
#' @rdname norm_feat-methods
#' @aliases norm_feat
setMethod("norm_feat", "siamcat", function(siamcat, verbose=1) {
    accessSlot(siamcat, "norm_feat", verbose)
})


################################################################################
#' @title Retrieve the normalized features from a SIAMCAT object
#'
#' @description Function to retrieve the normalized features from a SIAMCAT
#' object
#'
#' @usage get.norm_feat.matrix(siamcat)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#' containing normalized features
#'
#' @details The function returns the normalized features as matrix.
#' See \link{normalize.features} for more details.
#'
#' @return A matrix containing the normalized features
#'
#' @export
#'
#' @examples
#' data(siamcat_example)
#' feat.norm <- get.norm_feat.matrix(siamcat_example)
#' feat.norm[1:3, 1:3]
get.norm_feat.matrix <- function(siamcat) {
    return(norm_feat(siamcat)$norm.feat)
}

################################################################################
#' @title Retrieve the list of parameters for feature normalization from a
#' SIAMCAT object
#'
#' @description Function to retrieve the list of parameters for feature
#' normalization
#'
#' @usage norm_params(siamcat, verbose=1)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#' containing normalized features
#'
#' @param verbose integer, if the slot is empty, should a message be printed?
#'      values can be either \code{0} (no output) or \code{1} (print message)
#'
#' @details The function returns the list of normalization parameters used in
#' the feature normalization procedure. See \link{normalize.features} for
#' more details.
#'
#' @return A list of normalization parameters or \code{NULL}
#'
#' @export
#'
#' @rdname norm_params-methods
#'
#' @docType methods
#'
#' @examples
#' data(siamcat_example)
#' temp <- norm_params(siamcat_example)
#' names(temp)
setGeneric("norm_params", function(siamcat, verbose=1)
    standardGeneric("norm_params"))
#' @rdname norm_params-methods
#' @aliases norm_params
setMethod("norm_params", "siamcat", function(siamcat, verbose=1) {
    temp <-  siamcat@norm_feat$norm.param
    if (length(temp) == 0) temp <- NULL
    if (is.null(temp) & verbose > 0){
        message("Parameters for normalization are empty")
    }
    return(temp)
})

################################################################################
#' @title Retrieve the data split from a SIAMCAT object
#'
#' @description Function to retrieve the data split stored in the
#' \code{data_split} slot within a SIAMCAT object
#'
#' @usage data_split(siamcat, verbose=1)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#' containing a data split
#'
#' @param verbose integer, if the slot is empty, should a message be printed?
#'      values can be either \code{0} (no output) or \code{1} (print message)
#'
#' @details The function returns a list containing information about the data
#' split. See \link{create.data.split} for more details.
#'
#' @return A list containing the data split information or \code{NULL}
#'
#' @export
#'
#' @rdname data_split-methods
#'
#' @docType methods
#'
#' @examples
#' data(siamcat_example)
#' temp <- data_split(siamcat_example)
#' names(temp)
setGeneric("data_split", function(siamcat, verbose=1)
    standardGeneric("data_split"))
#' @rdname data_split-methods
#' @aliases data_split
setMethod("data_split", "siamcat", function(siamcat, verbose=1) {
    accessSlot(siamcat, "data_split", verbose)
})

###############################################################################
#' @title Retrieve the information stored in the \code{model_list} slot within
#' a SIAMCAT object
#'
#' @description Function to retrieve the information stored in the
#' \code{model_list} slot within a SIAMCAT object
#'
#' @usage model_list(siamcat, verbose=1)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#' that contains trained models
#'
#' @param verbose integer, if the slot is empty, should a message be printed?
#'      values can be either \code{0} (no output) or \code{1} (print message)
#'
#' @details The function will return a list containing the information stored
#' in the \code{model_list} slot of a SIAMCAT object. This list contains:
#' \itemize{
#' \item \code{models} - list of trained models
#' \item \code{model_type} - machine learning method used for training
#' \item \code{feature_type} - string describing on which type of features the
#' models were trained
#' }
#'
#' @return The list stored in the \code{model_list} slot of the SIAMCAT object
#' or \code{NULL}
#'
#' @export
#' @keywords internal
#'
#' @rdname model_list-methods
#'
#' @docType methods
#'
#' @examples
#' data(siamcat_example)
#' temp <- model_list(siamcat_example)
#' names(temp)
setGeneric("model_list", function(siamcat, verbose=1)
    standardGeneric("model_list"))
#' @rdname model_list-methods
#' @aliases model_list
setMethod("model_list", "siamcat", function(siamcat, verbose=1) {
    accessSlot(siamcat, "model_list", verbose)
})


###############################################################################
#' @title Retrieve list of trained models from a SIAMCAT object
#'
#' @description Function to retrieve the list of trained models
#'
#' @usage models(siamcat, verbose=1)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#' that contains trained models
#'
#' @param verbose integer, if the slot is empty, should a message be printed?
#'      values can be either \code{0} (no output) or \code{1} (print message)
#'
#' @details The function extracts the list of trained models.
#'
#' @return The list of models or \code{NULL}
#'
#' @export
#'
#' @rdname models-methods
#'
#' @docType methods
#'
#' @examples
#' data(siamcat_example)
#' temp <- models(siamcat_example)
#' temp[[1]]
setGeneric("models", function(siamcat, verbose=1)
    standardGeneric("models"))
#' @rdname models-methods
#' @aliases models
setMethod("models", "siamcat", function(siamcat, verbose=1) {
    temp <- siamcat@model_list$models
    if (length(temp) == 0) temp <- NULL
    if (is.null(temp) & verbose > 0){
        message("Model slot is empty")
    }
    return(temp)
})

###############################################################################
#' @title Retrieve the machine learning method from a SIAMCAT object
#'
#' @description Function to retrieve information on which type of machine
#' learning method was used for model training
#'
#' @usage model_type(siamcat, verbose=1)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#' that contains trained models
#'
#' @param verbose integer, if the slot is empty, should a message be printed?
#'      values can be either \code{0} (no output) or \code{1} (print message)
#'
#' @details The function extracts the information on which type of machine
#' learning method was used for model training.
#'
#' @return The string describing the machine learning method or \code{NULL}
#'
#' @export
#'
#' @rdname model_type-methods
#'
#' @docType methods
#'
#' @examples
#' data(siamcat_example)
#' model_type(siamcat_example)
setGeneric("model_type", function(siamcat, verbose=1)
    standardGeneric("model_type"))
#' @rdname model_type-methods
#' @aliases model_type
setMethod("model_type", "siamcat", function(siamcat, verbose=1) {
    temp <- siamcat@model_list$model.type
    if (length(temp) == 0) temp <- NULL
    if (is.null(temp) & verbose > 0){
        message("Model type is empty")
    }
    return(temp)
})


###############################################################################
#' @title Retrieve the feature type used for model training from a SIAMCAT
#' object
#'
#' @description Function to retrieve information on which type of features the
#' models were trained
#'
#' @usage feature_type(siamcat, verbose=1)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#' that contains trained models
#'
#' @param verbose integer, if the slot is empty, should a message be printed?
#'      values can be either \code{0} (no output) or \code{1} (print message)
#'
#' @details The function extracts the information on which type of features
#' the models were trained.
#'
#' @return The string describing type of feature used for the model training
#' or \code{NULL}
#'
#' @export
#'
#' @rdname feature_type-methods
#'
#' @docType methods
#'
#' @examples
#' data(siamcat_example)
#' feature_type(siamcat_example)
setGeneric("feature_type", function(siamcat, verbose=1)
    standardGeneric("feature_type"))
#' @rdname feature_type-methods
#' @aliases feature_type
setMethod("feature_type", "siamcat", function(siamcat, verbose=1) {
    temp <- siamcat@model_list$feature.type
    if (length(temp) == 0) temp <- NULL
    if (is.null(temp) & verbose > 0){
        message("Feature type is empty")
    }
    return(temp)
})

###############################################################################
#' @title Retrieve the weight matrix from a SIAMCAT object
#'
#' @description Function to retrieve the feature weights from a SIAMCAT object
#'
#' @usage weight_matrix(siamcat, verbose=1)
#'
#' @param siamcat (Required). An instance of \link{siamcat-class}
#' that contains trained models
#'
#' @param verbose integer, if the slot is empty, should a message be printed?
#'      values can be either \code{0} (no output) or \code{1} (print message)
#'
#' @details The function extracts the feature weights from all trained models
#' acorss all cross-validation folds and repeats.
#'
#' @return A matrix containing the feature weights or \code{NULL}
#'
#' @export
#'
#' @rdname weight_matrix-methods
#'
#' @docType methods
#'
#' @examples
#' data(siamcat_example)
#' temp <- weight_matrix(siamcat_example)
#' temp[1:3, 1:3]

setGeneric("weight_matrix", function(siamcat, verbose=1)
    standardGeneric("weight_matrix"))
#' @rdname weight_matrix-methods
#' @aliases weight_matrix
setMethod("weight_matrix", "siamcat", function(siamcat, verbose=1) {
    temp <- siamcat@model_list$models
    if (length(temp) == 0){
        if (verbose > 0){
            message('Weight matrix is empty')
        }
        return(NULL)
    }
    if (is(temp[[1]], "WrappedModel")){
        stop("Attention:\n",
                "This SIAMCAT object seems to have been constructed with ",
                "version 1.x, based on 'mlr'.\nYour current SIAMCAT version ",
                "has been upgraded to use 'mlr3' internally.\nPlease consider ",
                "re-training your SIAMCAT object or downgrading your SIAMCAT ",
                "version in order to continue.")

    }

    feat.type <- feature_type(siamcat)
    if (feat.type == 'original'){
        feat <- get.orig_feat.matrix(siamcat)
    } else if (feat.type == 'filtered'){
        feat <- get.filt_feat.matrix(siamcat)
    } else if (feat.type == 'normalized'){
        feat <- get.norm_feat.matrix(siamcat)
    }

    weight.mat <- matrix(NA, nrow=nrow(feat),
        ncol=length(temp), dimnames=list(rownames(feat),
        paste0('Model_', seq_along(temp))))
    for (i in seq_along(temp)){
        m.idx <- match(names(temp[[i]]$features),
            make.names(rownames(weight.mat)))
        weight.mat[m.idx, i] <- temp[[i]]$features
    }
    weight.mat[is.na(weight.mat)] <- 0

    return(weight.mat)
})


###############################################################################
#' @title Retrieve the matrix of feature weights from a SIAMCAT object
#'
#' @description Function to extract the feature weights from a SIAMCAT object
#'
#' @usage feature_weights(siamcat, verbose=1)
#'
#' @param siamcat (Required). A \link{siamcat-class} object
#' that contains trained models
#'
#' @param verbose integer, if the slot is empty, should a message be printed?
#'      values can be either \code{0} (no output) or \code{1} (print message)
#'
#' @details The function extracts the weight matrix from all trained models
#' (see \link{weight_matrix}) and computes several metrics on the feature
#' weights: \itemize{
#' \item \code{mean.weight} - mean weight across trained models
#' \item \code{median.weight} - median weight across trained models
#' \item \code{sd.weight} - standard deviation of the weight across trained
#' models
#' \item \code{mean.rel.weight} - mean \strong{relative} weight across
#' trained models (each model is normalized by the absolute of all weights)
#' \item \code{median.rel.weight} - median \strong{relative} weight across
#' trained models
#' \item \code{sd.rel.weight} - standard deviation of the \strong{relative}
#' weight across trained models
#' \item \code{percentage} - percentage of models in which this feature was
#' selected (i.e. non-zero)}
#'
#' @return A dataframe containing mean/median feature weight and
#' additional info or \code{NULL}
#'
#' @export
#'
#' @rdname feature_weights-methods
#'
#' @docType methods
#'
#' @examples
#' data(siamcat_example)
#' temp <- feature_weights(siamcat_example)
#' head(temp)
setGeneric("feature_weights", function(siamcat, verbose=1)
    standardGeneric("feature_weights"))
#' @rdname feature_weights-methods
#' @aliases feature_weights
setMethod("feature_weights", "siamcat", function(siamcat, verbose=1) {
    feat.weights <- weight_matrix(siamcat, verbose=0)
    if (is.null(feat.weights)) {
        if(verbose > 0){
            message("Feature weights are empty")
        }
        return(NULL)
    }
    feat.weights <- data.frame(
        mean.weight=rowMeans(feat.weights, na.rm=TRUE),
        median.weight=rowMedians(feat.weights, na.rm=TRUE),
        sd.weight=rowSds(feat.weights, na.rm=TRUE),
        mean.rel.weight=rowMeans(t(t(feat.weights)/
            colSums(abs(feat.weights), na.rm=TRUE)), na.rm=TRUE),
        median.rel.weight=rowMedians(t(t(feat.weights)/
            colSums(abs(feat.weights), na.rm=TRUE)), na.rm=TRUE),
        sd.rel.weight=rowSds(t(t(feat.weights)/
            colSums(abs(feat.weights), na.rm=TRUE)), na.rm=TRUE),
        percentage=rowMeans(feat.weights != 0, na.rm=TRUE)
    )
    return(feat.weights)
})


###############################################################################
#' @title Retrieve the prediction matrix from a SIAMCAT object
#'
#' @description Function to retrieve the prediction matrix from a SIAMCAT
#' object
#'
#' @usage pred_matrix(siamcat, verbose=1)
#'
#' @param siamcat (Required). A \link{siamcat-class} object
#'     that contains a prediction matrix
#'
#' @param verbose integer, if the slot is empty, should a message be printed?
#'      values can be either \code{0} (no output) or \code{1} (print message)
#'
#' @details The functions returns a matrix containing the predictions for all
#' samples across the different cross-validation repeats. See
#' \link{make.predictions} for more information.
#'
#' @return A matrix containing predictions or \code{NULL}
#'
#' @export
#'
#' @rdname pred_matrix-methods
#'
#' @docType methods
#'
#' @examples
#' data(siamcat_example)
#' temp <- pred_matrix(siamcat_example)
#' head(temp)
setGeneric("pred_matrix", function(siamcat, verbose=1)
    standardGeneric("pred_matrix"))
#' @rdname pred_matrix-methods
#' @aliases pred_matrix
setMethod("pred_matrix", "siamcat", function(siamcat, verbose=1) {
    accessSlot(siamcat, "pred_matrix", verbose)
})

###############################################################################
#' @title Retrieve the evaluation metrics from a SIAMCAT object
#'
#' @description Function to retrieve the evaluation metrics from a SIAMCAT
#' object
#'
#' @usage eval_data(siamcat, verbose=1)
#'
#' @param siamcat (Required). A \link{siamcat-class} object
#' that contains evaluation data
#'
#' @param verbose integer, if the slot is empty, should a message be printed?
#' values can be either \code{0} (no output) or \code{1} (print message)
#'
#' @return The list of evaluation data or \code{NULL}
#'
#' @details The functions returns a list containing the evaluation metrics from
#' a SIAMCAT object. See \link{evaluate.predictions} for more information on
#' evaluation data.
#'
#' @export
#'
#' @rdname eval_data-methods
#'
#' @docType methods
#'
#' @examples
#' data(siamcat_example)
#' temp <- eval_data(siamcat_example)
#' names(temp)
#' temp$auroc
setGeneric("eval_data", function(siamcat, verbose=1)
    standardGeneric("eval_data"))
#' @rdname eval_data-methods
#' @aliases eval_data
setMethod("eval_data", "siamcat", function(siamcat, verbose=1) {
    accessSlot(siamcat, "eval_data", verbose)
})
