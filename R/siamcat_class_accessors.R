#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between Microbial
### Communities And host phenoTypes
### EMBL Heidelberg 2012-2018 GNU GPL 3.0

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
  if( is.component.class(siamcat) ){
    # If physeq is a component class, might return as-is. Depends on slot.
    if( inherits(siamcat, get.component.classes()[slot]) ){
      # if slot-name matches, return physeq as-is.
      out = siamcat
    } else {
      # If slot/component mismatch, set out to NULL. Test later if this is an error.			
      out = NULL
    }
  } else if(!slot %in% slotNames(siamcat) ){
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
#' Retrieve model_list from object.
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
#'  that contains a label or instance of \link{data_split-class}.
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
#'  that contains a label or instance of \link{data_split-class}.
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
# Return as-is if already a otu_table object
#' @rdname physeq-methods
setMethod("physeq", "phyloseq", function(siamcat){ return(siamcat) })

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

#' Access features in siamcat@phylose@otu_table
#' @title get.features.matrix
#' @name get.features.matrix
#' @description Function to access features in siamcat@phylose@otu_table
#' @param siamcat an object of class \link{siamcat-class}t
#' @return Features as a matrix
#' @export
#' @examples
#'  data(siamcat_example)
#'  feat <- get.features.matrix(siamcat_example)
get.features.matrix <- function(siamcat) {
    return(matrix(siamcat@phyloseq@otu_table, nrow = nrow(siamcat@phyloseq@otu_table), ncol = ncol(siamcat@phyloseq@otu_table), 
        dimnames = list(rownames(siamcat@phyloseq@otu_table), colnames(siamcat@phyloseq@otu_table))))
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

#' @title Filter samples from \code{siamcat@label}
#' @description This functions filters \code{siamcat@label}.
#' @param siamcat an object of class \link{siamcat-class}
#' @param ids names of samples to be left in the \code{siamcat@label}
#' @param verbose control output: \code{0} for no output at all, \code{1} for more
#'  information about progress and success, defaults to \code{1}
#' @keywords filter.label
#' @export
#' @return siamcat an object of class \link{siamcat-class}
#' @examples
#'
#'  data(siamcat_example)
#'  # simple working example
#'  siamcat_filtered <- filter.label(siamcat_example, ids=c(1:10))
#'
filter.label <- function(siamcat, ids, verbose = 1) {
    labels_new <- new("label", label = siamcat@label@label[ids], header = siamcat@label@header, info = siamcat@label@info, 
        positive.lab = siamcat@label@positive.lab, negative.lab = siamcat@label@negative.lab, n.lab = siamcat@label@n.lab, 
        p.lab = siamcat@label@p.lab)
    labels_new@n.idx <- labels_new@label == labels_new@negative.lab
    labels_new@p.idx <- labels_new@label == labels_new@positive.lab
    
    if (verbose > 0) 
        cat("Keeping labels of", length(labels_new@label), "sample(s).\n")
    
    siamcat@label <- labels_new
    return(siamcat)
}

# based on https://github.com/joey711/phyloseq/blob/master/R/show-methods.R
#' @title Show method for siamcat class object
#' @rdname show-methods
#' @return none
#' @keywords internal
setMethod("show", "siamcat", function(object) {
    cat("siamcat-class object", fill = TRUE)
    if (!is.null(object@label)) 
        cat(paste("label()                label:           ", sum(object@label@n.idx), object@label@n.lab, "and", sum(object@label@p.idx), 
            object@label@p.lab, "samples", sep = " "), fill = TRUE)
    if (length(object@norm_param)) {
        cat(paste("norm_param()           norm_param:       Features normalized using", object@norm_param$norm.method, 
            sep = " "), fill = TRUE)
    }
    if (length(object@data_split@num.folds)) {
        cat(paste("data_split()            data_split:       ", object@data_split@num.resample, "cv rounds with", object@data_split@num.folds, 
            "folds", sep = " "), fill = TRUE)
    }
    if (length(object@model_list@model.type)) {
        cat(paste("model_list()            model_list:       ", length(object@model_list@models), object@model_list@model.type, 
            "models", sep = " "), fill = TRUE)
    }
    if (nrow(object@pred_matrix)) {
        cat(paste("pred_matrix()           pred_matrix:       Predictions for", nrow(object@pred_matrix), "samples from", 
            ncol(object@pred_matrix), "cv rounds", sep = " "), fill = TRUE)
    }
    if (length(object@eval_data)) {
        cat(paste("eval_data()             eval_data:         Average AUC:", round(object@eval_data$auc.average[[1]], 
            3), sep = " "), fill = TRUE)
    }
    
    # print otu_table (always there).
    cat("\ncontains phyloseq-class experiment-level object @phyloseq:", fill = TRUE)
    cat(paste("phyloseq@otu_table()   OTU Table:         [ ", ntaxa(otu_table(object@phyloseq)), " taxa and ", nsamples(otu_table(object@phyloseq)), 
        " samples ]", sep = ""), fill = TRUE)
    
    # print Sample Data if there
    if (!is.null(sample_data(object@phyloseq, FALSE))) {
        cat(paste("phyloseq@sam_data()    Sample Data:       [ ", dim(sample_data(object@phyloseq))[1], " samples by ", 
            dim(sample_data(object@phyloseq))[2], " sample variables ]", sep = ""), fill = TRUE)
    }
    
    # print tax Tab if there
    if (!is.null(tax_table(object@phyloseq, FALSE))) {
        cat(paste("phyloseq@tax_table()   Taxonomy Table:    [ ", dim(tax_table(object@phyloseq))[1], " taxa by ", 
            dim(tax_table(object@phyloseq))[2], " taxonomic ranks ]", sep = ""), fill = TRUE)
    }
    
    # print tree if there
    if (!is.null(phy_tree(object@phyloseq, FALSE))) {
        cat(paste("phyloseq@phy_tree()    Phylogenetic Tree: [ ", ntaxa(phy_tree(object@phyloseq)), " tips and ", phy_tree(object@phyloseq)$Nnode, 
            " internal nodes ]", sep = ""), fill = TRUE)
    }
    
    # print refseq summary if there
    if (!is.null(refseq(object@phyloseq, FALSE))) {
        cat(paste("phyloseq@refseq()      ", class(refseq(object@phyloseq))[1], ":      [ ", ntaxa(refseq(object@phyloseq)), 
            " reference sequences ]", sep = ""), fill = TRUE)
    }
})
