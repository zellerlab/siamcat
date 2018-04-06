#!/usr/bin/Rscript
###
# SIAMCAT -  Statistical Inference of Associations between
#   Microbial Communities And host phenoTypes
# R flavor
# EMBL Heidelberg 2012-2018
# GNU GPL 3.0
###

#' Build siamcat-class objects from their components.
#' @title siamcat
#' @name siamcat
#' @description Function to construct an object of class \link{siamcat-class}
#' @param ... list of arguments needed in order to construct a SIAMCAT object
#' @return A new \link{siamcat-class} object
#' @export
#' @examples
#' # example with package data
#' fn.in.feat  <- system.file("extdata", "feat_crc_study-pop-I_N141_tax_profile_mocat_bn_specI_clusters.tsv",
#'   package = "SIAMCAT")
#' fn.in.label <- system.file("extdata", "label_crc_study-pop-I_N141_tax_profile_mocat_bn_specI_clusters.tsv",
#'   package = "SIAMCAT")
#' fn.in.meta  <- system.file("extdata", "num_metadata_crc_study-pop-I_N141_tax_profile_mocat_bn_specI_clusters.tsv",
#'   package = "SIAMCAT")
#'
#' feat  <- read.features(fn.in.feat)
#' label <- read.labels(fn.in.label)
#' meta  <- read.meta(fn.in.meta)
#' siamcat <- siamcat(feat, label, meta)
siamcat <- function(...){
  arglist   <- list(...)

  # Remove names from arglist. Will replace them based on their class
  names(arglist) <- NULL

  # ignore all but component data classes.
  component_classes <- get.component.classes("both")

  for(argNr in 1:length(arglist)){
    classOfArg <- class(arglist[[argNr]])[1]
    if(classOfArg%in%names(component_classes)){
      names(arglist)[argNr] <- component_classes[classOfArg]
    }
  }

  if(is.null(arglist$phyloseq)){
    arglistphyloseq <- arglist[sapply(names(arglist),
      is.component.class, "phyloseq")]
    arglist$phyloseq <- do.call("new", c(list(Class="phyloseq"),
      arglistphyloseq))
  }
  arglist     <- arglist[sapply(names(arglist), is.component.class, "siamcat")]
  sc          <- do.call("new", c(list(Class="siamcat"), arglist))
  return(sc)
}

# source: https://github.com/joey711/phyloseq/blob/master/R/phyloseq-class.R
#' Show the component objects classes and slot names.
#' @keywords internal
#' @return list of component classes
get.component.classes <- function(class){
  # define classes vector
  # the names of component.classes needs to be the slot names to match getSlots / splat
  component.classes.siamcat <- c("model_list", "orig_feat", "label", "norm_param", "data_split","phyloseq") #slot names
  names(component.classes.siamcat) <- c("model_list", "orig_feat", "label","norm_param", "data_split", "phyloseq") #class names

  component.classes.phyloseq <- c("otu_table", "sam_data", "phy_tree", "tax_table", "refseq") #slot names
  names(component.classes.phyloseq) <- c("otu_table", "sample_data", "phylo", "taxonomyTable", "XStringSet") #class names

  if(class=="siamcat"){
    return(component.classes.siamcat)
  }else if(class=="phyloseq"){
    return(component.classes.phyloseq)
  }else if(class=="both"){
    return(c(component.classes.siamcat,component.classes.phyloseq))
  }
}

# Returns TRUE if x is a component class, FALSE otherwise.
# This shows up over and over again in data infrastructure
#' @keywords internal
is.component.class = function(x,class){
  x%in%get.component.classes(class)
}

#' Reset features in siamcat@phylose@otu_table to those in siamcat@orig_feat
#' @title reset.features
#' @name reset.features
#' @description Function reset features in siamcat@phylose@otu_table to those in
#' siamcat@orig_feat in an object of class \link{siamcat-class}
#' @param siamcat an object of class \link{siamcat-class}t
#' @return A new \link{siamcat-class} object
#' @export
#' @examples
#'  data(siamcat_example)
#'  siamcat_example <- reset.features(siamcat_example)
reset.features <- function(siamcat){
    siamcat@phyloseq@otu_table <- siamcat@orig_feat
    return(siamcat)
}

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
get.features <- function(siamcat){
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
get.phyloseq <- function(siamcat){
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
get.features.matrix <- function(siamcat){
    return(matrix(siamcat@phyloseq@otu_table,nrow=nrow(siamcat@phyloseq@otu_table), ncol=ncol(siamcat@phyloseq@otu_table),
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
get.model_list <- function(siamcat){
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
get.model.type <- function(siamcat){
    return(model.type=siamcat@model_list@model.type)
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
get.models <- function(siamcat){
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
get.eval_data <- function(siamcat){
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
get.pred_matrix <- function(siamcat){
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
get.label.label <- function(siamcat){
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
get.label <- function(siamcat){
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
get.label.list <- function(siamcat){
    return(list(label = siamcat@label@label, header = siamcat@label@header,
                info=siamcat@label@info, positive.lab=siamcat@label@positive.lab,
                negative.lab=siamcat@label@negative.lab,n.idx=siamcat@label@n.idx,
                p.idx=siamcat@label@p.idx, n.lab=siamcat@label@n.lab,
                p.lab=siamcat@label@p.lab))
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
get.label.info <- function(siamcat){
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
filter.label <- function(siamcat,ids, verbose=1){
  labels_new      <- new("label", label = siamcat@label@label[ids],
                         header = siamcat@label@header,
                         info = siamcat@label@info,
                         positive.lab = siamcat@label@positive.lab,
                         negative.lab = siamcat@label@negative.lab,
                         n.lab = siamcat@label@n.lab,
                         p.lab = siamcat@label@p.lab)
  labels_new@n.idx <- labels_new@label == labels_new@negative.lab
  labels_new@p.idx <- labels_new@label == labels_new@positive.lab

  if (verbose > 0) cat('Keeping labels of',length(labels_new@label),'sample(s).\n')

  siamcat@label <- labels_new
  return(siamcat)
}

# based on https://github.com/joey711/phyloseq/blob/master/R/show-methods.R
#' @title Show method for siamcat class object
#' @rdname show-methods
#' @return none
#' @keywords internal
setMethod("show", "siamcat", function(object){
  cat("siamcat-class object", fill=TRUE)
  if(!is.null(object@label)) cat(paste("label()                label:           ",
  sum(object@label@n.idx), object@label@n.lab, "and",
  sum(object@label@p.idx),object@label@p.lab,"samples", sep = " "),
  fill = TRUE)
  if(length(object@norm_param)){
    cat(paste("norm_param()           norm_param:       Features normalized using",
              object@norm_param$norm.method, sep = " "), fill = TRUE)
  }
  if(length(object@data_split@num.folds)){
    cat(paste("data_split()            data_split:       ",
    object@data_split@num.resample,"cv rounds with",
              object@data_split@num.folds,"folds", sep = " "), fill = TRUE)
  }
  if(length(object@model_list@model.type)){
    cat(paste("model_list()            model_list:       ",
    length(object@model_list@models),
              object@model_list@model.type,"models", sep = " "), fill = TRUE)
  }
  if(nrow(object@pred_matrix)){
    cat(paste("pred_matrix()           pred_matrix:       Predictions for",
    nrow(object@pred_matrix),"samples from",
              ncol(object@pred_matrix),"cv rounds", sep = " "), fill = TRUE)
  }
  if(length(object@eval_data)){
    cat(paste("eval_data()             eval_data:         Average AUC:",
    round(object@eval_data$auc.average[[1]],3), sep = " "), fill = TRUE)
  }

  # print otu_table (always there).
  cat("\ncontains phyloseq-class experiment-level object @phyloseq:", fill=TRUE)
  cat(paste("phyloseq@otu_table()   OTU Table:         [ ",
    ntaxa(otu_table(object@phyloseq)), " taxa and ",
    nsamples(otu_table(object@phyloseq)), " samples ]", sep = ""), fill = TRUE)

  # print Sample Data if there
  if(!is.null(sample_data(object@phyloseq, FALSE))){
    cat(paste("phyloseq@sam_data()    Sample Data:       [ ",
    dim(sample_data(object@phyloseq))[1], " samples by ",
              dim(sample_data(object@phyloseq))[2],
              " sample variables ]", sep = ""), fill = TRUE)
  }

  # print tax Tab if there
  if(!is.null(tax_table(object@phyloseq, FALSE))){
    cat(paste("phyloseq@tax_table()   Taxonomy Table:    [ ",
    dim(tax_table(object@phyloseq))[1], " taxa by ",
              dim(tax_table(object@phyloseq))[2],
              " taxonomic ranks ]", sep = ""), fill = TRUE)
  }

  # print tree if there
  if(!is.null(phy_tree(object@phyloseq, FALSE))){
    cat(paste("phyloseq@phy_tree()    Phylogenetic Tree: [ ",
    ntaxa(phy_tree(object@phyloseq)), " tips and ",
              phy_tree(object@phyloseq)$Nnode,
              " internal nodes ]", sep = ""),
        fill = TRUE
    )
  }

  # print refseq summary if there
  if(!is.null(refseq(object@phyloseq, FALSE))){
    cat(paste("phyloseq@refseq()      ",
    class(refseq(object@phyloseq))[1], ":      [ ",
    ntaxa(refseq(object@phyloseq)), " reference sequences ]", sep = ""),
    fill=TRUE)
  }
})
