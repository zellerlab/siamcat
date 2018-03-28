#!/usr/bin/Rscript
###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
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

#based on https://github.com/joey711/phyloseq/blob/master/R/show-methods.R
#' @rdname show-methods
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
