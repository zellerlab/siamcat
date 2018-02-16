###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# EMBL Heidelberg 2012-2017
# GNU GPL 3.0
###

#' @title Validate samples in siamcat@phyloseq@sam_datadata, labels, and features
#' @description This function checks if labels are available for all samples in features. Additionally validates siamcat@phyloseq@sam_datadata, if given
#' @param siamcat an object of class \link{siamcat}
#' @keywords SIAMCAT validate.data
#' @export
#' @return an object of class \link{siamcat}
validate.data <- function(siamcat, verbose=1){
  if(verbose>2) cat("+ starting validate.data\n")
  s.time <- proc.time()[3]
  # Check if labels are available for all samples in features
  if(verbose>2) cat("+++ checking if labels are available for all samples in features\n")
  if (length(siamcat@label@label) == dim(siamcat@phyloseq@otu_table)[2]){
    stopifnot(all(names(siamcat@label@label) %in% colnames(siamcat@phyloseq@otu_table)) && all(colnames(siamcat@phyloseq@otu_table) %in% names(siamcat@label@label)))
    # if of the same length, everything should match and be in the same order
    m <- match(names(siamcat@label@label), colnames(siamcat@phyloseq@otu_table))
    siamcat@phyloseq@otu_table <- siamcat@phyloseq@otu_table[,m]
    stopifnot(all(names(siamcat@label@label) == colnames(siamcat@phyloseq@otu_table)))

  } else if (length(siamcat@label@label) >= dim(siamcat@phyloseq@otu_table)[2]) {
    # if there are more labels than samples in features, remove them in labels
    stopifnot(all(colnames(siamcat@phyloseq@otu_table) %in% names(siamcat@label@label)))

    ids <- colnames(siamcat@phyloseq@otu_table)[colnames(siamcat@phyloseq@otu_table)%in%names(siamcat@label@label)]
    # create new labels object with reduced labels
    siamcat <- filter.label(siamcat,ids)

  } else if (length(siamcat@label@label) <= dim(siamcat@phyloseq@otu_table)[2]) {
    # if there are more samples in features, remove them and keep only the ones for which labels are available
    stopifnot(all(names(siamcat@label@label) %in% colnames(siamcat@phyloseq@otu_table)))

    cat('Warning: Removing', dim(siamcat@phyloseq@otu_table)[2]-length(siamcat@label@label),'sample(s) for which no labels are available.\n')

    m <- match(names(siamcat@label@label), colnames(siamcat@phyloseq@otu_table))
    siamcat@phyloseq@otu_table <- siamcat@phyloseq@otu_table[,m]

  }

  # Check for sample number in the different classes
  if(verbose>2) cat("+++ Checking sample number per class\n")
  for (i in siamcat@label@info$class.descr){
    if(sum(siamcat@label@label==i) <= 5) stop("Data set has only",sum(siamcat@label@label==i), "training examples of class",i," This is not enough for SIAMCAT to proceed")
    if (sum(siamcat@label@label==i) < 10){
      cat("Data set has only",sum(siamcat@label@label==i), "training examples of class",i," . Note that a dataset this small/skewed is not necessarily suitable for analysis in this pipe line." )
    }
  }

  # if siamcat@phyloseq@sam_datadata is available, check for overlap in labels
  if (!is.null(siamcat@phyloseq@sam_data)) {
    if(verbose>2) cat("+++ Check for overlap between labels and metadata\n")
    if (length(siamcat@label@label) == dim(siamcat@phyloseq@sam_data)[1]){
      stopifnot(all(names(siamcat@label@label) %in% rownames(siamcat@phyloseq@sam_data)) && all(rownames(siamcat@phyloseq@sam_data) %in% names(siamcat@label@label)))
      m    <- match(names(siamcat@label@label), rownames(siamcat@phyloseq@sam_data))
      siamcat@phyloseq@sam_data <- siamcat@phyloseq@sam_data[m,]
      stopifnot(all(names(siamcat@label@label) == rownames(siamcat@phyloseq@sam_data)))
    } else if (length(siamcat@label@label) <= dim(siamcat@phyloseq@sam_data)[1]) {
      stopifnot(all(names(siamcat@label@label) %in% rownames(siamcat@phyloseq@sam_data)))
      m <- match(names(siamcat@label@label), rownames(siamcat@phyloseq@sam_data))
      cat('Warning: Removing siamcat@phyloseq@sam_datadata information for', dim(siamcat@phyloseq@sam_data)[1]-length(siamcat@label@label),'superfluous samples.\n')
      siamcat@phyloseq@sam_data <- siamcat@phyloseq@sam_data[m,]
    } else {
      stop('Metadata is not available for all samples!')
    }
  }
siamcat@orig_feat <- otu_table(siamcat@phyloseq)
 e.time <- proc.time()[3]
if(verbose>1) cat("+ finished validate.data in",e.time-s.time,"s\n")
if(verbose==1)cat("Data succesfully validated\n")
return(siamcat)
}