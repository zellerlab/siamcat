#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title create a lefse input file from siamcat object
#' 
#' @description This function creates a lefse input file from siamcat object
#' 
#' @param siamcat object of class \link{siamcat-class}
#' 
#' @param filename name of the input file to which data will be save
#' 
#' @keywords siamcat.to.lefse
#' 
#' @return nothing but data is written to a file
#' 
#' @export
#' 
siamcat.to.lefse <- function(siamcat, filename="siamcat_output.txt") {
    feat   <- get.features.matrix(siamcat)
    label  <- label(siamcat)
    labelD <- label$label
    labelD[label$p.idx] <- label$p.lab
    labelD[label$n.idx] <- label$n.lab

    results <- rbind(labelD,
                        colnames(feat),
                        feat)

    rownames(results) <- c("label",
                            "sample_id",
                            rownames(feat))

    write.table(results,
                file = filename,   
                quote = FALSE, 
                sep = '\t', 
                row.names = TRUE,
                col.names = FALSE)
}

#' @title read an input file in a LEfSe input format
#' 
#' @description This reads an input file in a LEfSe input format
#'  
#' @param filename name of the input file in a LEfSe input format
#' 
#' @keywords read.lefse
#' 
#' @return a list with two elements: \itemize{
#'     \item \code{feat} a features matrix (just as returned by 
#'     \link{read.features} function
#'     \item \code{meta} a metadate matrix (just as returned by 
#'     \link{read.meta} function
#' 
#' @export
#' 
read.lefse <- function(filename="data.txt", n.meta = 1) {

    lefse <- read.csv("/Users/zych/Data/SIAMCAT datasets/lefseCRC.txt", 
                        sep = "\t", header = FALSE, stringsAsFactors = F)
    
    meta          <- lefse[1:n.meta,]
    samples.names <- lefse[n.meta+1,]
    feat          <- lefse[(n.meta+2):nrow(lefse),]
    
    rownames(meta) <- meta[,1]
    meta <- meta[,-1]
    colnames(meta) <- samples.names[,-1]
    meta <- sample_data(as.data.frame(t(meta)))
    
    rownames(feat) <- feat[,1]
    feat <- feat[,-1]
    colnames(feat) <- samples.names[,-1]
    feat <- apply(feat,c(1,2),as.numeric)
    feat <- otu_table(feat,taxa_are_rows = TRUE)
    
    return(list(feat = feat, meta = meta))
}