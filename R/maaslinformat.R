#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title create a MaAsLin input file from siamcat object
#' 
#' @description This function creates a MaAsLin merged PCL single input file 
#' from siamcat object
#' 
#' @param siamcat object of class \link{siamcat-class}
#' 
#' @param filename name of the input file to which data will be save
#' 
#' @keywords internal
#' 
#' @return nothing but data is written to a file
#' 
#' @examples
#' 
#' data(siamcat_example)
#' siamcat.to.maaslin(siamcat_example)
#' 
#' @export
#' 
siamcat.to.maaslin <- function(siamcat, filename="siamcat_output.pcl") {
    feat   <- get.features.matrix(siamcat)
    label  <- label(siamcat)
    meta   <- t(meta(siamcat))
    labelD <- label$label
    labelD[label$p.idx] <- label$p.lab
    labelD[label$n.idx] <- label$n.lab

    results <- rbind(labelD,
                        meta,
                        colnames(feat),
                        feat)

    rownames(results) <- c("label",
                            rownames(meta),
                            "sample_id",
                            rownames(feat))

    write.table(results,
                file = filename,   
                quote = FALSE, 
                sep = '\t', 
                row.names = TRUE,
                col.names = FALSE)
}