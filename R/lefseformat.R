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
#' @keywords internal
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