#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0

#' @title Summarize features
#' 
#' @description Summarize features
#'     
#' @param siamcat object of class \link{siamcat-class}
#' 
#' @param siamcat at which level to summarize (g__ = genus)
#' 
#' @param verbose control output: \code{0} for no output at all, \code{1}
#'     for only information about progress and success, \code{2} for normal
#'     level of information and \code{3} for full debug information,
#'     defaults to \code{1}
#'     
#' @keywords internal
#' 
#' @export
#' 
#' @return object of class \link{siamcat-class} with the
#'     slot \code{eval_data} filled
#'     
#' @examples
#'
#'     data(siamcat_example)
#'     # simple working example
#'     siamcat_evaluated <- summarize.features(siamcat_example)
#'
summarize.features <- function(siamcat, level = "g__", 
  keep_original_names = TRUE, verbose=1){
  
  if (verbose > 1)
    message("+ starting summarize.features")
  
  s.time <- proc.time()[3]
  
  if (verbose > 2) message("+++ summarizing on a level: ",level,"\n")
  feat <- get.features.matrix(siamcat)
  rowSplitted <- strsplit(rownames(feat),"\\.")
  genera    <- NULL
  origNames   <- NULL
  
  for(row in rowSplitted){
    levelPos <- grep(level,row)
    if(length(levelPos)>0){
      genera <- c(genera,row[levelPos])
      origNames<- c(origNames,paste(row[1:levelPos],collapse="."))
    }
  }
  
  generaSplitted <- strsplit(genera,level)
  generaSpVec    <- do.call(rbind,generaSplitted)[,2]
  generaSpVecMap <- cbind(generaSpVec,origNames)
  generaSpVecMap <- generaSpVecMap[!duplicated(generaSpVecMap[,1]),]
  summarized     <- NULL
  
  for(genus in generaSpVecMap[,1]){
    
    curRows <- feat[grep(paste0(level,genus),rownames(feat)),]
    
    if(!is.null(dim(curRows))){
      
      summarized <- rbind(summarized,colSums(curRows))
      
    }else{
      
      summarized <- rbind(summarized,curRows)
      
    }
    
  }
  
  if (verbose > 2) message("+++ summarized features table contains: ",
    nrow(summarized)," features\n")
  
  rownames(summarized) <- generaSpVecMap[,1]
  if(any(rownames(summarized)=="")){
    summarized           <- summarized[-which(rownames(summarized)==""),]
  }

  if(keep_original_names){
      rownames(summarized) <-  generaSpVecMap[,2]
  }

  colnames(summarized) <- colnames(feat)
  
  
  e.time <- proc.time()[3]
  if (verbose > 1)
    message(paste(
      "+ finished summarize.features in",
      formatC(e.time - s.time, digits = 3),
      "s"
    ))
  if (verbose == 1)
    message("Summarized features successfully.")
  features(siamcat) <- otu_table(summarized,taxa_are_rows = TRUE)
  return(siamcat)
}
