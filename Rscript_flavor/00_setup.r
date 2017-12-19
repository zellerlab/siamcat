#!/usr/bin/Rscript
###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# RScript flavor
#
# written by Georg Zeller
# with additions by Nicolai Karcher and Konrad Zych
# EMBL Heidelberg 2012-2017
#
# version 0.2.0
# file last updated: 09.11.2017
# GNU GPL 3.0
###

package.list <- c("RColorBrewer",
                  "beanplot",
                  "glmnet",
                  "LiblineaR",
                  "pROC",
                  "optparse",
                  "colorRamps",
                  "gelnet",
                  "mlr",
                  "PRROC")

# script can take the location of the siamcat package as argument
# by default it is assumed to be located in './SIAMCAT_0.2.0.tar.gz'
args = commandArgs(trailingOnly = TRUE)
package.path <- if(length(args)==0) "./SIAMCAT_0.2.0.tar.gz" else args[1]

notInst      <- which(!package.list%in%installed.packages())
if(length(notInst)>0) install.packages(package.list[notInst], repos="http://cran.uni-muenster.de")

if(!"SIAMCAT"%in%installed.packages()) install.packages(package.path, repos=NULL, type="source")
