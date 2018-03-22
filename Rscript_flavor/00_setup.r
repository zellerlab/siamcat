#!/usr/bin/Rscript
###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# EMBL Heidelberg 2012-2018
# GNU GPL 3.0
###
<<<<<<< HEAD

package.list <- c("RColorBrewer",
                  "beanplot",
                  "glmnet",
                  "LiblineaR",
                  "pROC",
                  "optparse",
                  "colorRamps",
                  "gelnet",
                  "mlr",
                  "PRROC",
                  "knitr",
                  "testthat",
                  "rmarkdown")

# script can take the location of the siamcat package as argument: Rscript 00_setup.r /path/to/SIAMCAT_0.2.0.tar.gz
# by default it is assumed to be located in './SIAMCAT_0.2.0.tar.gz'
args = commandArgs(trailingOnly = TRUE)
package.path <- if(length(args)==0) "./SIAMCAT_0.3.1.tar.gz" else args[1]

notInst      <- which(!package.list%in%installed.packages())
if(length(notInst)>0) install.packages(package.list[notInst], repos="http://cran.uni-muenster.de")

if(!"SIAMCAT"%in%installed.packages()) install.packages(package.path, repos=NULL, type="source")
=======
install.packages('session')
library(session)
restore.session(file="../dependencies/testImage.RData")
if(!"SIAMCAT"%in%installed.packages()) install.packages("./SIAMCAT_0.4.0.tar.gz", repos=NULL, type="source")
>>>>>>> development
