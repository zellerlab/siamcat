#!/usr/bin/Rscript
###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# EMBL Heidelberg 2012-2018
# GNU GPL 3.0
###

dependecies.list <- list.files("../dependencies/",pattern = "*.gz")
sapply(dependecies.list, install.packages,repos=NULL, type="source")
