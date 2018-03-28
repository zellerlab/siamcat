#!/usr/bin/Rscript
###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# RScript flavor
#
# written by Georg Zeller
# with additions by Jakob Wirbel, Nicolai Karcher and Konrad Zych
# EMBL Heidelberg 2012-2017
#
# version 0.2.0
# file last updated: 30.11.2017
# GNU GPL 3.0
###

suppressMessages(library('optparse'))
suppressMessages(library('SIAMCAT'))

# define arguments
option_list    <- list(
  make_option('--metadata_in',     type='character', help='Input file containing meta-data'),
  make_option('--label_in',        type='character', help='Input file containing labels'),
  make_option('--feat_in',         type='character',               help='Input file containing features'),
  make_option('--plot',            type='character', help='Output pdf file which will contain resulting plots')
  )

# parse arguments
opt          <- parse_args(OptionParser(option_list=option_list))
# print parameters of the run
cat("=== 03_confounder_check.r\n")
cat("=== Paramaters of the run:\n\n")
cat('metadata_in  =', opt$metadata_in, '\n')
cat('label_in     =', opt$label_in, '\n')
cat('feat_in      =', opt$feat_in, '\n')
cat('plot         =', opt$plot, '\n')
cat('\n')

start.time   <- proc.time()[1]


### read label and meta- data
feat  <- read.features(opt$feat_in)
label <- read.labels(opt$label_in)
meta  <- read.meta(opt$metadata_in)
siamcat <- construct.siamcat(feat,label,meta)
### Start core function
check.confounders(siamcat, fn.plot=opt$plot)
### End core function


cat('\nSuccessfully analyzed meta-data for potential confounding in ', proc.time()[1] - start.time, ' seconds\n', sep='')
