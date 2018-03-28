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
# file last updated: 26.06.2017
# GNU GPL 3.0
###

### parse commandline arguments
suppressMessages(library('optparse'))
suppressMessages(library('SIAMCAT'))
suppressMessages(library('methods'))

# define arguments
  option_list = list(
    make_option('--feat_in',         type='character',                     help='Input file containing features'),
    make_option('--mlr_models_list', type='character',                     help='Input RData file containing the trained models'),
    make_option('--label_in',        type='character',                     help='Input file containing labels'),
    make_option('--data_split',      type='character',                     help='Input file containing data_split object'),
    make_option('--pred',            type='character', default="pred.tsv", help='Output file to which predictions will be written')
)
# parse arguments
opt            <- parse_args(OptionParser(option_list=option_list))
# print parameters of the run
cat("=== 10_plm_predictor.r\n")
cat("=== Paramaters of the run:\n\n")
cat('feat_in         =', opt$feat_in,      '\n')
cat('mlr_models_list =', opt$mlr_models_list,'\n')
cat('label_in        =', opt$label_in,     '\n')
cat('data_split        =', opt$data_split, '\n')
cat('pred            =', opt$pred,           '\n')
cat('\n')


feat       <- read.features(opt$feat_in)
label      <- read.labels(opt$label_in)
siamcat    <- construct.siamcat(feat,label)

start.time   <- proc.time()[1]
load(opt$mlr_models_list)
siamcat@model_list <- model_list

load(opt$data_split)
siamcat@data_split <- data_split


siamcat <- make.predictions(siamcat)


### save prediction
pred.header <- paste('#Predictions for ', label@positive.lab, ':', label@p.lab,
  ' [', label@header, ']', sep='')
write(pred.header, file=opt$pred, append=FALSE)
#print(pred$pred)
suppressWarnings(write.table(siamcat@pred_matrix, file=opt$pred, quote=FALSE, sep='\t', row.names=TRUE, col.names=NA, append=TRUE))
cat('\nSaved all predictions\n')

cat('\nSuccessfully made preictions with the model in ' , proc.time()[1] - start.time,
    ' seconds\n', sep='')
