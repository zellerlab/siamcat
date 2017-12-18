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
    make_option(c('-s', '--srcdir'), type='character',                     help='Source directory of this and other utility scripts'),
    make_option('--feat_in',         type='character',                     help='Input file containing features'),
    make_option('--mlr_models_list', type='character',                     help='Input RData file containing the trained models'),
    make_option('--label_in',        type='character', default='NULL',     help='Input file containing labels'),
    make_option('--test_sets',       type='character', default='NULL',     help='Input file specifying which examples to use for testing'),
    make_option('--pred',            type='character', default="pred.tsv", help='Output file to which predictions will be written'),
    make_option('--model_matrix',    type='character',                     help='Input file containing information to rebuild models'))

# parse arguments
opt            <- parse_args(OptionParser(option_list=option_list))
# print parameters of the run
cat("=== 10_plm_predictor.r\n")
cat("=== Paramaters of the run:\n\n")
cat('srcdir          =', opt$srcdir,        '\n')
cat('feat_in         =', opt$feat_in,      '\n')
cat('mlr_models_list =', opt$fn.mlr_models_list,'\n')
cat('label_in        =', opt$label_in,     '\n')
cat('test_sets       =', opt$test_sets,    '\n')
cat('pred            =', opt$pred,           '\n')
cat('model_matrix    =', opt$model_matrix,      '\n')
cat('\n')

source.dir <- appendDirName(opt$srcdir)
# optional parameters will be reset to NULL if specified as 'NULL', 'NONE' or 'UNKNOWN'
if (is.null(opt$test_sets)) {
  opt$test_sets = NULL
  cat('fn.test.sample not specified: applying model(s) on whole data set\n')
}

feat  <- read.features(opt$feat_in)

if (is.null(opt$label_in)) {
  cat('fn.test.label not specified: skipping evaluation\n')
}else{
  label      <- read.labels(opt$label_in, feat)
}

model.mat <- read.table(file=opt$model_matrix, sep='\t', header = TRUE, stringsAsFactors=FALSE, row.names = 1, check.names=FALSE, quote='')


start.time   <- proc.time()[1]
load(opt$mlr_models_list)

pred <- make.predictions(feat=feat,
                         label=label,
                         data.split=opt$test_sets,
                         models.list=plm.out$models.list,
                         model.mat=model.mat)


### save prediction
pred.header <- paste('#Predictions for ', label$positive.label, ':', label$p.lab,
  ' [', label$header, ']', sep='')
write(pred.header, file=opt$pred, append=FALSE)
#print(pred$pred)
if (length(unique(names(pred$pred))) < length(pred$pred)) {
  suppressWarnings(write.table(pred$mat, file=opt$pred, quote=FALSE, sep='\t', row.names=TRUE, col.names=NA, append=TRUE))
} else {
  write.table(pred$pred, file=opt$pred, quote=FALSE, sep='\t', row.names=TRUE, col.names=FALSE, append=TRUE)
}
cat('\nSaved all predictions\n')

cat('\nSuccessfully made preictions with the model in ' , proc.time()[1] - start.time,
    ' seconds\n', sep='')
