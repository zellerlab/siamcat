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
    make_option('--fn.mlr_models_list',           type='character',                     help='Input RData file containing the trained models'),
    make_option('--label_in',        type='character', default='NULL',     help='Input file containing labels'),
    make_option('--test_sets',       type='character', default='NULL',     help='Input file specifying which examples to use for testing'),
    make_option('--pred',            type='character', default="pred.tsv", help='Output file to which predictions will be written'),
    make_option('--model_matrix',    type='character',                     help='Input file containing information to rebuild models'),
    make_option('--model_type',      type='character', default='lasso',    help='Which plm was trained?'),
    make_option('--hyperpars',       type='character',                     help='Input file containing the hyper-parameters'))

# parse arguments
opt            <- parse_args(OptionParser(option_list=option_list))
source.dir     <- opt$srcdir
fn.test.feat   <- opt$feat_in
fn.mlr_models_list  <- opt$fn.mlr_models_list
fn.test.label  <- opt$label_in
fn.test.sample <- opt$test_sets
fn.pred        <- opt$pred
model.matrix   <- opt$model_matrix
hyperpars      <- opt$hyperpars
model.type     <- opt$model_type

# print parameters of the run
cat("=== Paramaters of the run:\n\n")
cat('source.dir         =', source.dir,        '\n')
cat('fn.test.feat       =', fn.test.feat,      '\n')
cat('fn.mlr_models_list =', fn.mlr_models_list,'\n')
cat('fn.test.label      =', fn.test.label,     '\n')
cat('fn.test.sample     =', fn.test.sample,    '\n')
cat('fn.pred            =', fn.pred,           '\n')
cat('model.matrix       =', model.matrix,      '\n')
cat('hyperpars          =', hyperpars,         '\n')
cat('model.type         =', model.type,        '\n')
cat('\n')

source.dir <- appendDirName(source.dir)
# optional parameters will be reset to NULL if specified as 'NULL', 'NONE' or 'UNKNOWN'
if (is.null(fn.test.sample) || toupper(fn.test.sample)=='NULL' || toupper(fn.test.sample)=='NONE' || toupper(fn.test.sample)=='UNKNOWN') {
  fn.test.sample = NULL
  cat('fn.test.sample not specified: applying model(s) on whole data set\n')
}

feat  <- read.features(fn.test.feat)

if (is.null(fn.test.label) || toupper(fn.test.label)=='NULL' || toupper(fn.test.label)=='NONE' || toupper(fn.test.label)=='UNKNOWN') {
  fn.test.label = NULL
  cat('fn.test.label not specified: skipping evaluation\n')
}else{

  label        <- read.labels(fn.test.label, feat)
}

# LASSO model (coefficients)
#model = NULL
#model$W = read.table(file=fn.model, sep='\t', header=TRUE, row.names=1, stringsAsFactors=FALSE, check.names=FALSE, quote='')
#num.runs = ncol(model$W)
#stopifnot(nrow(model$W) == nrow(feat))

# Read in model matrix
model.mat <- read.table(file=model.matrix, sep='\t', header = TRUE, stringsAsFactors=FALSE, row.names = 1, check.names=FALSE, quote='')

#model$W = model$W[1:dim(model$W)[1]-1,]
# parse model header
#con = file(fn.model, 'r')
#model$header <- readLines(con, 1)
#close(con)
#model$header <- parse.model.header(model$header)
start.time   <- proc.time()[1]
load(fn.mlr_models_list)
num.runs = length(models.list)

### read test data and the trained model(s)
# features
pred <- plm.predictor(feat, label, test.samples=fn.test.sample, models.list, model.mat, hyperpars, model.type)

### save prediction
### save prediction
pred.header <- paste('#Predictions for ', label$positive.label, ':', label$p.lab,
  ' [', label$header, ']', sep='')
write(pred.header, file=fn.pred, append=FALSE)
#print(pred$pred)
if (length(unique(names(pred$pred))) < length(pred$pred)) {
  suppressWarnings(write.table(pred$mat, file=fn.pred, quote=FALSE, sep='\t', row.names=TRUE, col.names=NA, append=TRUE))
} else {
  write.table(pred$pred, file=fn.pred, quote=FALSE, sep='\t', row.names=TRUE, col.names=FALSE, append=TRUE)
}
cat('\nSaved all predictions\n')

cat('\nSuccessfully applied ', model.type, ' model in ' , proc.time()[1] - start.time,
    ' seconds\n', sep='')
