###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# RScript flavor
#
# written by Georg Zeller
# with additions by Nicolai Karcher and Konrad Zych
# EMBL Heidelberg 2012-2017
#
# version 0.2.0
# file last updated: 25.06.2017
# GNU GPL 3.0
###

### parse commandline arguments
suppressMessages(library('optparse'))
suppressMessages(library('SIAMCAT'))

r.seed          <- 223311

# define arguments
option_list     <- list(
  make_option(c('-s', '--srcdir'), type='character',                 help='Source directoyr of utility scripts'),
  make_option('--label_in',        type='character',                 help='Input file containing labels'),
  make_option('--train_sets',      type='character',                 help='Output file containing training sets'),
  make_option('--test_sets',       type='character',                 help='Output file containing test sets'),
  make_option('--num_folds',       type='integer',   default=10,     help='Number of cross-validation folds (i.e. subsets, needs to be >= 2)'),
  make_option('--resample',        type='integer',   default=0,      help='Resampling rounds (values <= 1 deactivate resampling)'),
  make_option('--stratify',        type='logical',   default=TRUE,   help='Should cross-validation be stratified such that an approx. 
                                                                           equal proportion of positive examples are contained in each subset 
                                                                           (only for binary labels)?'),
  make_option('--inseparable',     type='character', default='NULL', help=''),
  make_option('--metadata_in',     type='character',                 help='Input file containing metadata (only required if argument 
                                                                           \'inseparable\' is specified)')
)

# parse arguments
opt            <- parse_args(OptionParser(option_list=option_list))

source.dir     <- opt$srcdir
fn.in.label    <- opt$label_in
fn.train.folds <- opt$train_sets
fn.test.folds  <- opt$test_sets
num.folds      <- opt$num_folds
num.resample   <- opt$resample
stratify       <- opt$stratify
inseparable    <- opt$inseparable
fn.in.meta     <- opt$metadata_in

# print parameters of the run
cat("=== Paramaters of the run:\n\n")
cat('source.dir     =', source.dir, '\n')
cat('fn.in.label    =', fn.in.label, '\n')
cat('fn.train.folds =', fn.train.folds, '\n')
cat('fn.test.folds  =', fn.test.folds, '\n')
cat('num.folds      =', num.folds, '\n')
cat('num.resample   =', num.resample, '\n')
cat('stratify       =', stratify, '\n')
cat('inseparable    =', inseparable, '\n')
cat('fn.in.meta     =', fn.in.meta, '\n')
cat('\n')
if (substr(source.dir, nchar(source.dir), nchar(source.dir)) != '/') {
  source.dir = paste(source.dir, '/', sep='')
}

start.time  <- proc.time()[1]
set.seed(r.seed)


label         <- read.labels(fn.in.label)

### Core function sourced from the library
training.data <- data.splitter(label = label, num.folds, num.resample, stratify, inseparable, fn.in.meta)


write('#Cross-validation training folds', file=fn.train.folds, append=FALSE)
write('#Cross-validation test folds',     file=fn.test.folds, append=FALSE)

for (r in 1:training.data$num.resample) {
  for (f in 1:training.data$num.folds) {
    # append training and test examples, one line per fold
    fold.name = paste('>cv_fold', ifelse(num.resample>1, paste(f, '_rep', r, sep=''), as.character(f)), sep='')
    write.table(t(c(fold.name, training.data$training.folds[[r]][[f]])), file=fn.train.folds, quote=FALSE, 
                sep='\t', row.names=FALSE, col.names=FALSE, append=TRUE)
    write.table(t(c(fold.name, training.data$test.folds[[r]][[f]])),     file=fn.test.folds,  quote=FALSE, 
                sep='\t', row.names=FALSE, col.names=FALSE, append=TRUE)
  }
}
if (stratify) {
  cat('\nSuccessfully created data split for ', num.folds, '-fold stratified cross-validation', sep='')
} else {
  cat('Successfully created data split for ', num.folds, '-fold cross-validation', sep='')
}
if (num.resample > 1) {
  cat(' with ', num.resample, ' times repeated resampling\n', sep='')
} else {
  cat('\n')
}