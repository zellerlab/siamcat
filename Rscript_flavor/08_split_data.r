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

### parse commandline arguments
suppressMessages(library('optparse'))
suppressMessages(library('SIAMCAT'))

r.seed          <- 223311

# define arguments
option_list     <- list(
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
# print parameters of the run
cat("=== 08_data_splitter.r\n")
cat("=== Paramaters of the run:\n\n")
cat('label_in    =', opt$label_in, '\n')
cat('train_sets  =', opt$train_sets, '\n')
cat('test_sets   =', opt$test_sets, '\n')
cat('num_folds   =', opt$num_folds, '\n')
cat('resample    =', opt$resample, '\n')
cat('stratify    =', opt$stratify, '\n')
cat('inseparable =', opt$inseparable, '\n')
cat('metadata_in =', opt$metadata_in, '\n')
cat('\n')

start.time  <- proc.time()[1]
set.seed(r.seed)


label         <- read.labels(opt$label_in)

### Core function sourced from the library
training.data <- data.splitter(label = label,
                               num.folds=opt$num_folds,
                               num.resample=opt$resample, 
                               stratify=opt$stratify, 
                               inseparable=opt$inseparable,
                               meta=opt$metadata_in)

write('#Cross-validation training folds', file=opt$train_sets, append=FALSE)
write(paste0('#num.folds:\n#',opt$num_folds), file=opt$train_sets, append=TRUE)
write('#Cross-validation test folds',     file=opt$test_sets, append=FALSE)

for (r in 1:training.data$num.resample) {
  for (f in 1:training.data$num.folds) {
    # append training and test examples, one line per fold
    fold.name = paste('>cv_fold', ifelse(opt$resample>1, paste(f, '_rep', r, sep=''), as.character(f)), sep='')
    write.table(t(c(fold.name, training.data$training.folds[[r]][[f]])), file=opt$train_sets, quote=FALSE, 
                sep='\t', row.names=FALSE, col.names=FALSE, append=TRUE)
    write.table(t(c(fold.name, training.data$test.folds[[r]][[f]])),     file=opt$test_sets,  quote=FALSE, 
                sep='\t', row.names=FALSE, col.names=FALSE, append=TRUE)
  }
}
if (opt$stratify) {
  cat('\nSuccessfully created data split for ', opt$num_folds, '-fold stratified cross-validation', sep='')
} else {
  cat('Successfully created data split for ', opt$num_folds, '-fold cross-validation', sep='')
}
if (opt$resample > 1) {
  cat(' with ', opt$resample, ' times repeated resampling\n', sep='')
} else {
  cat('\n')
}