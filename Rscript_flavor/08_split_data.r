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
  make_option('--feat_in',         type='character',                 help='Input file containing features'),
  make_option('--data_split',      type='character',                 help='Output file containing data_split object'),
  make_option('--num_folds',       type='integer',   default=10,     help='Number of cross-validation folds (i.e. subsets, needs to be >= 2)'),
  make_option('--resample',        type='integer',   default=0,      help='Resampling rounds (values <= 1 deactivate resampling)'),
  make_option('--stratify',        type='logical',   default=TRUE,   help='Should cross-validation be stratified such that an approx. 
                                                                           equal proportion of positive examples are contained in each subset 
                                                                           (only for binary labels)?'),
  make_option('--inseparable',     type='character', default=NULL,   help=''),
  make_option('--metadata_in',     type='character',                 help='Input file containing metadata (only required if argument 
                                                                           \'inseparable\' is specified)'),
  make_option('--subdivide_train_set', type='logical',   default=FALSE,  help='The train data will be subdevided into different files by folds and resamplings.
                                                                              This enables to parallelization of the model training.')
)

# parse arguments
opt            <- parse_args(OptionParser(option_list=option_list))
# print parameters of the run
cat("=== 08_data_splitter.r\n")
cat("=== Paramaters of the run:\n\n")
cat('label_in    =', opt$label_in, '\n')
cat('feat_in     =', opt$feat_in, '\n')
cat('data_split  =', opt$data_split, '\n')
cat('num_folds   =', opt$num_folds, '\n')
cat('resample    =', opt$resample, '\n')
cat('stratify    =', opt$stratify, '\n')
cat('inseparable =', opt$inseparable, '\n')
cat('metadata_in =', opt$metadata_in, '\n')
cat('subdivide_train_set =', opt$subdivide_train_set, '\n')
cat('\n')

start.time  <- proc.time()[1]
set.seed(r.seed)


feat  <- read.features(opt$feat_in)
label <- read.labels(opt$label_in, feat)
if(!is.null(opt$inseparable)){
  meta    <- read.meta(opt$metadata_in)
  siamcat <- siamcat(feat,label,meta)
}else{
  siamcat <- siamcat(feat,label)
}
### Core function sourced from the library
siamcat <- data.splitter(siamcat,
                               num.folds=opt$num_folds,
                               num.resample=opt$resample, 
                               stratify=opt$stratify, 
                               inseparable=opt$inseparable)
# write headers:

data_split <- siamcat@data_split
save(data_split,file=opt$data_split)

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