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
suppressMessages(library('methods'))
suppressMessages(library('SIAMCAT'))

# parameters that cannot be configured externally # TODO
C.vec        <- 10^seq(-2,3,length=6+5+9)   # values for C (regularization strength) to be tested in model selection
r.seed       <- 2133                       # initiliazing the pseudo-random number generator to get reproducible results
DEBUG.CHECKS <- FALSE                # performs additional checks (asserting that custom log-reg prediction code is correct)


# define arguments
 option_list <- list(
    make_option(c('-s', '--srcdir'),   type='character',                  help='Source directory of this and other utility scripts'),
    make_option('--label_in',          type='character',                  help='Input file containing labels'),
    make_option('--feat_in',           type='character',                  help='Input file containing features'),
    make_option('--train_sets',        type='character', default='NULL',  help='Input file specifying which examples to use for training'),
    make_option('--model',             type='character',                  help='Text to which the trained models will be written'),
    make_option('--mlr_models_list',   type='character',                  help='Output RData file to which the object with trained models will be written'),
    make_option('--num_folds',         type='integer',   default=5,       help='Number of cross-validation folds for model selection 
    	                                                                        (i.e. subsets, needs to be >= 2)'),
    make_option('--stratify',          type='logical',   default=TRUE,    help='Should cross-validation for model selection be stratified 
    	                                                                        such that an approx. equal proportion of positive examples 
    	                                                                        are contained in each subset (only for binary labels)?'),
    make_option('--sel_criterion',     type='character', default='auroc', help='Evaluation criterion for model selection (options: \'acc\',
    	                                                                        \'auroc\', \'auprc\', \'auroc.2\')'),
    make_option('--min_nonzero_coeff', type='integer',   default=1,       help='Minimum number of non-zero coefficients required for a model 
    	                                                                        to be considered in model selection'),
    make_option('--model_matrix',      type='character',                  help='Output file containing information to rebuild models in 
    	                                                                        plm_predictor.r'),
    make_option('--model_type',        type='character', default='lasso', help='Which plm should be trained?'),
    make_option('--hyperpars',         type='character',                  help='Output file containing the hyper-parameters')
)

# parse arguments
opt               <- parse_args(OptionParser(option_list=option_list))
source.dir        <- opt$srcdir
fn.train.label    <- opt$label_in
fn.train.feat     <- opt$feat_in
fn.model          <- opt$model
fn.mlr_models_list<- opt$mlr_models_list
fn.train.sample   <- opt$train_sets
num.folds         <- opt$num_folds
stratify          <- opt$stratify
modsel.crit       <- opt$sel_criterion
min.nonzero.coeff <- opt$min_nonzero_coeff
model.matrix      <- opt$model_matrix
model.type        <- opt$model_type
hyper.params      <- opt$hyperpars

# print parameters of the run
cat("=== Paramaters of the run:\n\n")
cat('source.dir        =', source.dir, '\n')
cat('fn.train.label    =', fn.train.label, '\n')
cat('fn.train.feat     =', fn.train.feat, '\n')
cat('fn.train.sample   =', fn.train.sample, '\n')
cat('fn.model          =', fn.model, '\n')
cat('fn.mlr_models_list=', fn.mlr_models_list, '\n')
cat('num.folds         =', num.folds, '\n')
cat('stratify          =', stratify, '\n')
cat('modsel.crit       =', modsel.crit, '\n')
cat('min.nonzero.coeff =', min.nonzero.coeff, '\n')
cat('model.type        =', model.type, '\n')
cat('model.matrix      =', model.matrix, '\n')
cat('hyper.params      =', hyper.params, '\n')
cat('\n')

source.dir <- appendDirName(source.dir)

# optional parameters will be reset to NULL if specified as 'NULL', 'NONE' or 'UNKNOWN'
if (toupper(fn.train.sample)=='NULL' || toupper(fn.train.sample)=='NONE' || toupper(fn.train.sample)=='UNKNOWN') {
  fn.train.sample = NULL
  cat('fn.train.sample not specified: using whole data set for training\n')
}

start.time <- proc.time()[1]
set.seed(r.seed)



### read training data
# features
feat         <- read.features(fn.train.feat)
label        <- read.labels(fn.train.label, feat)

plm.out <- plm.trainer(feat = feat, label = label, data.split = fn.train.sample, stratify = stratify, modsel.crit = modsel.crit, 
                       min.nonzero.coeff=min.nonzero.coeff, model.type = model.type)


runParams <- list(num.folds=num.folds, stratify=stratify, modsel.crit=modsel.crit, min.nonzero.coeff=min.nonzero.coeff)
save(runParams,file="runParams.RData")
write.table(plm.out$out.matrix,   file=model.matrix, quote = FALSE, sep='\t', row.names=TRUE, col.names=NA)
write.table(plm.out$hyperpar.mat, file=hyper.params, quote = FALSE, sep='\t', row.names=TRUE, col.names=NA)
models.list  <- plm.out$models.list
print(models.list)
save(models.list, file=fn.mlr_models_list)

### save models
write(plm.out$model.header, file=fn.model, append=FALSE)
suppressWarnings(write.table(plm.out$W.mat, file=fn.model, quote=FALSE, sep='\t', row.names=TRUE, col.names=NA, append=TRUE))
# suppressWarnings(write.table(hyperpar.mat, file=hyper.params, quote=FALSE, sep='\t', row.names=TRUE, col.names=NA))
cat('Saved all trained models.\n')

cat('\nSuccessfully built ', plm.out$num.runs, ' LASSO models in ', proc.time()[1] - start.time,
    ' seconds\n', sep='')
