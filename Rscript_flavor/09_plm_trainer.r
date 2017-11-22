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
    make_option('--feat_in',           type='character',                  help='Input file containing features'),
    make_option('--label_in',          type='character',                  help='Input file containing labels'),
    make_option('--method',                type='character', default='lasso', help='class of learner, directly passed to mlr::makeLearner'),
    make_option('--train_sets',        type='character', default='NULL',  help='Input file specifying which examples to use for training'),
    make_option('--mlr_models_list',   type='character',                  help='Output RData file to which the object with trained models will be written'),
    make_option('--stratify',          type='logical',   default=TRUE,    help='Should cross-validation for model selection be stratified 
    	                                                                        such that an approx. equal proportion of positive examples 
    	                                                                        are contained in each subset (only for binary labels)?'),
    make_option('--sel_criterion',     type='character', default='auroc', help='Evaluation criterion for model selection (options: \'acc\',
    	                                                                        \'auroc\', \'auprc\', \'auroc.2\')'),
    make_option('--min_nonzero_coeff', type='integer',   default=1,       help='Minimum number of non-zero coefficients required for a model 
    	                                                                        to be considered in model selection'),
    make_option('--model',             type='character',                  help='Text file to which the trained models will be written'),

    make_option('--model_matrix',      type='character',                  help='Output file containing information to rebuild models in 
    	                                                                        plm_predictor function')
)

 opt         <- parse_args(OptionParser(option_list=option_list))

# print parameters of the run
cat("=== 09_plm_trainer.r\n")
cat("=== Paramaters of the run:\n\n")
cat('srcdir            =', opt$srcdir,  '\n')
cat('feat_in           =', opt$feat_in, '\n')
cat('label_in          =', opt$label_in, '\n')
cat('method            =', opt$method, '\n')
cat('train_sets        =', opt$train_sets, '\n')
cat('mlr_models_list   =', opt$mlr_models_list, '\n')
cat('stratify          =', opt$stratify, '\n')
cat('sel_criterion     =', opt$sel_criterion, '\n')
cat('min_nonzero_coeff =', opt$min_nonzero_coeff, '\n')
cat('model             =', opt$model, '\n')
cat('model_matrix      =', opt$model_matrix, '\n')
cat('\n')



opt$srcdir <- appendDirName(opt$srcdir)

# optional parameters will be reset to NULL if specified as 'NULL', 'NONE' or 'UNKNOWN'
if (toupper(opt$train_sets)=='NULL' || toupper(opt$train_sets)=='NONE' || toupper(opt$train_sets)=='UNKNOWN') {
  opt$train_sets= NULL
  cat('fn.train.sample not specified: using whole data set for training\n')
}

start.time <- proc.time()[1]
set.seed(r.seed)



### read training data
# features
feat         <- read.features(opt$feat_in)
label        <- read.labels(opt$label_in, feat)

plm.out <- train.model(feat = feat, label = label,  method = opt$method, data.split=opt$train_sets, stratify = TRUE, 
                       modsel.crit  = opt$sel_criterion,  min.nonzero.coeff = opt$min_nonzero_coeff)


write.table(plm.out$out.matrix, file = opt$model_matrix, quote = FALSE, sep='\t', row.names=TRUE, col.names=NA)
models.list  <- plm.out$models.list
save(models.list, file=opt$mlr_models_list)

### save models
suppressWarnings(write.table(plm.out$W.mat, file=opt$model , quote=FALSE, sep='\t', row.names=TRUE, col.names=NA, append=FALSE))
# suppressWarnings(write.table(hyperpar.mat, file=hyper.params, quote=FALSE, sep='\t', row.names=TRUE, col.names=NA))
cat('Saved all trained models.\n')

cat('\nSuccessfully built ', plm.out$num.runs, ' LASSO models in ', proc.time()[1] - start.time,
    ' seconds\n', sep='')
