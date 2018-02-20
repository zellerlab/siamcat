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
suppressMessages(library('session'))
restore.session(file="../dependencies/testImage.RData")
suppressMessages(library('SIAMCAT'))

# parameters that cannot be configured externally # TODO
C.vec        <- 10^seq(-2,3,length=6+5+9)   # values for C (regularization strength) to be tested in model selection
r.seed       <- 2133                       # initiliazing the pseudo-random number generator to get reproducible results
DEBUG.CHECKS <- FALSE                # performs additional checks (asserting that custom log-reg prediction code is correct)


# define arguments
 option_list <- list(
    make_option('--feat_in',           type='character',                  help='Input file containing features'),
    make_option('--label_in',          type='character',                  help='Input file containing labels'),
    make_option('--method',            type='character', default='lasso', help='class of learner, directly passed to mlr::makeLearner'),
    make_option('--train_sets',        type='character', default='NULL',  help='Input file specifying which examples to use for training'),
    make_option('--mlr_models_list',   type='character',                  help='Output RData file to which the object with trained models will be written'),
    make_option('--stratify',          type='logical',   default=TRUE,    help='Should cross-validation for model selection be stratified
    	                                                                        such that an approx. equal proportion of positive examples
    	                                                                        are contained in each subset (only for binary labels)?'),
    make_option('--sel_criterion',     type='character', default='auc',   help='Evaluation criterion for model selection (options: \'acc\',
    	                                                                        \'auc\', \'auprc\', \'f1\')'),
    make_option('--min_nonzero_coeff', type='integer',   default=1,       help='Minimum number of non-zero coefficients required for a model
    	                                                                        to be considered in model selection'),
    make_option('--param_set',          type='character', default=NULL,    help='a list of extra parameters for mlr run, may contain: cost - for lasso_ll and ridge_ll;
                                                                                 alpha for enet and ntree, mtry for RandomForrest')
)

 opt         <- parse_args(OptionParser(option_list=option_list))

# print parameters of the run
cat("=== 09_plm_trainer.r\n")
cat("=== Paramaters of the run:\n\n")
cat('feat_in           =', opt$feat_in, '\n')
cat('label_in          =', opt$label_in, '\n')
cat('method            =', opt$method, '\n')
cat('train_sets        =', opt$train_sets, '\n')
cat('mlr_models_list   =', opt$mlr_models_list, '\n')
cat('stratify          =', opt$stratify, '\n')
cat('sel_criterion     =', opt$sel_criterion, '\n')
cat('min_nonzero_coeff =', opt$min_nonzero_coeff, '\n')
cat('param_set         =', opt$param_set, '\n')
cat('\n')


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
siamcat      <- siamcat(feat,label)
num.runs      <- 1
con           <- file(opt$train_sets, 'r')
input         <- readLines(con)
close(con)
num.folds    <- 2
fold.name = list()
fold.exm.idx = list()
num.folds     <- as.numeric(strsplit(input[[3]],"#")[[1]][2])
for (i in 1:length(input)) {
        l               <- input[[i]]
        if (substr(l, 1, 1) != '#') {
          num.runs                 <- num.runs + 1
          s                        <- unlist(strsplit(l, '\t'))
          fold.name[[num.runs]]    <- substr(s[1], 2, nchar(s[1]))
          ### Note that the %in%-operation is order-dependend.
          fold.exm.idx[[num.runs]] <- which(names(label@label) %in% as.vector(s[2:length(s)]))
          cat(fold.name[[num.runs]], 'contains', length(fold.exm.idx[[num.runs]]), paste0(mode, 'ing'), 'examples\n')
        }
}
siamcat@dataSplit <- list(fold.name = fold.name,fold.exm.idx = fold.exm.idx, num.runs = num.runs, num.folds = num.folds)

siamcat  <- train.model(siamcat
                       method = opt$method,
                       data.split=opt$train_sets,
                       stratify = opt$stratify,
                       modsel.crit  = opt$sel_criterion,
                       min.nonzero.coeff = opt$min_nonzero_coeff,
                       param.set = opt$param_set)

modelsList <- siamcat@modelList
save(modelsList , file=opt$mlr_models_list)
cat('\n++++++++++++++++++++\nSuccessfully trained models in ', proc.time()[1] - start.time,
    ' seconds\n++++++++++++++++++++\n', sep='')
