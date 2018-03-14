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

# define arguments
  option_list = list(
  #make_option('--pkgdir', type='character', help='Source directory of dataprep'),
  make_option('--label_in',              type='character',                help='Input file containing labels'),
  make_option('--pred',               type='character',                help='Input file containing the trained classification model(s)'),
  make_option('--plot',               type='character',                help='Output file for plotting'),
  make_option('--write_eval_results', type='logical',   default=FALSE, help='Should calculated parameters be written into tab-delimited file? (Necessary for generation of test files'),
  make_option('--output_results',     type='character', default="eval_results.tsv",  help='Output file containing evaluation results (only necessary when write_eval_results is set to TRUE')
  )

# parse arguments
opt            <- parse_args(OptionParser(option_list=option_list))
cat("=== 11_model_evaler.r\n")
cat("=== Paramaters of the run:\n\n")
cat('label_in           =', opt$label_in, '\n')
cat('pred               =', opt$pred, '\n')
cat('plot               =', opt$plot, '\n')
cat('output_results     =', opt$output_results, '\n')
cat('\n')

### If variable source.dir does not end with "/", append "/" to end of source.dir
start.time <- proc.time()[1]
label      <- read.labels(opt$label_in)
siamcat <- siamcat(label)

pred <- read.table(file=opt$pred, sep='\t', header=TRUE, row.names=1, check.names=FALSE, comment.char="#")
pred <- as.matrix(pred)
siamcat@predMatrix <- pred

siamcat <-  evaluate.predictions(siamcat)

evaluation.model.plot(siamcat, fn.plot=opt$plot)


# The following code writes calculated auroc and aupr- values into a file for testing.
if (opt$write_eval_results == TRUE){
  # Testing only makes sense if dim(pred)[2] > 1
  if(ncol(pred) == 1){
    if(!is.null(siamcat@evalData$auc.all)) write.table(t(aucs), file=opt$output_results, quote=FALSE, sep='\t', col.names=FALSE, append=FALSE, row.names="auroc values")
    write.table(t(siamcat@evalData$auc.average), file=opt$output_results, quote=FALSE, sep='\t', col.names=FALSE, append=TRUE, row.names="auprc values")
  }else{
    cat("Only one prediction available, ignoring the write_eval_results option.\n")
  }

}

cat('\nSuccessfully evaluated predictions in ', proc.time()[1] - start.time,
    ' seconds\n', sep='')
