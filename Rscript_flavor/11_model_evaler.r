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
  make_option(c('-s', '--srcdir'),    type='character',                help='Source directory of this and other utility scripts'),
  make_option('--label',              type='character',                help='Input file containing labels'),
  make_option('--pred',               type='character',                help='Input file containing the trained classification model(s)'),
  make_option('--plot',               type='character',                help='Output file for plotting'),
  make_option('--write_eval_results', type='logical',   default=FALSE, help='Should calculated parameters be written into tab-delimited file? (Necessary for generation of test files'),
  make_option('--output_results',     type='character', default="eval_results.tsv",  help='Output file containing evaluation results (only necessary when write_eval_results is set to TRUE'),
  make_option('--model_type',         type='character',                help='Model type')
  )

# parse arguments
opt            <- parse_args(OptionParser(option_list=option_list))
source.dir     <- opt$srcdir
fn.label       <- opt$label
fn.pred        <- opt$pred
fn.plot        <- opt$plot
write.out      <- opt$write_eval_results
output.results <- opt$output_results
model.type     <- opt$model_type

cat("=== Paramaters of the run:\n\n")
cat('source.dir     =', source.dir, '\n')
cat('fn.label       =', fn.label, '\n')
cat('fn.pred        =', fn.pred, '\n')
cat('fn.plot        =', fn.plot, '\n')
cat('write.out      =', write.out, '\n')
cat('output.results =', output.results, '\n')
cat('\n')

### If variable source.dir does not end with "/", append "/" to end of source.dir
source.dir <- appendDirName(source.dir)
start.time <- proc.time()[1]
label      <- read.labels(fn.label)


pred <- read.table(file=fn.pred, sep='\t', header=FALSE, row.names=1, check.names=FALSE, quote='')
if (ncol(pred) > 1) {
  pred <- read.table(file=fn.pred, sep='\t', header=TRUE, row.names=1, check.names=FALSE, quote='')
}
pred <- as.matrix(pred)

eval.data <-  eval.result(label, pred)
pdf(fn.plot)
evaluation.model.plot(label, pred, eval.data, model.type)
tmp <- dev.off()

# The following code writes calculated auroc and aupr- values into a file for testing.
if (write.out == TRUE){
  # Testing only makes sense if dim(pred)[2] > 1
  if(ncol(pred) == 1){
    write.table(t(aucs), file=output.results, quote=FALSE, sep='\t', col.names=FALSE, append=FALSE, row.names="auroc values")
    write.table(t(aucspr), file=output.results, quote=FALSE, sep='\t', col.names=FALSE, append=TRUE, row.names="auprc values")
  }else{
    cat("Only one prediction available, ignoring the write_eval_results option.\n")
  }

}

cat('\nSuccessfully evaluated predictions in ', proc.time()[1] - start.time,
    ' seconds\n', sep='')
