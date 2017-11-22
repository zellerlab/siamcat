###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# RScript flavor
#
# written by Georg Zeller
# with additions by Nicolai Karcher and Konrad Zych
# EMBL Heidelberg 2012-2017
#
# version 0.2.0
# file last updated: 19.05.2017
# GNU GPL 3.0
###

suppressMessages(library('optparse'))
suppressMessages(library('SIAMCAT'))

# define arguments
option_list <- list(
  make_option('--feat_in', type='character', help='Input file containing features'),
  make_option('--metadata_in', type='character', help='Input file containing metadata'),
  make_option('--pred_names', type='character', help='names (comma-separated list) of the metavariables to be added to the feature matrix as predictors'),
  make_option('--std_meta', type='logical', default=TRUE, help='Shall added (metadata) features be standardized?'),
  make_option('--feat_out', type='character', help='Output file to which features after selection are written')
)


# parse arguments
opt         <- parse_args(OptionParser(option_list=option_list))
fn.in.feat  <- opt$feat_in
fn.in.meta  <- opt$metadata_in
fn.out.feat <- opt$feat_out
pred.names  <- opt$pred_names
std.meta    <- opt$std_meta

cat("=== 07_meta_predictor_adder.r\n")
cat("=== Paramaters of the run:\n\n")
cat('fn.in.feat  =', fn.in.feat, '\n')
cat('fn.in.meta  =', fn.in.meta, '\n')
cat('fn.out.feat =', fn.out.feat, '\n')
cat('pred.names  =', pred.names, '\n')
cat('std.meta    =', std.meta, '\n')
cat('\n')

start.time  <- proc.time()[1]


#### read features and metadata
feat        <- read.features(fn.in.feat)
meta        <- read.meta(fn.in.meta)
stopifnot(all(colnames(feat) == rownames(meta)))

pred.names <- strsplit(pred.names, ',', fixed=TRUE)[[1]]
feat       <- add.meta.pred(feat, meta, pred.names, std.meta)

### write combined feature table
write.table(feat, file=fn.out.feat, quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)

cat('\nSuccessfully added metadata to features in ', proc.time()[1] - start.time, ' seconds\n', sep='')
