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

# parameters that cannot be specified in the interface
# consens.thres = 0.5 # Exposed
norm.models <- FALSE
# TODO determine these automatically
z.score.lim <- c(-3,3)
fc.lim      <- c(-5,5)
detect.lim  <- 10^-8

# define arguments
option_list  <-  list(
  make_option(c('-s', '--srcdir'), type='character', help='Source directory of this and other utility scripts'),
  make_option('--label', type='character', help='Input file containing labels'),
  make_option('--feat', type='character', help='Input file containing features'),
  make_option('--origin_feat', type='character', help='Input file containing unnormalized/filtered features'),
  make_option('--meta', type='character', default='NULL', help='Input file containing metadata'),
  make_option('--model', type='character', help='Input file containing the trained classification model(s)'),
  make_option('--pred', type='character', help='Input file containing the trained classification model(s)'),
  #make_option('--plot', type='character', help='Output file for plotting'),
  make_option('--col_scheme', type='character', default='RdYlBu', help='Color scheme'),
  make_option('--heatmap_type', type='character', default='zscore', help='which metric should be used to plot feature changes in heatmap? (zscore|fc)'),
  make_option('--consens_thres', type='double', default=0.5, help='specifies the minimal ratio of models incorporating a feature to include it into heatmap'),
  make_option('--plot', type='character', help='Output file for plotting')
)


# parse arguments
opt            <- parse_args(OptionParser(option_list=option_list))
source.dir     <- opt$srcdir
fn.label       <- opt$label
fn.feat        <- opt$feat
fn.origin.feat <- opt$origin_feat
fn.meta        <- opt$meta
fn.model       <- opt$model
fn.pred        <- opt$pred
fn.plot        <- opt$plot
color.scheme   <- opt$col_scheme
heatmap.type   <- opt$heatmap_type
consens.thres  <- opt$consens_thres

cat("=== 12_model_interpretor.r\n")
cat("=== Paramaters of the run:\n\n")
cat('source.dir     =', source.dir, '\n')
cat('fn.feat        =', fn.feat, '\n')
cat('fn.origin.feat =', fn.origin.feat, '\n')
cat('fn.label       =', fn.label, '\n')
cat('fn.model       =', fn.model, '\n')
cat('fn.pred        =', fn.pred, '\n')
cat('fn.meta        =', fn.meta, '\n')
cat('fn.plot        =', fn.plot, '\n')
cat('color.scheme   =', color.scheme, '\n')
cat('consens.thres  =', consens.thres, '\n')
cat('\n')

### If variable source.dir does not end with "/", append "/" to end of source.dir
source.dir <- appendDirName(source.dir)
# TODO (remove last check)
# optional parameters will be reset to NULL if specified as 'NULL', 'NONE' or 'UNKNOWN'
if (toupper(fn.meta)=='NULL' || toupper(fn.meta)=='NONE' || toupper(fn.meta)=='UNKNOWN' || fn.meta=='NULL.tsv') {
  fn.meta = NULL
  cat('fn.meta not given: no metadata to display\n')
}
start.time <- proc.time()[1]

### read features and labels
# features
feat        <- read.features(fn.feat)
origin.feat <- read.features(fn.origin.feat)
label       <- read.labels(fn.label,feat)



### load trained model(s)
model        <- NULL
model$W      <- read.table(file=fn.model, sep='\t', header=TRUE, row.names=1, stringsAsFactors=FALSE, check.names=FALSE, quote='')
stopifnot(nrow(model$W) == nrow(feat))
# parse model header
con          <- file(fn.model, 'r')
model.header <- readLines(con, 1)
close(con)
model.header <- parse.model.header(model.header)
# TODO compare label and model headers
#cat(label.header, '\n')
#cat(model.header$label.header, '\n')
#stopifnot(substr(label.header,2,length(label.header)) == model.header$label.header)

### load predictions
# TODO compare prediction and label headers
pred         <- read.table(file=fn.pred, sep='\t', header=FALSE, row.names=1, check.names=FALSE, quote='')
if (dim(pred)[2] > 1) {
  pred           <- read.table(file=fn.pred, sep='\t', header=TRUE, row.names=1, check.names=FALSE, quote='')
}
pred         <- as.matrix(pred)
#cat(dim(pred), '\n')

### make sure that label and prediction are in the same order
stopifnot(all(names(label$label) %in% rownames(pred)) && all(rownames(pred) %in% names(label$label)))
m            <- match(names(label$label), rownames(pred))
pred         <- pred[m,,drop=FALSE]
stopifnot(all(names(label$label) == rownames(pred)))

### load metadata if provided
meta         <- read.meta(fn.meta)
stopifnot(all(names(label$label) == rownames(meta)))

pdf(fn.plot, paper='special', height=8.27, width=11.69) # format: A4 landscape
interpretor.model.plot(feat, label, meta, model, pred, color.scheme, consens.thres)
tmp <- dev.off()
cat('\nSuccessfully interpreted model in ', proc.time()[1] - start.time,
    ' seconds\n', sep='')
