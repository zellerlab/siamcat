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
option_list  <-  list(
  make_option(c('-s', '--srcdir'), type='character',                   help='Source directory of this and other utility scripts'),
  make_option('--label_in',        type='character',                   help='Input file containing labels'),
  make_option('--feat_in',         type='character',                   help='Input file containing features'),
  make_option('--origin_feat',     type='character',                   help='Input file containing unnormalized/filtered features'),
  make_option('--metadata_in',     type='character', default=NULL,   help='Input file containing metadata'),
  make_option('--mlr_models_list', type='character',                   help='Input RData file containing the output from train.model'),
  make_option('--pred',            type='character',                   help='Input file containing the trained classification model(s)'),
  make_option('--col_scheme',      type='character', default='RdYlBu', help='Color scheme'),
  make_option('--heatmap_type',    type='character', default='zscore', help='which metric should be used to plot feature changes in heatmap?
                                                                             (zscore|fc)'),
  make_option('--consens_thres',   type='double',    default=0.5,      help='specifies the minimal ratio of models incorporating a feature
                                                                             to include it into heatmap'),
  make_option('--plot',            type='character',                   help='Output file for plotting')
)


# parse arguments
opt            <- parse_args(OptionParser(option_list=option_list))
cat("=== 12_model_interpretor.r\n")
cat("=== Paramaters of the run:\n\n")
cat('srcdir          =', opt$srcdir, '\n')
cat('feat            =', opt$feat, '\n')
cat('origin_feat.    =', opt$origin_feat, '\n')
cat('metadata_in     =', opt$metadata_in, '\n')
cat('label           =', opt$label, '\n')
cat('mlr_models_list =', opt$mlr_models_list, '\n')
cat('pred            =', opt$pred, '\n')
cat('plot            =', opt$plot, '\n')
cat('col_scheme      =', opt$col_scheme, '\n')
cat('consens_thres   =', opt$consens_thres, '\n')
cat('heatmap_type    =', opt$heatmap_type, '\n')
cat('\n')

### If variable source.dir does not end with "/", append "/" to end of source.dir
source.dir <- appendDirName(opt$srcdir)
# TODO (remove last check)
# optional parameters will be reset to NULL if specified as 'NULL', 'NONE' or 'UNKNOWN'


### read features and labels
start.time <- proc.time()[1]
feat        <- read.features(opt$feat)
label       <- read.labels(opt$label,feat)
origin.feat <- read.features(opt$origin_feat)
if (is.null(opt$metadata_in)) {
  cat('meta not given: no metadata to display\n')
} else {
  meta         <- read.meta(opt$metadata_in)
  stopifnot(all(names(label$label) == rownames(meta)))
}




### load trained model(s)
load(opt$mlr_models_list) ##loads plm.out
# model        <- NULL
# model$W      <- read.table(file=opt$model, sep='\t', header=TRUE, row.names=1, stringsAsFactors=FALSE, check.names=FALSE, quote='')
# stopifnot(nrow(model$W) == nrow(feat))
# parse model header

pred         <- read.table(file=opt$pred, sep='\t', header=TRUE, row.names=1, check.names=FALSE, comment.char="#")
pred         <- as.matrix(pred)

### make sure that label and prediction are in the same order
stopifnot(all(names(label$label) %in% rownames(pred)) && all(rownames(pred) %in% names(label$label)))
m            <- match(names(label$label), rownames(pred))
pred         <- pred[m,,drop=FALSE]
stopifnot(all(names(label$label) == rownames(pred)))


interpretor.model.plot(feat=feat,
                       label=label,
                       meta=meta,
                       model=models.list,
                       pred=pred,
                       color.scheme=opt$col_scheme,
                       consens.thres=opt$consens_thres,
                       heatmap.type=opt$heatmap_type,
                       fn.plot=opt$plot)

cat('\nSuccessfully interpreted model in ', proc.time()[1] - start.time,
    ' seconds\n', sep='')
