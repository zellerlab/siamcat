###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# RScript flavor
#
# written by Georg Zeller
# with additions by Nicolai Karcher and Konrad Zych
# EMBL Heidelberg 2012-2017
#
# version 0.2.0
# file last updated: 13.06.2017
# GNU GPL 3.0
###

suppressMessages(library('optparse'))
suppressMessages(library('SIAMCAT'))

# define arguments
option_list <- list(
  make_option(c('-s', '--srcdir'), type='character',                   help='Source directory of this and other utility scripts'),
  make_option('--label_in',        type='character',                   help='Input file containing labels'),
  make_option('--feat_in',         type='character',                   help='Input file containing features'),
  make_option('--plot',            type='character',                   help='Output pdf file which will contain resulting plots'),
  make_option('--col_scheme',      type='character', default='RdYlBu', help='Color scheme'),
  make_option('--alpha',           type='double',    default=0.05,     help='Significance level: only features with p-values < alpha will be reported'),
  make_option('--min_fc',          type='double',    default=0,        help='Fold-change cutoff: only features with absolute log-10 fold change > min_fc will be reported'),
  make_option('--mult_test',       type='character', default='fdr',    help='Method to correct for multiple testing (one of \"fdr\", \"holm\", \"bonferroni\", \"BY\", or \"none\")'),
  make_option('--detect_limit',    type='double',    default=10^-8,    help='Lower detection limit for feature values (for log-transform and plots)'),
  make_option('--max_show',        type='integer',   default=50,       help='Maximum number of significant features to be shown in result plots'),
  make_option('--plot_type',       type='character', default='bean',   help = 'One of \"quantile.box\", \"box\", \"bean\", \"quantile.rect\"')
)

# parse arguments
opt          <- parse_args(OptionParser(option_list=option_list))
source.dir   <- opt$srcdir
fn.in.label  <- opt$label_in
fn.in.feat   <- opt$feat_in
fn.plot      <- opt$plot
color.scheme <- opt$col_scheme
alpha        <- opt$alpha
min.fc       <- opt$min_fc
mult.corr    <- opt$mult_test
detect.lim   <- opt$detect_limit
max.show     <- opt$max_show
plot.type    <- opt$plot_type

cat("=== 05_association_check.r\n")
cat("=== Paramaters of the run:\n\n")
cat('source.dir   =', source.dir, '\n')
cat('fn.in.label  =', fn.in.label, '\n')
cat('fn.in.feat   =', fn.in.feat, '\n')
cat('fn.plot      =', fn.plot, '\n')
cat('color.scheme =', color.scheme, '\n')
cat('alpha        =', alpha, '\n')
cat('min.fc       =', min.fc, '\n')
cat('mult.corr    =', mult.corr, '\n')
cat('detect.lim   =', detect.lim, '\n')
cat('max.show     =', max.show, '\n')
cat('plot.type    =', plot.type, '\n')
cat('\n')

### If variable source.dir does not end with "/", append "/" to end of source.dir
source.dir <- appendDirName(source.dir)
start.time <- proc.time()[1]

feat  <- read.features(fn.in.feat)
label <- read.labels(fn.in.label, feat)
check.associations(feat, label, fn.plot, color.scheme, alpha, min.fc, mult.corr, detect.lim, max.show, plot.type)


cat('\nSuccessfully analyzed statistically significant associations between individual features and labels in ', proc.time()[1] - start.time, ' seconds\n', sep='')
