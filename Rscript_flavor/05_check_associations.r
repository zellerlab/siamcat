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
cat("=== 05_association_check.r\n")
cat("=== Paramaters of the run:\n\n")
cat('srcdir       =', opt$srcdir, '\n')
cat('label_in     =', opt$label_in, '\n')
cat('feat_in      =', opt$feat_in, '\n')
cat('plot         =', opt$plot, '\n')
cat('col_scheme   =', opt$col_scheme, '\n')
cat('alpha        =', opt$alpha, '\n')
cat('min_fc       =', opt$min_fc, '\n')
cat('mult_test    =', opt$mult_test, '\n')
cat('detect_limit =', opt$detect_limit, '\n')
cat('max_show     =', opt$max_show, '\n')
cat('plot_type    =', opt$plot_type, '\n')
cat('\n')

### If variable source.dir does not end with "/", append "/" to end of source.dir
source.dir <- appendDirName(source.dir)
start.time <- proc.time()[1]

feat  <- read.features(fn.in.feat)
label <- read.labels(fn.in.label, feat)
check.associations(feat = feat, label = label, fn.plot = opt$plot, color.scheme = opt$col_scheme, alpha = opt$alpha, 
                   mult.corr = opt$mult_test, detect.lim = opt$detect_limit, max.show = opt$max_show,
                   plot.type = opt$plot_type)


cat('\nSuccessfully analyzed statistically significant associations between individual features and labels in ', proc.time()[1] - start.time, ' seconds\n', sep='')
