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

suppressMessages(library('optparse'))
suppressMessages(library('SIAMCAT'))

# define arguments
option_list = list(
  # make_option('--pkgdir', type='character', help='Source directory of dataprep'),
  make_option('--feat_in',         type='character',                    help='Input file containing features'),
  make_option('--param_out',       type='character',                    help='Output file to which normalization parameters will be written'),
  make_option('--method',          type='character', default='log.std', help='Normalization method (either \"rank.unit\", \"rank.std\", \"log.std\" or \"log.unit\")'),
  make_option('--log_n0',          type='double',    default=10^-8,     help='Pseudocount that is added before log-transformation'),
  make_option('--sd_min_quantile', type='double',    default=0.1,       help='Quantile of the distribution of standard deviation of all feature that will be added to the denominator during standardization of each feature in order to avoid underestimation (only for metod==\"log.std\")'),
  make_option('--vector_norm',     type='integer',   default=2,         help='Vector norm to use (either 1 or 2, only for method==\"log.unit\")'),
  make_option('--norm_margin',     type='integer',   default=1,         help='Margin for normalization (only for method==\"log.unit\"), can be either 1 (for features), 2 (for samples), or 3 (for global normalization)'),
  make_option('--feat_out',        type='character',                    help='Output file to which features after normalization are written')
)

# parse arguments
opt = parse_args(OptionParser(option_list=option_list))
cat("=== 06_generic_normalizer.r\n")
cat("=== Paramaters of the run:\n\n")
cat('opt$feat_in         =', opt$feat_in, '\n')
cat('opt$feat_out        =', opt$feat_out, '\n')
cat('opt$param_out       =', opt$param_out, '\n')
cat('opt$method          =', opt$method, '\n')
cat('opt$log_n0          =', opt$log_n0, '\n')
cat('opt$sd_min_quantile =', opt$sd_min_quantile, '\n')
cat('opt$vector_norm     =', opt$vector_norm, '\n')
cat('opt$norm_margin     =', opt$norm_margin, '\n')
cat('\n')

start.time = proc.time()[1]

### read feature data
feat  <- read.features(opt$feat_in)
### Start core function
normalized.data <- normalize.feat(feat = feat,
                      norm.method = opt$method,
                      norm.param = list(
                      log.n0      = opt$log_n0,
                      sd.min.q    = opt$sd_min_quantile,
                      n.p         = opt$vector_norm,
                      norm.margin = opt$norm_margin))

### End core function


### write output
write.table(normalized.data$feat, file=opt$feat_out, quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)

### write parameters
write('#normalization parameters', opt$param_out, append=FALSE, sep='')
for (p in 1:length(normalized.data$par)) {
  write(paste('###par', names(normalized.data$par)[p], mode(normalized.data$par[[p]]), sep=':'), opt$param_out, append=TRUE, sep='')
  write(normalized.data$par[[p]], opt$param_out, ncolumns=1, append=TRUE, sep='')
}


cat('\nSuccessfully normalized feature data in ', proc.time()[1] - start.time, ' seconds\n', sep='')
