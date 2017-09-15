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
option_list <- make_normalizer_options()
option_list[[length(option_list)+1]] <- make_option('--feat_out', type='character', help='Output file to which features after normalization are written')

# parse arguments
opt = parse_args(OptionParser(option_list=option_list))
cat("=== Paramaters of the run:\n\n")
cat('opt$feat_in         =', opt$feat_in, '\n')
cat('opt$feat_out        =', opt$feat_out, '\n')
cat('opt$param_out       =', opt$param_out, '\n')
cat('opt$method          =', opt$method, '\n')
cat('opt$log_n0          =', opt$log_n0, '\n')
cat('opt$sd_min_quantile =', opt$sd_min_quantile, '\n')
cat('opt$vector_norm     =', opt$vector_norm, '\n')
cat('opt$norm_sample     =', opt$norm_sample, '\n')
cat('opt$norm_feature    =', opt$norm_feature, '\n')
cat('opt$norm_global     =', opt$norm_global, '\n')
cat('\n')

start.time = proc.time()[1]

### read feature data
feat  <- read.features(opt$feat_in)
### Start core function
normalized.data <- normalize.feat(feat = feat,
                      norm.method = opt$method, 
                      log.n0      = opt$log_n0, 
                      sd.min.q    = opt$sd_min_quantile, 
                      n.p         = opt$vector_norm, 
                      n.sample    = opt$norm_sample, 
                      n.feature   = opt$norm_feature, 
                      n.global    = opt$norm_global)

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
