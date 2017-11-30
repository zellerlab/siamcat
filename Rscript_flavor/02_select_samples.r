###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# RScript flavor
#
# written by Georg Zeller
# with additions by Nicolai Karcher and Konrad Zych
# EMBL Heidelberg 2012-2017
#
# version 0.2.0
# file last updated: 14.06.2017
# GNU GPL 3.0
###

### parse commandline arguments
suppressMessages(library('optparse'))
suppressMessages(library('SIAMCAT'))
# define arguments
option_list <- list(
    make_option(c('-s', '--srcdir'), type='character',               help='Source directory of this and other utility scripts'),
    make_option('--metadata_in',     type='character',               help='Input file containing meta-data'),
    make_option('--metadata_out',    type='character',               help='Output file to which meta-data after selection is written'),
    make_option('--label_in',        type='character',               help='Input file containing labels'),
    make_option('--label_out',       type='character',               help='Output file to which labels after selection are written'),
    make_option('--feat_in',         type='character',               help='Input file containing features'),
    make_option('--feat_out',        type='character',               help='Output file to which features after selection are written'),
    make_option('--filter_var',      type='character',               help='Meta-data variable on which filtering is based'),
    make_option('--allowed_range',   type='character', default=NULL, help='Range of allowed values (closed interval)'),
    make_option('--allowed_set',     type='character', default=NULL, help='Set of allowed values (comma-separated list)')
)
# parse arguments
opt           <- parse_args(OptionParser(option_list=option_list))
# print parameters of the run
cat("=== 02_sample_selector.r\n")
cat("=== Paramaters of the run:\n\n")
cat('srcdir        =', opt$srcdir, '\n')
cat('metadata_in   =', opt$metadata_in, '\n')
cat('metadata_out  =', opt$metadata_out, '\n')
cat('label_in      =', opt$label_in, '\n')
cat('label_out     =', opt$label_out, '\n')
cat('feat_in       =', opt$feat_in, '\n')
cat('feat_out      =', opt$feat_out, '\n')
cat('filter_var    =', opt$filter_var, '\n')
cat('allowed_range =', opt$allowed_range, '\n')
cat('allowed_set   =', opt$allowed_set, '\n')

source.dir  <- appendDirName(opt$srcdir)
start.time  <- proc.time()[1]


### read label, feature and meta- data
feat  <- read.features(opt$feat_in)
label <- read.labels(opt$label_in, feat)
meta  <- read.meta(opt$metadata_in)
stopifnot(all(names(label$label) == rownames(meta)))


### select samples fulfilling the filter criteria
# (i.e. samples having certain metadata values)
results     <-  select.samples(meta=meta, 
                               feat=feat, 
                               label=label$label, 
                               filter=opt$filter_var, 
                               allowed.range=opt$allowed_range, 
                               allowed.set=opt$allowed_set)


### write label, feature and meta-data with selected sample set
# labels
write(label$header,        file=opt$label_out, append=FALSE)
write.table(t(as.matrix(results$label)), file=opt$label_out, quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE, append=TRUE)
# features
write.table(results$feat,  file=opt$feat_out, quote=FALSE,  sep='\t', row.names=TRUE, col.names=TRUE)
# meta-data
write.table(results$meta,  file=opt$metadata_out, quote=FALSE,  sep='\t', row.names=TRUE, col.names=TRUE)


cat('\nSuccessfully selected samples in ', proc.time()[1] - start.time, ' seconds\n', sep='')
