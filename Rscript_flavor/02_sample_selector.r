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

source.dir    <- opt$srcdir
fn.in.meta    <- opt$metadata_in
fn.out.meta   <- opt$metadata_out
fn.in.label   <- opt$label_in
fn.out.label  <- opt$label_out
fn.in.feat    <- opt$feat_in
fn.out.feat   <- opt$feat_out
filter        <- opt$filter_var
allowed.range <- opt$allowed_range


# print parameters of the run
cat("=== Paramaters of the run:\n\n")
cat('source.dir   =', source.dir, '\n')
cat('fn.in.meta   =', fn.in.meta, '\n')
cat('fn.out.meta  =', fn.out.meta, '\n')
cat('fn.in.label  =', fn.in.label, '\n')
cat('fn.out.label =', fn.out.label, '\n')
cat('fn.in.feat   =', fn.in.feat, '\n')
cat('fn.out.feat  =', fn.out.feat, '\n')
cat('filter       =', filter, '\n')

source.dir  <- appendDirName(source.dir)
start.time  <- proc.time()[1]


### read label, feature and meta- data
# features
# reading in the files
feat  <- read.features(fn.in.feat)
label <- read.labels(fn.in.label, feat)
meta  <- read.meta(fn.in.meta)
stopifnot(all(names(label$label) == rownames(meta)))


### select samples fulfilling the filter criteria
# (i.e. samples having certain metadata values)
results     <-  select.samples(meta=meta, feat=feat, label=label$label, filter=filter, allowed.range=allowed.range)



### write label, feature and meta-data with selected sample set
# labels
write(label$header,        file=fn.out.label, append=FALSE)
write.table(t(as.matrix(results$label)), file=fn.out.label, quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE, append=TRUE)
# features
write.table(results$feat,  file=fn.out.feat, quote=FALSE,  sep='\t', row.names=TRUE, col.names=TRUE)
# meta-data
write.table(results$meta,  file=fn.out.meta, quote=FALSE,  sep='\t', row.names=TRUE, col.names=TRUE)


cat('\nSuccessfully selected samples in ', proc.time()[1] - start.time, ' seconds\n', sep='')
