#!/usr/bin/Rscript
###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# RScript flavor
#
# written by Georg Zeller
# with additions by Jakob Wirbel, Nicolai Karcher and Konrad Zych
# EMBL Heidelberg 2012-2018
# GNU GPL 3.0
###

### parse commandline arguments
suppressMessages(library('optparse'))
suppressMessages(library('SIAMCAT'))
# define arguments
option_list <- list(
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
cat('metadata_in   =', opt$metadata_in, '\n')
cat('metadata_out  =', opt$metadata_out, '\n')
cat('label_in      =', opt$label_in, '\n')
cat('label_out     =', opt$label_out, '\n')
cat('feat_in       =', opt$feat_in, '\n')
cat('feat_out      =', opt$feat_out, '\n')
cat('filter_var    =', opt$filter_var, '\n')
cat('allowed_range =', opt$allowed_range, '\n')
cat('allowed_set   =', opt$allowed_set, '\n')

start.time  <- proc.time()[1]


### read label, feature and meta- data
feat  <- read.features(opt$feat_in)
label <- read.labels(opt$label_in)
meta  <- read.meta(opt$metadata_in)
siamcat <- siamcat(feat,label,meta)


### select samples fulfilling the filter criteria
# (i.e. samples having certain metadata values)
siamcat     <-  select.samples(siamcat,
                               filter=opt$filter_var,
                               allowed.range=opt$allowed_range,
                               allowed.set=opt$allowed_set)


### write label, feature and meta-data with selected sample set
# labels
write(label@header,        file=opt$label_out, append=FALSE)
write.table(t(as.matrix(siamcat@label@label)), file=opt$label_out, quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE, append=TRUE)
# features
write.table(siamcat@phyloseq@otu_table,  file=opt$feat_out, quote=FALSE,  sep='\t', row.names=TRUE, col.names=TRUE)
# meta-data
write.table(siamcat@phyloseq@sam_data,  file=opt$metadata_out, quote=FALSE,  sep='\t', row.names=TRUE, col.names=TRUE)


cat('\nSuccessfully selected samples in ', proc.time()[1] - start.time, ' seconds\n', sep='')
