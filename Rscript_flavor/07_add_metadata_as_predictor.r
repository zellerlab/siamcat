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
option_list <- list(
  make_option('--feat_in', type='character', help='Input file containing features'),
  make_option('--metadata_in', type='character', help='Input file containing metadata'),
  make_option('--pred_names', type='character', help='names (comma-separated list) of the metavariables to be added to the feature matrix as predictors'),
  make_option('--std_meta', type='logical', default=TRUE, help='Shall added (metadata) features be standardized?'),
  make_option('--feat_out', type='character', help='Output file to which features after selection are written')
)


# parse arguments
opt         <- parse_args(OptionParser(option_list=option_list))
cat("=== 07_meta_predictor_adder.r\n")
cat("=== Paramaters of the run:\n\n")
cat('feat_in     =', opt$feat_in, '\n')
cat('metadata_in =', opt$metadata_in, '\n')
cat('pred_names  =', opt$pred_names, '\n')
cat('std_meta    =', opt$std_meta, '\n')
cat('feat_out    =', opt$feat_out, '\n')
cat('\n')

start.time  <- proc.time()[1]

#### read features and metadata
feat        <- read.features(opt$feat_in)
meta        <- read.meta( opt$metadata_in)
siamcat <- construct.siamcat(feat,meta)

pred.names <- strsplit(opt$pred_names, ',', fixed=TRUE)[[1]]
siamcat     <- add.meta.pred(siamcat,
	                        pred.names=pred.names,
	                        std.meta=opt$std_meta)

### write combined feature table
write.table(siamcat@phyloseq@otu_table, file=opt$feat_out, quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)

cat('\nSuccessfully added metadata to features in ', proc.time()[1] - start.time, ' seconds\n', sep='')
