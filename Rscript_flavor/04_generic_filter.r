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

suppressMessages(library('optparse'))
suppressMessages(library('SIAMCAT'))

# define arguments
option_list = list(
  # make_option('--pkgdir', type='character', help='Source directory of dataprep'),
  make_option('--feat_in', type='character', help='Input file containing features'),
  make_option('--method', type='character', default='abundance', help='Filtering method (one of \"abundance\", \"cum.abundance\", or \"prevalence\")'),
  make_option('--cutoff', type='double', default=0.001, help='abundance / prevalence cutoff applied for filtering'),
  make_option('--recomp_prop', type='logical', default=FALSE, help='Should relative abundances be be recomputed?'),
  make_option('--rm_unmapped', type='logical', default=TRUE, help='Should the abundance of unmapped reads be removed?'),
  make_option('--feat_out', type='character', help='Output file to which features after selection are written')
)

# parse arguments
opt         <-  parse_args(OptionParser(option_list=option_list))

fn.in.feat    <- opt$feat_in
fn.out.feat   <- opt$feat_out
method        <- opt$method
cutoff        <- opt$cutoff
recomp.prop   <- opt$recomp_prop
rm.unmapped   <- opt$rm_unmapped

# print parameters of the run
cat("=== Paramaters of the run:\n\n")
cat('opt$feat_in     =', fn.in.feat, '\n')
cat('opt$feat_out    =', fn.out.feat, '\n')
cat('opt$method      =', method, '\n')
cat('opt$cutoff      =', cutoff, '\n')
cat('opt$recomp_prop =', recomp.prop, '\n')
cat('opt$rm_unmapped =', rm.unmapped, '\n')
cat('\n')

start.time  <- proc.time()[1]

#### read feature data and preprocess
feat        <- read.features(fn.in.feat)

### Start core function
filtered.data <- filter.feat(feat = feat, filter.method = method, cutoff = cutoff, recomp.prop = recomp.prop,
                   rm.unmapped = rm.unmapped)
### End core function

# write filtered feature table
write.table(filtered.data, file=fn.out.feat, quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)

cat('\nSuccessfully filtered feature data in ', proc.time()[1] - start.time, ' seconds\n', sep='')
