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
option_list    <- list(
  make_option(c('-s', '--srcdir'), type='character', help='Source directory of this and other utility scripts'),
  make_option('--metadata_in',     type='character', help='Input file containing meta-data'),
  make_option('--label_in',        type='character', help='Input file containing labels'),
  make_option('--plot',            type='character', help='Output pdf file which will contain resulting plots')
  )

# parse arguments
opt          <- parse_args(OptionParser(option_list=option_list))
# print parameters of the run
cat("=== 03_confounder_check.r\n")
cat("=== Paramaters of the run:\n\n")
cat('source.dir   =', opt$srcdir, '\n')
cat('metadata_in  =', opt$metadata_in, '\n')
cat('label_in     =', opt$label_in, '\n')
cat('plot         =', opt$plot, '\n')
cat('\n')

source.dir   <- appendDirName(opt$srcdir)
start.time   <- proc.time()[1]


### read label and meta- data
label        <- read.labels(opt$label_in)
meta         <- read.meta(opt$metadata_in)
if(any(names(label$label) != rownames(meta))) stop("Label names do not match metadata names! Stopping.\n")

m            <- match(names(label$label), rownames(meta))
meta         <- meta[m,]

### Start core function
confounder.check(meta = meta, label = label, fn.plot=opt$plot)
### End core function


cat('\nSuccessfully analyzed meta-data for potential confounding in ', proc.time()[1] - start.time, ' seconds\n', sep='')
