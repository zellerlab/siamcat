 ###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# RScript flavor
#
# written by Georg Zeller
# with additions by Nicolai Karcher and Konrad Zych
# EMBL Heidelberg 2012-2017
#
# version 0.2.0
# file last updated: 12.06.2017
# GNU GPL 3.0
###

suppressMessages(library('optparse'))
suppressMessages(library('SIAMCAT'))

# define arguments
option_list <- list(
  make_option(c('-s', '--srcdir'),  type='character', help='Source directory of this and other utility scripts'),
  make_option('--metadata_in',      type='character', help='Input file containing meta-data'),
  make_option('--metadata_out',     type='character', default = "meta_validated.tsv",
                                    help='Output file to which validated meta-data is written'),
  make_option('--label_in',         type='character', help='Input file containing labels'),
  make_option('--label_out',        type='character', default = "label_validated.tsv",
                                    help='Output file to which validated label data is written'),
  make_option('--feat_in',          type='character', help='Input file containing features'),
  make_option('--feat_out',         type='character', default = "feat_validated.tsv",
                                    help='Output file to which validated feature data is written')
)

# parse arguments
opt          <- parse_args(OptionParser(option_list=option_list))

# print parameters of the run
cat("=== 01_data_validator.r\n")
cat("=== Paramaters of the run:\n\n")
cat('srcdir       =', opt$srcdir, '\n')
cat('metadata_in  =', opt$metadata_in, '\n')
cat('metadata_out =', opt$metadata_out, '\n')
cat('label_in     =', opt$label_in, '\n')
cat('label_out    =', opt$label_out, '\n')
cat('feat_in      =', opt$feat_in, '\n')
cat('feat_out     =', opt$feat_out, '\n')
cat('\n')

source.dir <- appendDirName(opt$srcdir)
start.time <- proc.time()[1]

# reading in the files
feat  <- read.features(opt$feat_in)
label <- read.labels(opt$label_in,feat)
meta  <- read.meta(opt$metadata_in)



### Start Core function
validated.files <- validate.data(feat = feat, label = label, meta = meta )
### End Core function

### write validated label, feature and meta-data
# labels
if (!is.null(opt$label_out)){
  write.table(label$header,          file = opt$label_out,  quote = FALSE, sep = '\t', row.names = FALSE,
              col.names = FALSE, append = FALSE)
  write.table(t(as.matrix(validated.files$label$label)), file = opt$label_out,  quote = FALSE, sep = '\t', row.names = FALSE, ### a bit of a dirty hack?
              col.names = TRUE, append = TRUE)
}
write.table(validated.files$feat,    file = opt$feat_out,   quote = FALSE, sep = '\t', row.names = TRUE,
            col.names = TRUE)
if (!is.null(meta) && !is.null(opt$metadata_out)) {
write.table(validated.files$meta,  file=opt$metadata_out, quote = FALSE, sep = '\t', row.names = TRUE,
              col.names = NA)
}

cat('\nSuccessfully validated data in ', proc.time()[1] - start.time, ' seconds\n', sep='')
