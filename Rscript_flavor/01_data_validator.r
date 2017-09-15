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

source.dir   <- opt$srcdir
fn.in.meta   <- opt$metadata_in
fn.out.meta  <- opt$metadata_out
fn.in.label  <- opt$label_in
fn.out.label <- opt$label_out
fn.in.feat   <- opt$feat_in
fn.out.feat  <- opt$feat_out

# print parameters of the run
cat("=== Paramaters of the run:\n\n")
cat('source.dir   =', source.dir, '\n')
cat('fn.in.meta   =', fn.in.meta, '\n')
cat('fn.out.meta  =', fn.out.meta, '\n')
cat('fn.in.label  =', fn.in.label, '\n') 
cat('fn.out.label =', fn.out.label, '\n')
cat('fn.in.feat   =', fn.in.feat, '\n')
cat('fn.out.feat  =', fn.out.feat, '\n')
cat('\n')

source.dir <- appendDirName(source.dir)
start.time <- proc.time()[1]

# reading in the files
feat  <- read.features(fn.in.feat)
label <- read.labels(fn.in.label,feat)
meta  <- read.meta(fn.in.meta)



### Start Core function
validated.files <- validate.data(feat = feat, label = label$label, meta = meta )
### End Core function

### write validated label, feature and meta-data
# labels
if (!is.null(fn.out.label)){
  write.table(label$header,          file = fn.out.label,  quote = FALSE, sep = '\t', row.names = FALSE,
              col.names = FALSE, append = FALSE)
  write.table(t(as.matrix(validated.files$label)), file = fn.out.label,  quote = FALSE, sep = '\t', row.names = FALSE, ### a bit of a dirty hack?
              col.names = TRUE, append = TRUE)
}
write.table(validated.files$feat,    file = opt$feat_out,   quote = FALSE, sep = '\t', row.names = TRUE,
            col.names = TRUE)
if (!is.null(meta) && !is.null(opt$metadata_out)) {
  write.table(validated.files$meta,  file=opt$metadata_out, quote = FALSE, sep = '\t', row.names = TRUE,
              col.names = NA)
}

cat('\nSuccessfully validated data in ', proc.time()[1] - start.time, ' seconds\n', sep='')
