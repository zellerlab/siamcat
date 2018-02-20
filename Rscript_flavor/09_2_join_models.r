#!/usr/bin/Rscript
###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# RScript flavor
#
# written by Georg Zeller
# with additions by Jakob Wirbel, Nicolai Karcher and Konrad Zych
# EMBL Heidelberg 2012-2017
#
# version 0.3.1
# file last updated: 01.02.2018
# GNU GPL 3.0
###

### parse commandline arguments
load("../dependencies/testImage.RData")

# define arguments
option_list <- list(
  make_option('--common_prefix', type='character',  help='Name prefix that all input model files have in common.'),
  make_option('--joined_model_list', type='character',  help='Output which contains all the input models in a single file.')
)

opt         <- parse_args(OptionParser(option_list=option_list))

# print parameters of the run
cat("=== 09.2_joun_models.r\n")
cat("=== Paramaters of the run:\n\n")
cat('common_prefix=', opt$common_prefix, '\n')
cat('joined_model_list=', opt$joined_model_list, '\n')
cat('\n')

start.time <- proc.time()[1]

# read input files:
input.files <- dir('.', paste0("^", opt$common_prefix))

if ( length(input.files) < 1 ){
  cat( paste0( "Error: could not find files with prefix  \"", opt$common_prefix, "\"\n"))
} else {
  cat(paste(length(input.files), "files with the specified prefix were found.\n"))
}

# Load models from files:
input.models <- lapply(input.files, function(file) get(load(file)))

# Read the model type (from the first file, the last item in the list):
# It is assumed that all input models were build using the same method.
model.type <- input.models[[1]][[ length(input.models[[1]]) ]]

# create a joined model list:
models.list <- unlist(lapply(input.models, function(model) model[ 1:length(model)-1 ]), recursive = FALSE)
models.list$model.type <- model.type

# write joined model list to file:
save(models.list , file=opt$joined_model_list)


cat('\n++++++++++++++++++++\nSuccessfully joined models ', proc.time()[1] - start.time,
    ' seconds\n++++++++++++++++++++\n', sep='')
