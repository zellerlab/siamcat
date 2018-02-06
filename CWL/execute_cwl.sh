#!/bin/bash
# Examples how to execute the CWL workflow:

# arguments given to this script
yaml_job_file="$1"
job_name="$2"

location_of_this_script="$( cd "$(dirname "${BASH_SOURCE[0]}")/.." ; pwd -P )"
#cwl_workflow="${location_of_this_script}/CWL/workflows/siamcat_workflow.cwl"
cwl_workflow="${location_of_this_script}/CWL/workflows/siamcat_workflow_scatter_models.cwl"
# If no csv with metadata is present, use the alternative workflow:
#cwl_workflow="${location_of_this_script}/workflows/siamcat_workflow_no_meta.cwl"
working_dir="${PWD}/${job_name}"
output_dir="${working_dir}/out" 
temp_dir="${working_dir}/temp"
base_dir="${working_dir}/base"
cache_dir="${working_dir}/cache"
mkdir $working_dir
mkdir $output_dir $temp_dir $base_dir $cache_dir

# Use the reference implementation cwltool (no parallelizing and no usage of batch systems possible)
cwltool --no-container --debug --tmpdir-prefix=$temp_dir --basedir=$base_dir --outdir=$output_dir --cachedir=$cache_dir $cwl_workflow $yaml_job_file

# Use toil as cwl runner (allows parallelization):
#cwltoil --no-container --logDebug --jobStore=$working_dir/jobstore --workDir=$temp_dir --basedir=$base_dir --outdir=$output_dir  $cwl_workflow $yaml_job_file

# Apply cwltoil with a batch system, here for instance SLURM:
#cwltoil --no-container --batchSystem=Slurm --disableCaching --logDebug --jobStore=$working_dir/jobstore --workDir=$temp_dir --basedir=$base_dir --outdir=$output_dir $cwl_workflow $yaml_job_file

# NOTE: If you wish to use docker, please just remove the "--no-container" flag from above command lines.

# The workflow is tested for both cwltool and toil, however it should work with any other cwl runner.
# Please let us no if it doesn't work (and also if it does work).