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

# cwltool:
# --------
# (no parallelization allowed)

# ignore containers
cwltool --no-container --debug --tmpdir-prefix=$temp_dir --basedir=$base_dir --outdir=$output_dir --cachedir=$cache_dir $cwl_workflow $yaml_job_file

# use docker containers:
cwltoil --logDebug --jobStore=$working_dir/jobstore --workDir=$temp_dir --basedir=$base_dir --outdir=$output_dir  $cwl_workflow $yaml_job_file

# use udocker containers:
cwltoil --user-space-docker-cmd=/home/breuerk/udocker --logDebug --jobStore=$working_dir/jobstore --workDir=$temp_dir --basedir=$base_dir --outdir=$output_dir  $cwl_workflow $yaml_job_file

# cwltoil:
# --------
# (parallelization possible)

# Example using Slurm and udocker:
cwltoil --user-space-docker-cmd=/home/breuerk/udocker --batchSystem=Slurm --disableCaching --logDebug --jobStore=$working_dir/jobstore --workDir=$temp_dir --basedir=$base_dir --outdir=$output_dir $cwl_workflow $yaml_job_file

# The flag "--no-container" exists, too. Like for cwltool: if neither "--user-space-docker-cmd" nor "--no-container" is set, toil will try to run docker
