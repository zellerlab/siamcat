#!/bin/bash
# run in cwl_siamcat_playground docker container
# build images: cd ../../dockerfiles/cwl_siamcat_workflow && docker build -t "cwl_siamcat_playground:latest" .
# start container interactively: docker run -d  -it -v "/home/kbreuer/:/home/kbreuer/" --name cwl_siamcat_playground cwl_siamcat_playground
# attach to running container: docker attach --detach-keys ctrl-d cwl_siamcat_playground
# detach and leave running by pressing: crtl-d

rel_out_path="$1" 

working_dir="/home/kbreuer/cwl_test_dir/${rel_out_path}"
output_dir="${working_dir}/out" 
temp_dir="${working_dir}/temp"
base_dir="${working_dir}/base"
cache_dir="${working_dir}/cache"
cwl_dir="/home/kbreuer/siamcat-dev/CWL"
input_dir="/home/kbreuer/test_data_siamcat"
yml_job_filepath="${cwl_dir}/test_jobs/cancer_test_no_options.yml"
workflow_file="${cwl_dir}/workflows/siamcat_workflow_clean.cwl"

mkdir $working_dir
mkdir $output_dir $temp_dir $base_dir $cache_dir
echo "Workflow file: $workflow_file; job file: $yml_job_filepath"

#cwltool --debug --tmpdir-prefix=$temp_dir --basedir=$base_dir --outdir=$output_dir --cachedir=$cache_dir $workflow_file $yml_job_filepath
cwltoil --logDebug --maxCores 8 --jobStore=$working_dir/jobstore --workDir=$temp_dir --basedir=$base_dir --outdir=$output_dir  $workflow_file $yml_job_filepath







