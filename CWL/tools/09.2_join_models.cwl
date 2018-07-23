cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing: $(inputs.input_mlr_models)

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/siamcat:dev

baseCommand: 09.2_join_models.r

inputs:
  input_mlr_models:
    doc: |
      array of RData files containing the model lists 
      that should be contatenated 
    type:
      type: array
      items: File

arguments:
  - prefix: --joined_model_list
    valueFrom: |
      ${ return inputs.input_mlr_models[0].nameroot.replace(/_\d\d\dr_\d\d\df.*/i,".RData") }
    # name for merged output file
  - prefix: --common_prefix
    valueFrom: |
      ${ return inputs.input_mlr_models[0].nameroot.replace(/_\d\d\dr_\d\d\df.*/i,"") }
    position: 1
    # Prefix of files can be set to match everything
    # since only the model files are assumed 
    # to exist in the runtime working dir

outputs:
  joined_model:
    type: File
    outputBinding:
      glob: ${ return inputs.input_mlr_models[0].nameroot.replace(/_\d\d\dr_\d\d\df.*/i,".RData") }