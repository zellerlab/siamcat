cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/siamcat:dev

baseCommand: 09_train_models.r

inputs:
  feat_in:
    type: File
    inputBinding:
      prefix: --feat_in
      position: 1
  label_in:
    type: File?
    inputBinding:
      prefix: --label_in
      position: 1
  train_sets:
    type: File?
    inputBinding:
      prefix: --train_sets
      position: 1
  train_method:
    type: string?
    inputBinding:
      prefix: --method
      position: 2
  stratify:
    doc: TRUE or FALSE is allowed
    type: string?
    inputBinding:
      prefix: --stratify
      position: 2
  sel_criterion:
    type: string?
    inputBinding:
      prefix: --sel_criterion
      position: 2
  min_nonzero_coeff:
    type: int?
    inputBinding:
      prefix: --min_nonzero_coeff
      position: 2
  param_set:
    type: string?
    inputBinding:
      prefix: --param_set
      position: 2

arguments:
    - position: 3
      prefix: --mlr_models_list
      valueFrom: $(inputs.train_sets.nameroot)_model.RData

outputs:
  model:
    type: File
    outputBinding:
      glob: $(inputs.train_sets.nameroot)_model.RData
