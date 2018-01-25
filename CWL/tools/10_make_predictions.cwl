cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/siamcat:0.3.1

baseCommand: 10_make_predictions.r

inputs:
  feat_in:
    type: File
    inputBinding:
      prefix: --feat_in
      position: 1
  label_in:
    type: File
    inputBinding:
      prefix: --label_in
      position: 1
  test_sets:
    type: File?
    inputBinding:
      prefix: --test_sets
      position: 1
  model:
    type: File
    inputBinding:
      prefix: --mlr_models_list
      position: 1

arguments:
    - position: 2
      prefix: --pred
      valueFrom: $(inputs.feat_in.nameroot)_predictions.tsv

outputs:
  predictions:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_predictions.tsv
