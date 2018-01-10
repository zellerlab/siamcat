cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}

baseCommand: 10_make_predictions.r

inputs:
  feat_in:
    type: File
    inputBinding:
      prefix: --feat_in
      position: 2
  label_in:
    type: File?
    inputBinding:
      prefix: --label_in
      position: 2
  test_sets:
    type: File?
    inputBinding:
      prefix: --test_sets
      position: 2
  model_rdata:
    type: File
    inputBinding:
      prefix: --mlr_models_list
      position: 2
  model_matrix:
    type: File
    inputBinding:
      prefix: --model_matrix
      position: 2

arguments:
    - position: 2
      prefix: --pred
      valueFrom: $(inputs.feat_in.nameroot)_predictions.tsv

outputs:
  predictions:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_predictions.tsv
