cwlVersion: v1.0
class: CommandLineTool
baseCommand: Rscript
requirements:
  - class: InlineJavascriptRequirement

inputs:
  srcdir:
    type: string
    inputBinding:
      prefix: --srcdir
      position: 2
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
    - position: 0
      valueFrom: $(inputs.srcdir)/10_make_predictions.r
    - position: 2
      prefix: --pred
      valueFrom: $(inputs.feat_in.nameroot)_predictions.tsv

outputs:
  predictions:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_predictions.tsv
