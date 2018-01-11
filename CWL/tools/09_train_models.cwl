cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}

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
  method:
    type: string?
    inputBinding:
      prefix: --method
      position: 2
  stratify:
    type: boolean?
    inputBinding:
      prefix: --stratify
      valueFrom: $(self.toString())
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

arguments:
    - position: 3
      prefix: --model
      valueFrom: $(inputs.feat_in.nameroot)_model.tsv
    - position: 3
      prefix: --mlr_models_list
      valueFrom: $(inputs.feat_in.nameroot)_model.RData
    - position: 3
      prefix: --model_matrix
      valueFrom: $(inputs.feat_in.nameroot)_modelMatrix.txt

outputs:
  model_tsv:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_model.tsv
  model_rdata:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_model.RData
  model_matrix:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_modelMatrix.txt
