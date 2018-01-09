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
    type: File
    inputBinding:
      prefix: --label_in
      position: 2
  method:
    type: string?
    inputBinding:
      prefix: --method
      position: 2
  train_sets:
    type: File?
    inputBinding:
      prefix: --train_sets
      position: 2
  num_folds:
    type: int
    inputBinding:
      prefix: --num_folds
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
    - position: 0
      valueFrom: $(inputs.srcdir)/09_train_models.r
    - position: 2
      prefix: --model
      valueFrom: $(inputs.feat_in.nameroot)_model.tsv
    - position: 2
      prefix: --mlr_models_list
      valueFrom: $(inputs.feat_in.nameroot)_model.RData
    - position: 2
      prefix: --model_matrix
      valueFrom: $(inputs.feat_in.nameroot)_modelMatrix.txt

outputs:
  trained_model:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_model.tsv
  trained_model_rdata:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_model.RData
  trained_model_matrix:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_modelMatrix.txt
