cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}

baseCommand: 08_split_data.r

inputs:
  label_in:
    type: File
    inputBinding:
      prefix: --label_in
      position: 2
  metadata_in:
    type: File?
    inputBinding:
      prefix: --metadata_in
      position: 2
  num_folds:
    type: int?
    inputBinding:
      prefix: --num_folds
      position: 2
  resample:
    type: int?
    inputBinding:
      prefix: --resample
      position: 2
  stratify:
    type: boolean?
    inputBinding:
      position: 2
      prefix: --stratify
      valueFrom: $(self.toString())
  inseparable:
    type: string?
    inputBinding:
      prefix: --inseparable
      position: 2

arguments:
    - position: 2
      prefix: --train_sets
      valueFrom: $(inputs.feat_in.nameroot)_trainSets.tsv
    - position: 2
      prefix: --test_sets
      valueFrom: $(inputs.feat_in.nameroot)_testSets.tsv

outputs:
  train_sets_out:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_trainSets.tsv
  test_sets_out:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_testSets.tsv
