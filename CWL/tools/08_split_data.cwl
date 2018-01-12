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
      position: 1
  metadata_in:
    type: File?
    inputBinding:
      prefix: --metadata_in
      position: 1
  num_folds:
    type: int?
    inputBinding:
      prefix: --num_folds
      position: 2
  resample:
    type: int?
    default: 0
    inputBinding:
      prefix: --resample
      position: 2
  stratify:
    doc: TRUE or FALSE is allowed
    type: string?
    inputBinding:
      position: 2
      prefix: --stratify
  inseparable:
    type: string?
    inputBinding:
      prefix: --inseparable
      position: 2

arguments:
    - position: 3
      prefix: --train_sets
      valueFrom: $(inputs.label_in.nameroot)_trainSets.tsv
    - position: 3
      prefix: --test_sets
      valueFrom: $(inputs.label_in.nameroot)_testSets.tsv

outputs:
  train_sets_out:
    type: File
    outputBinding:
      glob: $(inputs.label_in.nameroot)_trainSets.tsv
  test_sets_out:
    type: File
    outputBinding:
      glob: $(inputs.label_in.nameroot)_testSets.tsv
