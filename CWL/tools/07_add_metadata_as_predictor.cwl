cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}

baseCommand: 07_add_metadata_as_predictor.r

inputs:
  feat_in:
    type: File
    inputBinding:
      prefix: --feat_in
      position: 1
  metadata_in:
    type: File
    inputBinding:
      prefix: --metadata_in
      position: 1
  pred_names:
    type: string
    inputBinding:
      prefix: --pred_names
      position: 2
  std_meta:
    type: boolean?
    inputBinding:
      prefix: --std_meta
      position: 3
      valueFrom: $(self.toString())

arguments:
  - prefix: --feat_out
    valueFrom: $(inputs.feat_in.nameroot)_metaAsPred.tsv
    position: 4
    
outputs:
  feat_out:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_metaAsPred.tsv
