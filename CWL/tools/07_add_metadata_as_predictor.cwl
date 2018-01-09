cwlVersion: v1.0
class: CommandLineTool
baseCommand: Rscript
requirements:
  - class: InlineJavascriptRequirement

inputs:
  feat_in:
    type: File
    inputBinding:
      prefix: --feat_in
      position: 2
  metadata_in:
    type: File
    inputBinding:
      prefix: --metadata_in
      position: 2
  pred_names:
    type: string
    inputBinding:
      prefix: --pred_names
      position: 2
  std_meta:
    type: boolean?
    inputBinding:
      prefix: --std_meta
      position: 2
      valueFrom: $(self.toString())

arguments:
  - position: 0
    valueFrom: $(inputs.srcdir)/07_add_metadata_as_predictor.r
  - prefix: --feat_out
    valueFrom: $(inputs.feat_in.nameroot)_metaAsPred.tsv
    position: 2
    
outputs:
  feat_out:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_metaAsPred.tsv
