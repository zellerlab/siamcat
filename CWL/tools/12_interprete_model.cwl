class: CommandLineTool
cwlVersion: v1.0

requirements:
  InlineJavascriptRequirement: {}

baseCommand: 12_interpret_model.r

inputs:
  feat_in:
    type: File
    inputBinding:
      position: 1
      prefix: '--feat'
  original_feat:
    type: File
    inputBinding:
      position: 1
      prefix: '--origin_feat'
  label_in:
    type: File
    inputBinding:
      position: 1
      prefix: '--label'
  metadata_in:
    type: File?
    inputBinding:
      position: 1
      prefix: '--meta'
  model_tsv:
    type: File
    inputBinding:
      position: 1
      prefix: '--model'
  predictions:
    type: File
    inputBinding:
      position: 1
      prefix: '--pred'
  color_scheme:
    type: string?
    inputBinding:
      position: 2
      prefix: '--col_scheme'
  consensus_threshold:
    type: float?
    inputBinding:
      position: 2
      prefix: '--consens_thres'
  heatmap_type:
    type: string?
    inputBinding:
      position: 2
      prefix: '--heatmap_type'

arguments:
  - position: 3
    prefix: '--plot'
    valueFrom: $(inputs.feat_in.nameroot)_model_plots.pdf

outputs:
  model_plots:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_model_plots.pdf