class: CommandLineTool
cwlVersion: v1.0

requirements:
  InlineJavascriptRequirement: {}

baseCommand: 12_interpret_model.r

inputs:
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
  feat_in:
    type: File
    inputBinding:
      position: 2
      prefix: '--feat'
  heatmap_type:
    type: string?
    inputBinding:
      position: 2
      prefix: '--heatmap_type'
  label_in:
    type: File
    inputBinding:
      position: 2
      prefix: '--label'
  metadata_in:
    type: File?
    inputBinding:
      position: 2
      prefix: '--meta'
  model_tsv:
    type: File
    inputBinding:
      position: 2
      prefix: '--model'
  original_feat:
    type: File
    inputBinding:
      position: 2
      prefix: '--origin_feat'
  predictions:
    type: File
    inputBinding:
      position: 2
      prefix: '--pred'

arguments:
  - position: 2
    prefix: '--plot'
    valueFrom: $(inputs.feat_in.nameroot)_model_plots.pdf

outputs:
  model_plots:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_model_plots.pdf