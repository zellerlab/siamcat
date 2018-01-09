class: CommandLineTool
cwlVersion: v1.0
baseCommand: Rscript
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
  srcdir:
    type: string
    inputBinding:
      position: 2
      prefix: '--srcdir'
outputs:
  model_plots:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_model_plots.pdf
arguments:
  - position: 0
    valueFrom: $(inputs.srcdir)/12_interpret_model.r
  - position: 2
    prefix: '--plot'
    valueFrom: $(inputs.feat_in.nameroot)_model_plots.pdf
requirements:
  - class: InlineJavascriptRequirement