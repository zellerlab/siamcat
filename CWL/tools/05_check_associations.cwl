cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/siamcat:0.3.1

baseCommand: 05_check_associations.r

inputs:
  feat_in:
    type: File
    inputBinding:
      position: 1
      prefix: --feat_in
  label_in:
    type: File
    inputBinding:
      position: 1
      prefix: --label_in
  mult_test:
    type: string?
    inputBinding:
      position: 2
      prefix: --mult_test
  alpha:
    type: float?
    inputBinding:
      position: 2
      prefix: --alpha
  min_fc:
    type: float?
    inputBinding:
      position: 2
      prefix: --min_fc
  max_show:
    type: float?
    inputBinding:
      position: 2
      prefix: --max_show
  detect_limit:
    type: float?
    inputBinding:
      position: 2
      prefix: --detect_limit
  sort_by:
    type: string?
    inputBinding:
      position: 2
      prefix: --sort_by
  assoc_col_scheme:
    type: string?
    inputBinding:
      position: 2
      prefix: --col_scheme
  plot_type:
    type: string?
    inputBinding:
      position: 2
      prefix: --plot_type

arguments:
  - position: 3
    prefix: --plot
    valueFrom: $(inputs.feat_in.nameroot)_assocToLabels.pdf

outputs:
  association_plots_out:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_assocToLabels.pdf
