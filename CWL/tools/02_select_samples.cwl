class: CommandLineTool
cwlVersion: v1.0

requirements:
  InlineJavascriptRequirement: {}

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/siamcat:dev

baseCommand: 02_select_samples.r

inputs:
  feat_in:
    type: File
    inputBinding:
      position: 1
      prefix: '--feat_in'
  label_in:
    type: File
    inputBinding:
      position: 1
      prefix: '--label_in'
  metadata_in:
    type: File?
    inputBinding:
      position: 1
      prefix: '--metadata_in'
  filter_variable:
    type: string
    inputBinding:
      position: 2
      prefix: '--filter_var'
  exclusive_parameters: 
    # either allowed_range or allowed_set can be provided
    # if both are provided allowed_set is ignored
    allowed_range:
      type: string
      inputBinding:
        position: 2
        prefix: '--allowed_range'
    allowed_set:
      type: string
      inputBinding:
        position: 2
        prefix: '--allowed_set'

arguments:
  - position: 3
    prefix: '--feat_out'
    valueFrom: $(inputs.feat_in.nameroot)_select.tsv
  - position: 3
    prefix: '--label_out'
    valueFrom: $(inputs.label_in.nameroot)_select_.tsv
  - position: 3
    valueFrom: |
      ${
        if (inputs.metadata_in){
          return [ "--metadata_out", inputs.metadata_in.nameroot + "_select.tsv" ];
        } else {
          return null;
        }
      }
requirements:
  - class: InlineJavascriptRequirement

outputs:
  selected_feat:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_select.tsv
  selected_label:
    type: File
    outputBinding:
      glob: $(inputs.label_in.nameroot)_select.tsv
  selected_metadata:
    type: File?
    outputBinding:
      glob: |
        ${
          if (inputs.metadata_in){
            return inputs.metadata_in.nameroot + "_select.tsv";
          } else {
            return null;
          }
        }