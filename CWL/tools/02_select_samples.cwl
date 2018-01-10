class: CommandLineTool
cwlVersion: v1.0

requirements:
  InlineJavascriptRequirement: {}

baseCommand: 02_select_samples.r

inputs:
  allowed_range:
    type: string?
    inputBinding:
      position: 2
      prefix: '--allowed_range'
  allowed_set:
    type: string?
    inputBinding:
      position: 2
      prefix: '--allowed_set'
  feat_in:
    type: File
    inputBinding:
      position: 2
      prefix: '--feat_in'
  filter_var:
    type: string
    inputBinding:
      position: 2
      prefix: '--filter_var'
  label_in:
    type: File
    inputBinding:
      position: 2
      prefix: '--label_in'
  metadata_in:
    type: File?
    inputBinding:
      position: 2
      prefix: '--metadata_in'

arguments:
  - position: 2
    prefix: '--feat_out'
    valueFrom: $(inputs.feat_in.nameroot)_select.tsv
  - position: 2
    prefix: '--label_out'
    valueFrom: $(inputs.label_in.nameroot)_select_.tsv
  - position: 2
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