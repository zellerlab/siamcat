class: CommandLineTool
cwlVersion: v1.0
baseCommand: Rscript
inputs:
  feat_in:
    type: File
    inputBinding:
      position: 2
      prefix: '--feat_in'
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
  srcdir:
    type: string
    inputBinding:
      position: 2
      prefix: '--srcdir'
outputs:
  validated_feat:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_valid.tsv
  validated_label:
    type: File
    outputBinding:
      glob: $(inputs.label_in.nameroot)_valid.tsv
  validated_metadata:
    type: File?
    outputBinding:
      glob: |
        ${
          if (inputs.metadata_in){
            return inputs.metadata_in.nameroot + "_valid.tsv";
          } else {
            return null;
          }
        }
arguments:
  - position: 0
    valueFrom: $(inputs.srcdir)/01_validate_data.r
  - position: 2
    prefix: '--feat_out'
    valueFrom: $(inputs.feat_in.nameroot)_valid.tsv
  - position: 2
    prefix: '--label_out'
    valueFrom: $(inputs.label_in.nameroot)_valid.tsv
  - position: 2
    valueFrom: |
      ${
        if (inputs.metadata_in){
          return [ "--metadata_out", inputs.metadata_in.nameroot + "_valid.tsv" ];
        } else {
          return null;
        }
      }
requirements:
  - class: InlineJavascriptRequirement