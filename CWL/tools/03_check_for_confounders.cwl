cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
  
baseCommand: 03_check_for_confounders.r

inputs:
  label_in:
    type: File
    inputBinding:
      prefix: --label_in
      position: 2
  metadata_in:
    type: File
    inputBinding:
      position: 2
      prefix: --metadata_in
      
arguments:
  - prefix: --plot
    position: 2
    valueFrom: $(inputs.metadata_in.nameroot)_metaCheck.pdf
    position: 2

outputs:
  confounders_plot:
    type: File
    outputBinding:
      glob: $(inputs.metadata_in.nameroot)_metaCheck.pdf
