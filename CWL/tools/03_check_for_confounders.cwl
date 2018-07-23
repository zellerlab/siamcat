cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
  
baseCommand: 03_check_for_confounders.r

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/siamcat:dev

inputs:
  label_in:
    type: File
    inputBinding:
      prefix: --label_in
      position: 1
  metadata_in:
    type: File
    inputBinding:
      position: 1
      prefix: --metadata_in
      
arguments:
  - prefix: --plot
    position: 2
    valueFrom: $(inputs.metadata_in.nameroot)_metaCheck.pdf

outputs:
  confounders_plot:
    type: File
    outputBinding:
      glob: $(inputs.metadata_in.nameroot)_metaCheck.pdf
