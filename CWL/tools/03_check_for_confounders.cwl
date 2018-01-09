cwlVersion: v1.0
class: CommandLineTool
baseCommand: Rscript
requirements:
  - class: InlineJavascriptRequirement

inputs:
  srcdir:
    type: string
    inputBinding:
      prefix: --srcdir
      position: 2
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
  - position: 0
    valueFrom: $(inputs.srcdir)/03_check_for_confounders.r
  - prefix: --plot
    position: 2
    valueFrom: $(inputs.metadata_in.nameroot)_metaCheck.pdf
    position: 2

outputs:
  plot:
    type: File
    outputBinding:
      glob: $(inputs.metadata_in.nameroot)_metaCheck.pdf
