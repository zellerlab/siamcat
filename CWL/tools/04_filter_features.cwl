cwlVersion: v1.0
class: CommandLineTool
baseCommand: Rscript
requirements:
  - class: InlineJavascriptRequirement

inputs:
  feat_in:
    type: File
    inputBinding:
      prefix: --feat_in
      position: 2
  filt_method:
    type: string
    inputBinding:
      prefix: --method
      position: 2
  cutoff:
    type: float
    inputBinding:
      prefix: --cutoff
      position: 2
  rm_unmapped:
    type: boolean
    inputBinding:
      prefix: --rm_unmapped
      position: 2
      valueFrom: $(self.toString())
  recomp_prop:
    type: boolean
    inputBinding:
      prefix: --recomp_prop
      position: 2
      valueFrom: $(self.toString())

arguments:
  - position: 0
    valueFrom: $(inputs.srcdir)/04_filter_features.r
  - prefix: --feat_out
    valueFrom: $(inputs.feat_in.nameroot)_filtered.tsv
    position: 2

outputs:
  filtered_feat:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_filtered.tsv
