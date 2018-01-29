cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
 
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/siamcat:dev

baseCommand: 04_filter_features.r

inputs:
  feat_in:
    type: File
    inputBinding:
      prefix: --feat_in
      position: 1
  filter_method:
    type: string?
    inputBinding:
      prefix: --method
      position: 2
  filter_cutoff:
    type: float?
    inputBinding:
      prefix: --cutoff
      position: 2
  rm_unmapped:
    doc: TRUE or FALSE is allowed
    type: string?
    inputBinding:
      prefix: --rm_unmapped
      position: 2
  recomp_prop:
    doc: TRUE or FALSE is allowed
    type: string?
    inputBinding:
      prefix: --recomp_prop
      position: 2

arguments:
  - prefix: --feat_out
    valueFrom: $(inputs.feat_in.nameroot)_filtered.tsv
    position: 3

outputs:
  filtered_feat:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_filtered.tsv
