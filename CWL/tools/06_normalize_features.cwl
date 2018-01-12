cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}

baseCommand: 06_normalize_features.r

inputs:
  feat_in:
    type: File
    inputBinding:
      prefix: --feat_in
      position: 1
  norm_method:
    type: string?
    inputBinding:
      prefix: --method
      position: 2
  log_n0:
    type: float?
    inputBinding:
      prefix: --log_n0
      position: 2
  sd_min_quantile:
    type: float?
    inputBinding:
      prefix: --sd_min_quantile
      position: 2
  norm_sample:
    doc: TRUE or FALSE is allowed
    type: string?
    inputBinding:
      prefix: --norm_sample
      position: 2
  norm_global:
    doc: TRUE or FALSE is allowed
    type: string?
    inputBinding:
      prefix: --norm_global
      position: 2
  vector_norm:
    type: int?
    inputBinding:
      prefix: --vector_norm
      position: 2
  norm_feature:
    doc: TRUE or FALSE is allowed
    type: string?
    inputBinding:
      prefix: --norm_feature
      position: 2
      
arguments:
  - position: 3
    prefix:  --feat_out
    valueFrom: $(inputs.feat_in.nameroot)_norm.tsv
  - position: 3
    prefix:  --param_out
    valueFrom: $(inputs.feat_in.nameroot)_normParam.txt
    
outputs:
  feat_out:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_norm.tsv
  normalization_parameters:
    type: File?
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_normParam.txt
