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

baseCommand: 10_2_merge_data.r

inputs:
  input_predictions:
    doc: |
      array of the predictions files
      to be merge
    type:
      type: array
      items: File

arguments:
  - prefix: --merged_pred
    valueFrom: |
      ${ inputs.input_predictions[0].nameroot.replace(/\d\d\d/i,".tsv") }
    # name for merged output file
  - prefix: --pred_prefix
    valueFrom: "*"
    position: 1
    # Prefix of files, can be set to "*"
    # since no other files are assumed 
    # to exist in the runtime working dir

outputs:
  merged_predictions:
    type: File
    outputBinding:
      glob: ${ inputs.input_predictions[0].nameroot.replace(/\d\d\d/i,".tsv") }