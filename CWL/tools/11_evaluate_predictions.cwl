cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}

baseCommand: 11_evaluate_predictions.r

inputs:
  label_in:
    type: File
    inputBinding:
      prefix: --label
      position: 2
  predictions:
    type: File
    inputBinding:
      prefix: --pred
      position: 2
  write_eval_results:
    type: boolean?
    inputBinding:
      prefix: --write_eval_results
      valueFrom: $(self.toString())
      position: 2

arguments:
  - position: 2
    prefix: --plot
    valueFrom:   $(inputs.feat_in.nameroot)_evalPlots.pdf
  - position: 2
    prefix: --output_results
    valueFrom:   $(inputs.feat_in.nameroot)_evalResults.txt

outputs:
  evaluation_plot:
    type: File
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_evalPlot.pdf
  evaluation_results:
    type: File?
    outputBinding:
      glob: $(inputs.feat_in.nameroot)_evalResults.txt
