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
  - position: 0
    valueFrom: $(inputs.srcdir)/11_evaluate_predictions.r
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
