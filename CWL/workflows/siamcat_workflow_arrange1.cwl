class: Workflow
cwlVersion: v1.0
id: siamcat_workflow
label: siamcat_workflow
inputs:
  - id: label_in
    type: File
    'sbg:x': -211.46906743520267
    'sbg:y': 215.17076736908456
  - id: feat_in
    type: File
    'sbg:x': -213.1508026123047
    'sbg:y': 354.48820821737314
  - id: metadata_in
    type: File?
    'sbg:x': -213.43544332479505
    'sbg:y': 65.19374084472656
outputs:
  - id: association_plots_out
    outputSource:
      - check_associations/association_plots_out
    type: File
    'sbg:x': 1920.3959628117993
    'sbg:y': -265.89943341310783
  - id: normalization_parameters_out
    outputSource:
      - normalize_features/normalization_parameters_out
    type: File
    'sbg:x': 1104.0839124465738
    'sbg:y': -203.8899465928698
  - id: evaluation_results
    outputSource:
      - evaluate_predictions/evaluation_results
    type: File?
    'sbg:x': 1931.515819217699
    'sbg:y': 423.6074051792846
  - id: evaluation_plot
    outputSource:
      - evaluate_predictions/evaluation_plot
    type: File
    'sbg:x': 1936.9710693359375
    'sbg:y': 560.4991510912977
  - id: model_plots
    outputSource:
      - interprete_model/model_plots
    type: File
    'sbg:x': 1924.0077850204946
    'sbg:y': -32.05451321281125
  - id: confounders_plot
    outputSource:
      - check_for_confounders/confounders_plot
    type: File
    'sbg:x': 1939.0229602762402
    'sbg:y': 771.1966557660979
  - id: model_tsv
    outputSource:
      - train_models/model_tsv
    type: File
    'sbg:x': 1105.4731946731363
    'sbg:y': -87.36692340320536
steps:
  - id: validate_data
    in:
      - id: feat_in
        source:
          - feat_in
      - id: label_in
        source:
          - label_in
      - id: metadata_in
        source:
          - metadata_in
    out:
      - id: validated_feat
      - id: validated_label
      - id: validated_metadata
    run: >-
      /media/sf_Dokumente/MicrobiomeEMBL/siamcat_cwl/CWL/tools/01_validate_data.cwl
    'sbg:x': 15.152003171656446
    'sbg:y': 213.78602506789176
  - id: filter_features
    in:
      - id: feat_in
        source:
          - select_samples/selected_feat
          - validate_data/validated_feat
    out:
      - id: filtered_feat
    run: >-
      /media/sf_Dokumente/MicrobiomeEMBL/siamcat_cwl/CWL/tools/04_filter_features.cwl
    'sbg:x': 249.00100933052406
    'sbg:y': 12.314643731955186
  - id: check_associations
    in:
      - id: feat_in
        source:
          - filter_features/filtered_feat
      - id: label_in
        source:
          - select_samples/selected_label
          - validate_data/validated_label
    out:
      - id: association_plots_out
    run: >-
      /media/sf_Dokumente/MicrobiomeEMBL/siamcat_cwl/CWL/tools/05_check_associations.cwl
    'sbg:x': 758.5679099044457
    'sbg:y': -262.3459206756455
  - id: normalize_features
    in:
      - id: feat_in
        source:
          - filter_features/filtered_feat
    out:
      - id: feat_out
      - id: normalization_parameters_out
    run: >-
      /media/sf_Dokumente/MicrobiomeEMBL/siamcat_cwl/CWL/tools/06_normalize_features.cwl
    'sbg:x': 546.7723766116638
    'sbg:y': 2.6954385965526964
  - id: split_data
    in:
      - id: label_in
        source:
          - validate_data/validated_label
      - id: metadata_in
        source:
          - validate_data/validated_metadata
    out:
      - id: test_sets_out
      - id: train_sets_out
    run: /media/sf_Dokumente/MicrobiomeEMBL/siamcat_cwl/CWL/tools/08_split_data.cwl
    'sbg:x': 617.3671975721693
    'sbg:y': 497.5342735119465
  - id: make_predictions
    in:
      - id: feat_in
        source:
          - normalize_features/feat_out
      - id: label_in
        source:
          - validate_data/validated_label
      - id: model_matrix
        source:
          - train_models/model_matrix
      - id: model_rdata
        source:
          - train_models/model_rdata
      - id: test_sets
        source:
          - split_data/test_sets_out
    out:
      - id: predictions
    run: >-
      /media/sf_Dokumente/MicrobiomeEMBL/siamcat_cwl/CWL/tools/10_make_predictions.cwl
    'sbg:x': 1480.8948640695066
    'sbg:y': 239.193603515625
  - id: train_models
    in:
      - id: feat_in
        source:
          - normalize_features/feat_out
      - id: label_in
        source:
          - validate_data/validated_label
      - id: train_sets
        source:
          - split_data/train_sets_out
    out:
      - id: model_matrix
      - id: model_rdata
      - id: model_tsv
    run: >-
      /media/sf_Dokumente/MicrobiomeEMBL/siamcat_cwl/CWL/tools/09_train_models.cwl
    'sbg:x': 885.5039230004554
    'sbg:y': 187.50314745710568
  - id: evaluate_predictions
    in:
      - id: label_in
        source:
          - validate_data/validated_label
      - id: predictions
        source:
          - make_predictions/predictions
      - id: write_eval_results
        default: true
    out:
      - id: evaluation_plot
      - id: evaluation_results
    run: >-
      /media/sf_Dokumente/MicrobiomeEMBL/siamcat_cwl/CWL/tools/11_evaluate_predictions.cwl
    'sbg:x': 1742.0773034613228
    'sbg:y': 478.7590748329334
  - id: interprete_model
    in:
      - id: feat_in
        source:
          - normalize_features/feat_out
      - id: label_in
        source:
          - validate_data/validated_label
      - id: metadata_in
        source:
          - validate_data/validated_metadata
      - id: model_tsv
        source:
          - train_models/model_tsv
      - id: original_feat
        source:
          - validate_data/validated_feat
      - id: predictions
        source:
          - make_predictions/predictions
    out:
      - id: model_plots
    run: >-
      /media/sf_Dokumente/MicrobiomeEMBL/siamcat_cwl/CWL/tools/12_interprete_model.cwl
    'sbg:x': 1741.5103701741173
    'sbg:y': -29.675089878030956
  - id: check_for_confounders
    in:
      - id: label_in
        source:
          - validate_data/validated_label
      - id: metadata_in
        source:
          - validate_data/validated_metadata
    out:
      - id: confounders_plot
    run: >-
      /media/sf_Dokumente/MicrobiomeEMBL/siamcat_cwl/CWL/tools/03_check_for_confounders.cwl
    'sbg:x': 490.946993293249
    'sbg:y': 753.4881214090526
