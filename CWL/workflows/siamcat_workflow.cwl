class: Workflow
cwlVersion: v1.0
id: siamcat_workflow
label: siamcat_workflow
inputs:
  - id: label_in
    type: File
    'sbg:x': 0
    'sbg:y': 114
  - id: feat_in
    type: File
    'sbg:x': 0
    'sbg:y': 221
  - id: metadata_in
    type: File?
    'sbg:x': 0
    'sbg:y': 7
outputs:
  - id: validated_label
    outputSource:
      - 02_select_samples/selected_label
    type: File
    'sbg:x': 699.37255859375
    'sbg:y': 0
  - id: plot
    outputSource:
      - 03_check_for_confounders/plot
    type: File
    'sbg:x': 858.3580932617188
    'sbg:y': 502.8577880859375
steps:
  - id: 01_validate_data
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
    'sbg:x': 150.046875
    'sbg:y': 100
  - id: 02_select_samples
    in:
      - id: feat_in
        source:
          - 01_validate_data/validated_feat
      - id: label_in
        source:
          - 01_validate_data/validated_label
      - id: metadata_in
        source:
          - 01_validate_data/validated_metadata
    out:
      - id: selected_feat
      - id: selected_label
      - id: selected_metadata
    run: >-
      /media/sf_Dokumente/MicrobiomeEMBL/siamcat_cwl/CWL/tools/02_select_samples.cwl
    'sbg:x': 381.709716796875
    'sbg:y': 100
  - id: 03_check_for_confounders
    in:
      - id: label_in
        source:
          - 02_select_samples/selected_label
      - id: metadata_in
        source:
          - 02_select_samples/selected_metadata
    out:
      - id: plot
    run: >-
      /media/sf_Dokumente/MicrobiomeEMBL/siamcat_cwl/CWL/tools/03_check_for_confounders.cwl
    'sbg:x': 601.7725628193322
    'sbg:y': 328.8736739720103
  - id: 04_filter_features
    in:
      - id: feat_in
        source:
          - 02_select_samples/selected_feat
    out:
      - id: filtered_feat
    run: >-
      /media/sf_Dokumente/MicrobiomeEMBL/siamcat_cwl/CWL/tools/04_filter_features.cwl
    'sbg:x': 698.3451901878254
    'sbg:y': 125.49263130664269
  - id: 05_check_associations
    in:
      - id: feat_in
        source:
          - 04_filter_features/filtered_feat
      - id: label_in
        source:
          - 02_select_samples/selected_label
    out:
      - id: association_plots_out
    run: >-
      /media/sf_Dokumente/MicrobiomeEMBL/siamcat_cwl/CWL/tools/05_check_associations.cwl
    'sbg:x': 875.2228110716128
    'sbg:y': 256.88105438867467
  - id: 06_normalize_features
    in:
      - id: feat_in
        source:
          - 04_filter_features/filtered_feat
    out:
      - id: feat_out
      - id: normalization_parameters_out
    run: >-
      /media/sf_Dokumente/MicrobiomeEMBL/siamcat_cwl/CWL/tools/06_normalize_features.cwl
    'sbg:x': 894.7427978515625
    'sbg:y': 25.95263188150811
  - id: 07_add_metadata_as_predictor
    in:
      - id: feat_in
        source:
          - 06_normalize_features/feat_out
      - id: metadata_in
        source:
          - 02_select_samples/selected_metadata
    out:
      - id: feat_out
    run: >-
      /media/sf_Dokumente/MicrobiomeEMBL/siamcat_cwl/CWL/tools/07_add_metadata_as_predictor.cwl
    'sbg:x': 1197.909912109375
    'sbg:y': 107
  - id: 08_split_data
    in: []
    out:
      - id: test_sets_out
      - id: train_sets_out
    run: /media/sf_Dokumente/MicrobiomeEMBL/siamcat_cwl/CWL/tools/08_split_data.cwl
    'sbg:x': 1200.3895263671875
    'sbg:y': 261
  - id: 09_train_models
    in:
      - id: train_sets
        source:
          - 08_split_data/train_sets_out
    out:
      - id: trained_model
      - id: trained_model_matrix
      - id: trained_model_rdata
    run: >-
      /media/sf_Dokumente/MicrobiomeEMBL/siamcat_cwl/CWL/tools/09_train_models.cwl
    'sbg:x': 1612.6225291005835
    'sbg:y': 204.79157740430864
  - id: 10_make_prediction
    in: []
    out:
      - id: predictions_out
    run: >-
      /media/sf_Dokumente/MicrobiomeEMBL/siamcat_cwl/CWL/tools/10_make_prediction.cwl
    'sbg:x': 1262.9136962890625
    'sbg:y': 436.7671203613281
  - id: 11_evaluate_predictions
    in: []
    out:
      - id: evaluation_plot
      - id: evaluation_results
    run: >-
      /media/sf_Dokumente/MicrobiomeEMBL/siamcat_cwl/CWL/tools/11_evaluate_predictions.cwl
    'sbg:x': 1515.3683987387767
    'sbg:y': 458.04526175098727
