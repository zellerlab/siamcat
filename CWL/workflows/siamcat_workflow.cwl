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
    'sbg:x': 1952.6685911746124
    'sbg:y': 725.331343086159
  - id: normalization_parameters_out
    outputSource:
      - normalize_features/normalization_parameters_out
    type: File
    'sbg:x': 1170.382378133225
    'sbg:y': 588.6885709128314
  - id: evaluation_results
    outputSource:
      - evaluate_predictions/evaluation_results
    type: File?
    'sbg:x': 1952.359835561269
    'sbg:y': 348.0480077107312
  - id: evaluation_plot
    outputSource:
      - evaluate_predictions/evaluation_plot
    type: File
    'sbg:x': 1949.9984400431752
    'sbg:y': 559.1964096073256
  - id: model_plots
    outputSource:
      - interprete_model/model_plots
    type: File
    'sbg:x': 1943.5490857442437
    'sbg:y': -29.44900868598708
  - id: confounders_plot
    outputSource:
      - check_for_confounders/confounders_plot
    type: File
    'sbg:x': 1941.2224255879937
    'sbg:y': -322.2702631625083
  - id: model_tsv
    outputSource:
      - train_models/model_tsv
    type: File
    'sbg:x': 1183.3419863384065
    'sbg:y': -202.1716046648809
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
    'sbg:x': 245.98484980263123
    'sbg:y': 460.4359285994099
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
    'sbg:x': 659.5949778238062
    'sbg:y': 725.4715164185877
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
    'sbg:x': 657.0314061482312
    'sbg:y': 465.6464459739973
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
    'sbg:x': 645.8034128823999
    'sbg:y': -34.249268468373465
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
    'sbg:x': 1363.3388679502434
    'sbg:y': 206.673560587478
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
    'sbg:x': 951.9482287088812
    'sbg:y': 176.0787371317328
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
    'sbg:x': 1692.0694580078125
    'sbg:y': 473.4021848041717
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
    'sbg:x': 1699.3190236410376
    'sbg:y': -29.37906087241761
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
    'sbg:x': 648.4335802713812
    'sbg:y': -322.994492847952
