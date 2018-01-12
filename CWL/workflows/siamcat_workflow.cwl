class: Workflow
cwlVersion: v1.0
id: siamcat_workflow
label: siamcat_workflow
inputs:

  # Main input:
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

# for step check_associations:
  - id: alpha
    type: float?
    'sbg:exposed': true
  - id: detect_limit
    type: float?
    'sbg:exposed': true
  - id: max_show
    type: float?
    'sbg:exposed': true
  - id: min_fc
    type: float?
    'sbg:exposed': true
  - id: mult_test
    type: string?
    'sbg:exposed': true
  - id: plot_type
    type: string?
    'sbg:exposed': true
  - id: sort_by
    type: string?
    'sbg:exposed': true
  - id: assoc_col_scheme
    type: string?
    'sbg:exposed': true


  # for step filter_features:
  - id: rm_unmapped
    type: string?
    'sbg:exposed': true
  - id: filter_cutoff
    type: float?
    'sbg:exposed': true
  - id: filter_method
    type: string?
    'sbg:exposed': true
  - id: recomp_prop
    type: string?
    'sbg:exposed': true

  # for step normalize_features:
  - id: norm_method
    type: string?
    'sbg:exposed': true
  - id: norm_margin
    type: int?
    'sbg:exposed': true
  - id: sd_min_quantile
    type: float?
    'sbg:exposed': true
  - id: vector_norm
    type: int?
    'sbg:exposed': true

  # for step split_data:
  - id: num_folds
    type: int?
    'sbg:exposed': true
  - id: resample
    type: int?
    'sbg:exposed': true
  - id: inseparable
    type: string?
    'sbg:exposed': true

  # for step train_models:
  - id: train_method
    type: string?
    'sbg:exposed': true
  - id: min_nonzero_coeff
    type: int?
    'sbg:exposed': true
  - id: sel_criterion
    type: string?
    'sbg:exposed': true

  # for step evaluate_predictions
  - id: write_eval_results
    type: string?
    'sbg:exposed': true

  # for step interprete_model:
  - id: consensus_threshold
    type: float?
    'sbg:exposed': true
  - id: heatmap_type
    type: string?
    'sbg:exposed': true
  - id: interp_col_scheme
    type: string?
    'sbg:exposed': true

  # for steps split_data and train_models:
  - id: stratify
    type: string?
    'sbg:exposed': true

outputs:
  - id: association_plots_out
    outputSource:
      - check_associations/association_plots_out
    type: File
    'sbg:x': 1952.6685911746124
    'sbg:y': 725.331343086159
  - id: normalization_parameters
    outputSource:
      - normalize_features/normalization_parameters
    type: File?
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
  - id: model
    outputSource:
      - train_models/model
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
    run: ../tools/01_validate_data.cwl
    'sbg:x': 15.152003171656446
    'sbg:y': 213.78602506789176
  - id: filter_features
    in:
      - id: filter_cutoff
        source:
          - filter_cutoff
      - id: feat_in
        source:
          - validate_data/validated_feat
      - id: filter_method
        source:
          - filter_method
      - id: recomp_prop
        source:
          - recomp_prop
      - id: rm_unmapped
        source:
          - rm_unmapped
    out:
      - id: filtered_feat
    run: ../tools/04_filter_features.cwl
    'sbg:x': 245.98484980263123
    'sbg:y': 460.4359285994099
  - id: check_associations
    in:
      - id: alpha
        source:
          - alpha
      - id: detect_limit
        source:
          - detect_limit
      - id: feat_in
        source:
          - validate_data/validated_feat
      - id: label_in
        source:
          - validate_data/validated_label
      - id: max_show
        source:
          - max_show
      - id: min_fc
        source:
          - min_fc
      - id: mult_test
        source:
          - mult_test
      - id: plot_type
        source:
          - plot_type
      - id: sort_by
        source:
          - sort_by
      - id: assoc_col_scheme
        source:
          - assoc_col_scheme
    out:
      - id: association_plots_out
    run: ../tools/05_check_associations.cwl
    'sbg:x': 659.5949778238062
    'sbg:y': 725.4715164185877
  - id: normalize_features
    in:
      - id: feat_in
        source:
          - filter_features/filtered_feat
      - id: norm_method
        source:
          - norm_method
      - id: sd_min_quantile
        source:
          - sd_min_quantile
      - id: vector_norm
        source:
          - vector_norm
    out:
      - id: feat_out
      - id: normalization_parameters
    run: ../tools/06_normalize_features.cwl
    'sbg:x': 657.0314061482312
    'sbg:y': 465.6464459739973
  - id: split_data
    in:
      - id: inseparable
        source:
          - inseparable
      - id: label_in
        source:
          - validate_data/validated_label
      - id: metadata_in
        source:
          - validate_data/validated_metadata
      - id: num_folds
        source:
          - num_folds
      - id: resample
        source:
          - resample
      - id: stratify
        source:
          - stratify
    out:
      - id: test_sets_out
      - id: train_sets_out
    run: ../tools/08_split_data.cwl
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
      - id: model
        source:
          - train_models/model
      - id: test_sets
        source:
          - split_data/test_sets_out
    out:
      - id: predictions
    run: ../tools/10_make_predictions.cwl
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
      - id: train_method
        source:
          - train_method
      - id: min_nonzero_coeff
        source:
          - min_nonzero_coeff
      - id: sel_criterion
        source:
          - sel_criterion
      - id: stratify
        source:
          - stratify
      - id: train_sets
        source:
          - split_data/train_sets_out
    out:
      - id: model
    run: ../tools/09_train_models.cwl
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
        source:
          - write_eval_results
    out:
      - id: evaluation_plot
      - id: evaluation_results
    run: ../tools/11_evaluate_predictions.cwl
    'sbg:x': 1692.0694580078125
    'sbg:y': 473.4021848041717
  - id: interprete_model
    in:
      - id: consensus_threshold
        source:
          - consensus_threshold
      - id: feat_in
        source:
          - normalize_features/feat_out
      - id: heatmap_type
        source:
          - heatmap_type
      - id: label_in
        source:
          - validate_data/validated_label
      - id: metadata_in
        source:
          - validate_data/validated_metadata
      - id: model
        source:
          - train_models/model
      - id: original_feat
        source:
          - validate_data/validated_feat
      - id: predictions
        source:
          - make_predictions/predictions
      - id: interp_col_scheme
        source:
          - interp_col_scheme
    out:
      - id: model_plots
    run: ../tools/12_interprete_model.cwl
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
    run: ../tools/03_check_for_confounders.cwl
    'sbg:x': 648.4335802713812
    'sbg:y': -322.994492847952
requirements:
  MultipleInputFeatureRequirement: {}
