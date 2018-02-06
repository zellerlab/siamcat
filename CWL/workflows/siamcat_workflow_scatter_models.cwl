class: Workflow
cwlVersion: v1.0
id: siamcat_workflow
label: siamcat_workflow

requirements:
  MultipleInputFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}

inputs:
  # Main input:
  label_in:
    type: File
    'sbg:x': -211.46906743520267
    'sbg:y': 215.17076736908456
  feat_in:
    type: File
    'sbg:x': -213.1508026123047
    'sbg:y': 354.48820821737314
  metadata_in:
    type: File?
    'sbg:x': -213.43544332479505
    'sbg:y': 65.19374084472656

# for step check_associations:
  alpha:
    type: float?
    'sbg:exposed': true
  detect_limit:
    type: float?
    'sbg:exposed': true
  max_show:
    type: float?
    'sbg:exposed': true
  min_fc:
    type: float?
    'sbg:exposed': true
  mult_test:
    type: string?
    'sbg:exposed': true
  plot_type:
    type: string?
    'sbg:exposed': true
  sort_by:
    type: string?
    'sbg:exposed': true
  assoc_col_scheme:
    type: string?
    'sbg:exposed': true


  # for step filter_features:
  rm_unmapped:
    type: string?
    'sbg:exposed': true
  filter_cutoff:
    type: float?
    'sbg:exposed': true
  filter_method:
    type: string?
    'sbg:exposed': true
  recomp_prop:
    type: string?
    'sbg:exposed': true

  # for step normalize_features:
  norm_method:
    type: string?
    'sbg:exposed': true
  norm_margin:
    type: int?
    'sbg:exposed': true
  sd_min_quantile:
    type: float?
    'sbg:exposed': true
  vector_norm:
    type: int?
    'sbg:exposed': true

  # for step split_data:
  num_folds:
    type: int?
    'sbg:exposed': true
  resample:
    type: int?
    'sbg:exposed': true
  inseparable:
    type: string?
    'sbg:exposed': true

  # for step train_models:
  train_method:
    type: string?
    'sbg:exposed': true
  min_nonzero_coeff:
    type: int?
    'sbg:exposed': true
  sel_criterion:
    type: string?
    'sbg:exposed': true

  # for step evaluate_predictions
  write_eval_results:
    type: string?
    'sbg:exposed': true

  # for step interprete_model:
  consensus_threshold:
    type: float?
    'sbg:exposed': true
  heatmap_type:
    type: string?
    'sbg:exposed': true
  interp_col_scheme:
    type: string?
    'sbg:exposed': true

  # for steps split_data and train_models:
  stratify:
    type: string?
    'sbg:exposed': true


steps:
  validate_data:
    in:
      feat_in:
        source: feat_in
      label_in:
        source: label_in
      metadata_in:
        source: metadata_in
    out:
      - validated_feat
      - validated_label
      - validated_metadata
    run: ../tools/01_validate_data.cwl
    'sbg:x': 15.152003171656446
    'sbg:y': 213.78602506789176

  check_for_confounders:
    in:
      label_in:
        source: validate_data/validated_label
      metadata_in:
        source: validate_data/validated_metadata
    out:
      - confounders_plot
    run: ../tools/03_check_for_confounders.cwl
    'sbg:x': 648.4335802713812
    'sbg:y': -322.994492847952

  filter_features:
    in:
      filter_cutoff:
        source: filter_cutoff
      feat_in:
        source: validate_data/validated_feat
      filter_method:
        source: filter_method
      recomp_prop:
        source: recomp_prop
      rm_unmapped:
        source: rm_unmapped
    out:
      - filtered_feat
    run: ../tools/04_filter_features.cwl
    'sbg:x': 245.98484980263123
    'sbg:y': 460.4359285994099

  check_associations:
    in:
      alpha:
        source: alpha
      detect_limit:
        source: detect_limit
      feat_in:
        source: validate_data/validated_feat
      label_in:
        source: validate_data/validated_label
      max_show:
        source: max_show
      min_fc:
        source: min_fc
      mult_test:
        source: mult_test
      plot_type:
        source: plot_type
      sort_by:
        source: sort_by
      assoc_col_scheme:
        source: assoc_col_scheme
    out:
      - association_plots_out
    run: ../tools/05_check_associations.cwl
    'sbg:x': 659.5949778238062
    'sbg:y': 725.4715164185877

  normalize_features:
    in:
      feat_in:
        source: filter_features/filtered_feat
      norm_method:
        source: norm_method
      sd_min_quantile:
        source: sd_min_quantile
      vector_norm:
        source: vector_norm
    out:
      - feat_out
      - normalization_parameters
    run: ../tools/06_normalize_features.cwl
    'sbg:x': 657.0314061482312
    'sbg:y': 465.6464459739973
    
  split_data:
    in:
      inseparable:
        source: inseparable
      label_in:
        source: validate_data/validated_label
      metadata_in:
        source: validate_data/validated_metadata
      num_folds:
        source: num_folds
      resample:
        source: resample
      stratify:
        source: stratify
      subdivide_train_set:
        valueFrom: "TRUE"
    out:
      - test_sets_out
      - train_sets_out
    run: ../tools/08_split_data.cwl
    'sbg:x': 645.8034128823999
    'sbg:y': -34.249268468373465

  train_models:
    scatter: [train_sets]
    scatterMethod: 'dotproduct'
    in:
      feat_in:
        source: normalize_features/feat_out
      label_in:
        source: validate_data/validated_label
      train_method:
        source: train_method
      min_nonzero_coeff:
        source: min_nonzero_coeff
      sel_criterion:
        source: sel_criterion
      stratify:
        source: stratify
      train_sets:
        source: split_data/train_sets_out
    out:
      - model
    run: ../tools/09_train_models.cwl
    'sbg:x': 951.9482287088812
    'sbg:y': 176.0787371317328

  join_models:
    in:
      input_mlr_models:
        source: train_models/model
    out:
      - joined_model
    run: ../tools/09.2_join_models.cwl
    'sbg:x': 1400.0
    'sbg:y': 206.673560587478

  make_predictions:
    in:
      feat_in:
        source: normalize_features/feat_out
      label_in:
        source: validate_data/validated_label
      model:
        source: join_models/joined_model
      test_sets:
        source: split_data/test_sets_out
    out:
      - predictions
    run: ../tools/10_make_predictions.cwl
    'sbg:x': 1363.3388679502434
    'sbg:y': 206.673560587478

  evaluate_predictions:
    in:
      label_in:
        source: validate_data/validated_label
      predictions:
        source: make_predictions/predictions
      write_eval_results:
        source: write_eval_results
    out:
      - evaluation_plot
      - evaluation_results
    run: ../tools/11_evaluate_predictions.cwl
    'sbg:x': 1692.0694580078125
    'sbg:y': 473.4021848041717

  interprete_model:
    in:
      consensus_threshold:
        source: consensus_threshold
      feat_in:
        source: normalize_features/feat_out
      heatmap_type:
        source: heatmap_type
      label_in:
        source: validate_data/validated_label
      metadata_in:
        source: validate_data/validated_metadata
      model:
        source: join_models/joined_model
      original_feat:
        source: validate_data/validated_feat
      predictions:
        source: make_predictions/predictions
      interp_col_scheme:
        source: interp_col_scheme
    out:
      - model_plots
    run: ../tools/12_interprete_model.cwl
    'sbg:x': 1699.3190236410376
    'sbg:y': -29.37906087241761

outputs:
  association_plots_out:
    outputSource: check_associations/association_plots_out
    type: File
    'sbg:x': 1952.6685911746124
    'sbg:y': 725.331343086159
  normalization_parameters:
    outputSource: normalize_features/normalization_parameters
    type: File?
    'sbg:x': 1170.382378133225
    'sbg:y': 588.6885709128314
  evaluation_results:
    outputSource: evaluate_predictions/evaluation_results
    type: File?
    'sbg:x': 1952.359835561269
    'sbg:y': 348.0480077107312
  evaluation_plot:
    outputSource: evaluate_predictions/evaluation_plot
    type: File
    'sbg:x': 1949.9984400431752
    'sbg:y': 559.1964096073256
  model_plots:
    outputSource: interprete_model/model_plots
    type: File
    'sbg:x': 1943.5490857442437
    'sbg:y': -29.44900868598708
  confounders_plot:
    outputSource: check_for_confounders/confounders_plot
    type: File
    'sbg:x': 1941.2224255879937
    'sbg:y': -322.2702631625083
  model:
    outputSource: join_models/joined_model
    type: File
    'sbg:x': 1183.3419863384065
    'sbg:y': -202.1716046648809