Citrin 10-member final model reproducibility (2026-11-04)
========================================================

1. Open the RStudio project: 20261029_rebuttal.Rproj.
2. Place citrin_final_ensemble_reproducibility_20261104.R in script/.
3. Confirm the original final-training model exists at:
   citrin_nested_ensemble_v1_1/final_training_model_NOT_VALIDATED.rds
   (or set MODEL_RDS_OVERRIDE at the start of the script).
4. Confirm raw_data.xlsx and sample_info.csv are in the usual project data/ folder.
5. Run:
   source("script/citrin_final_ensemble_reproducibility_20261104.R",encoding="UTF-8")

Output (created automatically):
   result/FinalModel_reproducibility/

Files include:
   final_training_model_FROZEN_NOT_VALIDATED.rds (exact-byte archived copy)
   consensus30_frozen_protein_list.csv
   diverse30_per_member_protein_lists.csv
   member_glmnet_coefficients_ALL_10.csv
   member_model_parameters.csv
   member_preprocessing_center_scale_ALL.csv
   frozen_ComBat_batch_shifts_ALL.csv (only for a ComBat-trained model)
   final_training_inner_thresholds.csv
   scoring_example_imputation_quantiles.csv
   scoring_example_member_calculations.csv
   scoring_example_feature_contributions.csv
   scoring_example_ensemble_summary.csv
   scoring_rule_methods_and_caveats.txt
   model_manifest.csv
   R_sessionInfo.txt

No classifier fitting or parameter selection occurs. Two WORKED EXAMPLE
specimens (SM_1 and WS_1) belonged to the original training dataset. These
are in-sample implementation examples, NOT independent validation and NOT
OOF predictions. The model RDS and original inputs are not modified.

Dependencies: R with openxlsx and glmnet installed.
If the frozen RDS is missing, the script stops instead of creating a new model.
