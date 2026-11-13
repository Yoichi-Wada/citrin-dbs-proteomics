LoF/LoF frozen ensemble scoring (output directory revised 2026-11-04)

1. Open 20261029_rebuttal.Rproj in RStudio.
2. Save the .R script into the project script/ folder.
3. Run:
   source("script/citrin_SLC25A13_LoFLoF_frozen_ensemble_scoring_20261104.R")

Only output location changed from the previous version.
Results will be written to the project root's existing result/ directory:
  result/LoFLoF_frozen_ensemble/
The subfolder will be created automatically.
Model loading and scoring implementation are unchanged.
Requires the frozen final_training_model_NOT_VALIDATED.rds created by the ensemble training.
