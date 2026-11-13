# Table S19. Application of the frozen final classifier to the 16 LoF/LoF specimens.

| Specimen | Genotype | Acquisition run | Run contained training specimens | SLC25A13 detected | Quantified proteins (of 5,199) | Consensus-30 probability | Consensus-30 call, full inner-OOF recall threshold (0.034) | Consensus-30 call, inner-OOF specificity threshold (0.210) | Diverse-30 probability | Diverse-30 call, full inner-OOF recall threshold (0.090) | Diverse-30 call, inner-OOF specificity threshold (0.517) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| SS_1 | LoF/LoF | Run1 | No (extrapolation) | No | 4367 | 0.058 | positive | negative | 0.145 | positive | negative |
| SS_2 | LoF/LoF | Run1 | No (extrapolation) | No | 4436 | 0.067 | positive | negative | 0.171 | positive | negative |
| SS_3 | LoF/LoF | Run1 | No (extrapolation) | No | 4642 | 0.068 | positive | negative | 0.150 | positive | negative |
| SS_4 | LoF/LoF | Run1 | No (extrapolation) | No | 4781 | 0.010 | negative | negative | 0.028 | negative | negative |
| SS_5 | LoF/LoF | Run1 | No (extrapolation) | No | 4630 | 0.021 | negative | negative | 0.112 | positive | negative |
| SS_6 | LoF/LoF | Run1 | No (extrapolation) | No | 4879 | 0.279 | positive | positive | 0.432 | positive | negative |
| SS_7 | LoF/LoF | Run1 | No (extrapolation) | No | 4762 | 0.035 | positive | negative | 0.096 | positive | negative |
| SS_8 | LoF/LoF | Run1 | No (extrapolation) | No | 4788 | 0.019 | negative | negative | 0.133 | positive | negative |
| SS_9 | LoF/LoF | Run2 | Yes | No | 4508 | 0.016 | negative | negative | 0.075 | negative | negative |
| SS_10 | LoF/LoF | Run2 | Yes | No | 4396 | 0.006 | negative | negative | 0.010 | negative | negative |
| SS_11 | LoF/LoF | Run2 | Yes | No | 4953 | 0.046 | positive | negative | 0.171 | positive | negative |
| SS_12 | LoF/LoF | Run2 | Yes | No | 4357 | 0.042 | positive | negative | 0.042 | negative | negative |
| SS_13 | LoF/LoF | Run3 | Yes | No | 4862 | 0.108 | positive | negative | 0.174 | positive | negative |
| SS_14 | LoF/LoF | Run3 | Yes | No | 4957 | 0.055 | positive | negative | 0.158 | positive | negative |
| SS_15 | LoF/LoF | Run3 | Yes | No | 4623 | 0.202 | positive | negative | 0.166 | positive | negative |
| SS_16 | LoF/LoF | Run3 | Yes | No | 4451 | 0.090 | positive | negative | 0.052 | negative | negative |

**Summary**

| Strategy | Group | Median probability (range) | Positive at full-recall threshold | Positive at specificity threshold |
|---|---|---|---|---|
| Consensus-30 | Reference: 37 controls, out-of-fold (nested CV) | 0.030 (0.005–0.569) | — | — |
| Consensus-30 | Reference: 7 LoF/Mis, out-of-fold (nested CV) | 0.131 (0.029–0.814) | — | — |
| Consensus-30 | LoF/LoF, all runs (n = 16), frozen model | 0.051 (0.006–0.279) | 11/16 | 1/16 |
| Consensus-30 | LoF/LoF, Run1 (n = 8), frozen model | 0.047 (0.010–0.279) | 5/8 | 1/8 |
| Consensus-30 | LoF/LoF, Run2 (n = 4), frozen model | 0.029 (0.006–0.046) | 2/4 | 0/4 |
| Consensus-30 | LoF/LoF, Run3 (n = 4), frozen model | 0.099 (0.055–0.202) | 4/4 | 0/4 |
| Diverse-30 | Reference: 37 controls, out-of-fold (nested CV) | 0.048 (0.006–0.508) | — | — |
| Diverse-30 | Reference: 7 LoF/Mis, out-of-fold (nested CV) | 0.152 (0.087–0.706) | — | — |
| Diverse-30 | LoF/LoF, all runs (n = 16), frozen model | 0.139 (0.010–0.432) | 11/16 | 0/16 |
| Diverse-30 | LoF/LoF, Run1 (n = 8), frozen model | 0.139 (0.028–0.432) | 7/8 | 0/8 |
| Diverse-30 | LoF/LoF, Run2 (n = 4), frozen model | 0.059 (0.010–0.171) | 1/4 | 0/4 |
| Diverse-30 | LoF/LoF, Run3 (n = 4), frozen model | 0.162 (0.052–0.174) | 3/4 | 0/4 |

The final model (consensus 30-protein panel with ten elastic-net members; secondary diverse-panel strategy) was trained on the 44 non-LoF/LoF specimens and applied without modification. Columns give the acquisition run, whether that run contained training specimens (Run 1 did not; scores for Run 1 are extrapolations), SLC25A13 detection (undetectable in all 16), the number of quantified proteins among the 5,199 eligible proteins, and for each strategy the ensemble-mean predicted probability and the calls at the two thresholds derived from the inner family-leave-one-family-out validation of the 44-specimen training cohort (full inner-OOF recall: 0.034 consensus / 0.090 diverse; inner-OOF specificity ≥ 0.90: 0.210 / 0.517).
For reference, the summary block gives the out-of-fold probabilities of the 37 controls and 7 LoF/Mis specimens in the nested cross-validation (medians 0.030 and 0.131 for consensus; 0.048 and 0.152 for diverse). The thresholds were estimated on the training cohort by inner validation, do not describe the final model's sensitivity or specificity on that cohort, and differ from the fold-specific thresholds used in the cross-validation; threshold counts are therefore not given for the reference groups. Missing values were imputed from each specimen's own 1% quantile with the training-derived scale and a fixed seed; probabilities were insensitive to the imputation seed (within-specimen SD ≤ 0.020 over five seeds). The reference scores were generated by fold-specific models in the nested cross-validation, whereas the LoF/LoF scores were generated by the final model trained on all 44 development specimens; no AUC, test or performance estimate is derived from comparisons between them. Each specimen was scored alone, so that its imputation draws do not depend on other specimens or on input order; scoring the 16 specimens together changed the probabilities by at most 0.031 and altered the full-recall call of one specimen under the diverse strategy (SS_7) and none under the consensus strategy.
