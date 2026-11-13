# Table S15. Protein-exclusion sensitivity analysis of the exploratory consensus-30 ensemble.

| Excluded protein(s) | AUC (95% CI) | Sensitivity (95% CI) | Specificity (95% CI) | Balanced accuracy (95% CI) | LoF/Mis families detected | Brier score | Distinct proteins in the 35 consensus panels |
|---|---|---|---|---|---|---|---|
| None (primary analysis) | 0.834 (0.602–0.971) | 7/7 (0.59–1.00) | 31/37 (0.68–0.94) | 0.919 (0.848–0.974) | 5/5 | 0.112 |  88 |
| PRL | 0.888 (0.752–0.982) | 6/7 (0.42–1.00) | 29/37 (0.62–0.90) | 0.820 (0.550–0.943) | 4/5 | 0.106 |  92 |
| MTTP | 0.903 (0.780–0.987) | 6/7 (0.42–1.00) | 29/37 (0.62–0.90) | 0.820 (0.541–0.943) | 4/5 | 0.101 |  99 |
| ARHGEF10 | 0.873 (0.748–0.966) | 7/7 (0.59–1.00) | 23/37 (0.45–0.78) | 0.811 (0.735–0.878) | 5/5 | 0.115 |  99 |
| SEC31A | 0.838 (0.652–0.952) | 6/7 (0.42–1.00) | 13/37 (0.20–0.53) | 0.604 (0.323–0.742) | 4/5 | 0.128 |  97 |
| PRL, MTTP, ARHGEF10 and SEC31A | 0.764 (0.500–0.953) | 7/7 (0.59–1.00) | 14/37 (0.22–0.55) | 0.689 (0.611–0.766) | 5/5 | 0.127 | 104 |

**Misclassified specimens by scenario**

| Scenario | False positives (controls called positive) | Missed LoF/Mis |
|---|---|---|
| None (primary analysis) | DD_1, DD_10, DD_19, DD_6, WS_9, WW_6 | none |
| PRL | DD_1, DD_10, DD_19, DD_2, DD_6, WW_3, WW_4, WW_6 | SM_6 |
| MTTP | DD_1, DD_10, DD_19, DD_6, DD_8, WS_4, WW_4, WW_6 | SM_6 |
| ARHGEF10 | DD_1, DD_10, DD_18, DD_19, DD_21, DD_6, DD_7, WS_3, WS_6, WS_9, WW_2, WW_3, WW_4, WW_6 | none |
| SEC31A | DD_1, DD_10, DD_12, DD_13, DD_16, DD_17, DD_18, DD_19, DD_2, DD_21, DD_4, DD_5, DD_6, DD_7, DD_8, DD_9, WS_10, WS_3, WS_6, WS_9, WW_2, WW_3, WW_4, WW_6 | SM_6 |
| PRL, MTTP, ARHGEF10 and SEC31A | DD_1, DD_10, DD_12, DD_13, DD_16, DD_18, DD_19, DD_20, DD_21, DD_5, DD_6, DD_7, WS_10, WS_3, WS_6, WS_7, WS_8, WS_9, WW_2, WW_3, WW_4, WW_5, WW_6 | none |

**Consensus-panel membership of selected proteins (outer folds of 35)**

| Protein | baseline | no_PRL | no_MTTP | no_ARHGEF10 | no_SEC31A | all4 |
|---|---|---|---|---|---|---|
| SEC31A | 35 | 35 | 35 | 35 |  0 |  0 |
| ARHGEF10 | 35 | 35 | 35 |  0 | 35 |  0 |
| NAPB | 35 | 35 | 35 | 35 | 35 | 35 |
| TESC | 34 | 34 | 34 | 32 | 34 | 34 |
| RAB2B | 35 | 35 | 35 | 35 | 35 | 35 |
| BLTP3B | 34 | 35 | 34 | 34 | 34 | 35 |
| DOCK11 | 35 | 35 | 35 | 33 | 35 | 35 |
| PRL | 33 |  0 | 34 | 31 | 33 |  0 |
| MTTP | 35 | 35 |  0 | 35 | 35 |  0 |
| CPNE1 | 34 | 34 | 34 | 34 | 34 | 34 |
| PCYOX1 | 33 | 34 | 33 | 33 | 33 | 34 |
| PPAT | 33 | 34 | 33 | 34 | 33 | 34 |

For each scenario the complete nested model-development procedure (training-only imputation, stability selection, consensus-panel formation, member fitting and inner-loop threshold derivation; family-grouped leave-one-family-out, 35 outer folds, ten members) was re-run after removing the indicated protein(s) from the 5,199-protein candidate matrix before any selection step. Values are family-held-out out-of-fold results of the primary strategy (consensus-30, complete-inner-sensitivity threshold) on the same 44 specimens (7 LoF/Mis from 5 families; 37 controls). AUC and balanced accuracy: family-cluster bootstrap 95% CI (2,000 replicates, conditional on the fitted models; the baseline intervals are those reported in Table 2, obtained with the same method); sensitivity and specificity: exact binomial 95% CI (specimen-level, not accounting for within-family correlation). The last column gives the number of distinct proteins that entered the consensus candidate panel in at least one of the 35 outer folds.
The scenarios (joint exclusion of the four proteins highlighted by the reviewer, and each protein individually) were fixed before the runs and are all reported. These are post hoc sensitivity analyses within the development cohort: they indicate whether a protein is individually indispensable for obtaining discrimination in this cohort, but do not establish its contribution in the primary model, do not show that the remaining signal is free of confounding, and are not independent validation. AUC, Brier score and threshold-based sensitivity and specificity are distinct aspects of performance: because the inner-loop threshold is re-derived in every scenario, threshold-based counts can change substantially with little change in AUC. Selection counts of core proteins that remain at 35/35 do not imply that the rest of the panel or the fitted coefficients were unchanged, since a removed protein's place in each fold's panel is taken by another protein and the stochastic imputation differs once a feature is removed. Sheet 'Per-specimen calls' gives each specimen's out-of-fold call and probability under every scenario; sheet 'Core proteins' gives, for selected proteins, the number of outer folds (of 35) in which they entered the consensus panel under each scenario.
