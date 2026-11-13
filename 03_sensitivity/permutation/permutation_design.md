# Label-permutation null for the primary consensus30 ensemble AUC (fixed design, 2026-11-07)

Decided before any permutation result was seen (user + ChatGPT review, 2026-11-07 15:30 JST).

- Test statistic: out-of-fold AUC of the consensus30 ensemble (primary). diverse30 AUC recorded as auxiliary only; the primary is fixed in advance.
- N = 100 permutations, perm_001..perm_100 of the existing constrained plan (`runs/plan/`, written by `make_permutation_plan.R`): the seven positive labels are assigned to five randomly chosen Run 3 families of sizes 1,1,1,2,2; measurements, family map and run composition unchanged. N will not be changed after seeing results. P = (b + 1) / 101.
- Pipeline: `citrin_nested_ensemble_familyLOFO_v1_1_permAUC.R` = production harness v1_1 with one switch (CITRIN_ENSEMBLE_SKIP_INNER=true) that skips the inner family-LOFO calibration, which is used only for decision thresholds and the adaptive choice. Outer training (preprocessing, stability selection, consensus-30 panel, 10 elastic-net members with seeds 2025-2034, lambda selection) is re-run from scratch inside every outer fold of every permutation; no panel, coefficient or threshold is carried over. Engine v2_2 unchanged (md5 fbff0fa82f0b2da1dafb553bc5564d27). No final model fit.
- Equivalence check (`equivalence_check/`): original labels, SKIP_INNER=true; the out-of-fold probabilities of consensus30 and diverse30 for all 44 specimens must equal the production run (the corresponding run folder (not distributed)) to numerical precision. Result recorded below.
- Each permutation has its own output directory (`runs/perm_k`) and family cache; nothing is shared between permutations.
- Launch: `launch_perms.sh` (NW concurrent single-worker runs) waits for the WWfix ablation queue, then runs perm_001..100; `collect_perm_auc.R` writes `permutation_null_AUC.csv` and `permutation_summary.csv`.
- The previously planned re-run of the single-model (select, seed 2025) permutation was cancelled (the corresponding run folder (not distributed) disabled, original kept with suffix .DISABLED_20261107_...). Its old-data result (AUC 0.637, P 0.26) will be removed from the Response.

## Equivalence check result (2026-11-07 15:44)
- Original labels, SKIP_INNER=true, 4 workers: all 44 out-of-fold probabilities of consensus30 and of diverse30 are bit-identical to the production run (max |difference| = 0.0; 44/44 exact for each); the consensus-30 panels are identical in all 35 outer folds (1,050 fold-protein pairs). AUC consensus30 0.833977, diverse30 0.880309, Brier 0.1119 / 0.1072, family AUC 0.773 / 0.847 — identical to production.
- Cost: 1.05 core-hours per full outer LOFO under current machine load (mean 1.8 min per fold), i.e. about 1/27 of the full harness. 100 permutations ≈ 105 core-hours.
- Launched `launch_perms.sh` with NW=8 (memory-safe alongside the 4 multiseed workers); it starts when the WWfix ablation queue has written "ablation compare exit".
