#!/bin/bash
# Exploratory: nested family-held-out ten-member ensemble re-trained on the 17 carrier / LoF-Mis specimens only
# (same harness as the primary analysis with the cohort-size checks taken from the environment; engine unchanged).
# Inputs: a data folder holding raw_data.xlsx, a sample_info.csv restricted to the 17 WS_/SM_ specimens and family_map.csv
#         (prepared by abundance_gate_two_step.R: results_rerun/carrier_vs_LoFMis_retraining_input/). family_map.csv is not distributed.
set -u
W="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"; REPO="$(cd "$W/../.." && pwd -P)"
IN="${CITRIN_SUBSET_DATA_DIR:-$W/results_rerun/carrier_vs_LoFMis_retraining_input}"
export CITRIN_DATA_DIR="$IN" CITRIN_FAMILY_FILE="${CITRIN_FAMILY_FILE:-$IN/family_map.csv}"
export CITRIN_ENSEMBLE_ENGINE="$REPO/01_nested_ensemble/citrin_44_nested_fixed_panel_v2_2.R"
export CITRIN_EXPECT_44=false CITRIN_ENSEMBLE_EXPECT_N=17 CITRIN_ENSEMBLE_EXPECT_FAMS=11 CITRIN_ENSEMBLE_MIN_INNER_N=10
export CITRIN_ENSEMBLE_OUTPUT="${CITRIN_ENSEMBLE_OUTPUT:-$W/results_rerun/carrier_vs_LoFMis_retraining}" CITRIN_ENSEMBLE_OUTER_LIMIT=0 CITRIN_ENSEMBLE_WORKERS=2
export CITRIN_ENSEMBLE_FIT_FINAL=false CITRIN_ENSEMBLE_COMBAT=no_combat
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
mkdir -p "$CITRIN_ENSEMBLE_OUTPUT/logs"; cd "$W"
for w in 1 2; do CITRIN_ENSEMBLE_WORKER=$w nice -n 10 Rscript "$W/citrin_nested_ensemble_familyLOFO_v1_1_subset.R" > "$CITRIN_ENSEMBLE_OUTPUT/logs/worker_$w.log" 2>&1 & done
wait
CITRIN_ENSEMBLE_WORKER=1 CITRIN_ENSEMBLE_AGGREGATE_ONLY=true Rscript "$W/citrin_nested_ensemble_familyLOFO_v1_1_subset.R" > "$CITRIN_ENSEMBLE_OUTPUT/logs/aggregate.log" 2>&1
echo "DONE $(date)" >> "$CITRIN_ENSEMBLE_OUTPUT/logs/aggregate.log"
