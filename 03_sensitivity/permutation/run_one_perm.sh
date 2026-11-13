#!/bin/bash
# One label permutation of the primary-ensemble AUC reference (see permutation_design.md).
# Usage: bash run_one_perm.sh <permutation index 1..100>
# Inputs: runs/plan/perm_XXX/ (sample_info.csv with permuted labels + raw_data.xlsx), created by make_permutation_plan.R;
#         family_map.csv (not distributed; set CITRIN_FAMILY_FILE or place it in data/).
set -u
W="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"; REPO="$(cd "$W/../.." && pwd -P)"
RUNS="${CITRIN_PERMUTATION_RUNS:-$W/runs}"
P=$(printf "perm_%03d" "$1")
export CITRIN_DATA_DIR="$RUNS/plan/$P" CITRIN_FAMILY_FILE="${CITRIN_FAMILY_FILE:-${CITRIN_EXTERNAL_DIR:-$REPO/data}/family_map.csv}"
export CITRIN_ENSEMBLE_ENGINE="$REPO/01_nested_ensemble/citrin_44_nested_fixed_panel_v2_2.R"
export CITRIN_ENSEMBLE_OUTPUT="$RUNS/$P" CITRIN_ENSEMBLE_SKIP_INNER=true CITRIN_ENSEMBLE_FIT_FINAL=false CITRIN_ENSEMBLE_COMBAT=no_combat
export CITRIN_ENSEMBLE_WORKERS=1 CITRIN_ENSEMBLE_WORKER=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
if [[ -f "$CITRIN_ENSEMBLE_OUTPUT/outer_performance_all_policies.csv" ]]; then echo "[$P] already done"; exit 0; fi
mkdir -p "$CITRIN_ENSEMBLE_OUTPUT/logs"
cd "$W"; start=$(date +%s); echo "[$P] start $(date)"
Rscript "$W/citrin_nested_ensemble_familyLOFO_v1_1_permAUC.R" > "$CITRIN_ENSEMBLE_OUTPUT/logs/worker.log" 2>&1 && \
CITRIN_ENSEMBLE_AGGREGATE_ONLY=true Rscript "$W/citrin_nested_ensemble_familyLOFO_v1_1_permAUC.R" > "$CITRIN_ENSEMBLE_OUTPUT/logs/aggregate.log" 2>&1
echo "[$P] exit $? after $(( $(date +%s) - start )) s"
