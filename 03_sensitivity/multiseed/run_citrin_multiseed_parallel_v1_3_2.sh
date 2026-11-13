#!/usr/bin/env bash
# Launcher for citrin_multiseed_verification_v1_3_2.R
#
# Differences from ChatGPT's v1.2 launcher:
#   - Detaches itself (nohup) and runs under `caffeinate -is`, so the Mac does
#     not go to sleep during the ~20 h run and closing the terminal is harmless.
#   - Engine paths default to the same directory as this launcher (override with
#     VERIFY_SCRIPT44 / VERIFY_SCRIPTR3).
#   - Same env-var protocol, same single final aggregation step as v1.2.
#
# Usage (first start and every resume are the same command):
#   export CITRIN_DATA_DIR=...; export CITRIN_FAMILY_FILE=...
#   [VERIFY_* overrides] bash run_citrin_multiseed_parallel_v1_3_2.sh
# Status at any time:
#   VERIFY_OUTPUT=<dir> VERIFY_WORKER=1 VERIFY_AGGREGATE_ONLY=1 Rscript citrin_multiseed_verification_v1_3_2.R
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HARNESS="$SCRIPT_DIR/citrin_multiseed_verification_v1_3_2.R"
PKG_DEFAULT="$SCRIPT_DIR"

: "${CITRIN_DATA_DIR:?Set CITRIN_DATA_DIR to the directory containing raw_data.xlsx and sample_info.csv}"
: "${CITRIN_FAMILY_FILE:?Set CITRIN_FAMILY_FILE to the verified family_map.csv}"
export VERIFY_SCRIPT44="${VERIFY_SCRIPT44:-$PKG_DEFAULT/citrin_44_nested_fixed_panel_v2_2.R}"
export VERIFY_SCRIPTR3="${VERIFY_SCRIPTR3:-$PKG_DEFAULT/citrin_Run3_nested_fixed_panel_v2_2.R}"

# Validate paths BEFORE detaching: otherwise the caller sees a success message
# even when every worker will fail immediately in the background.
command -v Rscript >/dev/null 2>&1 || { echo "Rscript not found in PATH" >&2; exit 2; }
for f in "$HARNESS" "$VERIFY_SCRIPT44" "$VERIFY_SCRIPTR3" \
         "$CITRIN_DATA_DIR/raw_data.xlsx" "$CITRIN_DATA_DIR/sample_info.csv" \
         "$CITRIN_FAMILY_FILE"; do
  [[ -f "$f" ]] || { echo "Required file not found: $f" >&2; exit 2; }
done

export VERIFY_SEEDS="${VERIFY_SEEDS:-2025,2026,2027,2028,2029,2030,2031,2032,2033,2034}"
export VERIFY_MODES="${VERIFY_MODES:-fixed30,select}"
export VERIFY_SUBSAMPLE="${VERIFY_SUBSAMPLE:-sample}"
export VERIFY_SCOPE="${VERIFY_SCOPE:-44,Run3}"
export VERIFY_N_BOOT="${VERIFY_N_BOOT:-0}"
export VERIFY_OUTER_LIMIT="${VERIFY_OUTER_LIMIT:-0}"
export VERIFY_OUTPUT="${VERIFY_OUTPUT:-$HOME/citrin_multiseed_v1_3_2_full}"
export VERIFY_WORKERS="${VERIFY_WORKERS:-6}"

# One BLAS/OpenMP thread per R process (6 workers on a 12-core machine).
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export MKL_NUM_THREADS=1

mkdir -p "$VERIFY_OUTPUT/logs"
# Only one parent launcher may supervise a given output directory at once.
# The detached child atomically creates this lock; do not remove a live lock.
if [[ -d "$VERIFY_OUTPUT/.launcher_lock" ]]; then
  echo "Launcher lock exists: $VERIFY_OUTPUT/.launcher_lock" >&2
  echo "Check whether a previous launcher is active before restarting." >&2
  exit 2
fi

# ---- detach: re-exec this script in the background under caffeinate ----
if [[ -z "${CITRIN_LAUNCHER_DETACHED:-}" ]]; then
  export CITRIN_LAUNCHER_DETACHED=1
  if command -v caffeinate >/dev/null 2>&1; then
    nohup caffeinate -is bash "${BASH_SOURCE[0]}" >>"$VERIFY_OUTPUT/logs/launcher.log" 2>&1 &
  else
    nohup bash "${BASH_SOURCE[0]}" >>"$VERIFY_OUTPUT/logs/launcher.log" 2>&1 &
  fi
  echo "Launched in background (pid $!) with $VERIFY_WORKERS workers."
  echo "  plan:     seeds=$VERIFY_SEEDS modes=$VERIFY_MODES subsample=$VERIFY_SUBSAMPLE scope=$VERIFY_SCOPE outer_limit=$VERIFY_OUTER_LIMIT"
  echo "  output:   $VERIFY_OUTPUT"
  echo "  logs:     $VERIFY_OUTPUT/logs/launcher.log, worker_N.log, aggregate.log"
  echo "  progress: ls $VERIFY_OUTPUT/runs/*/RUN_COMPLETE 2>/dev/null | wc -l"
  echo "  status:   VERIFY_OUTPUT=\"$VERIFY_OUTPUT\" VERIFY_WORKER=1 VERIFY_AGGREGATE_ONLY=1 Rscript \"$HARNESS\""
  exit 0
fi

# ---- detached part ----
LAUNCH_LOCK="$VERIFY_OUTPUT/.launcher_lock"
if ! mkdir "$LAUNCH_LOCK" 2>/dev/null; then
  echo "Another launcher is active or left a stale launcher lock: $LAUNCH_LOCK" >&2
  exit 2
fi
printf '%s\n' "$$" > "$LAUNCH_LOCK/pid"
trap 'rm -rf "$LAUNCH_LOCK"' EXIT
printf '[%s] Plan: seeds=%s modes=%s subsampling=%s scope=%s outer_limit=%s workers=%s output=%s\n' \
  "$(date)" "$VERIFY_SEEDS" "$VERIFY_MODES" "$VERIFY_SUBSAMPLE" "$VERIFY_SCOPE" \
  "$VERIFY_OUTER_LIMIT" "$VERIFY_WORKERS" "$VERIFY_OUTPUT"

pids=()
for ((w=1; w<=VERIFY_WORKERS; w++)); do
  VERIFY_WORKER="$w" VERIFY_AGGREGATE_ONLY=0 \
    Rscript --vanilla "$HARNESS" >>"$VERIFY_OUTPUT/logs/worker_${w}.log" 2>&1 &
  pids+=("$!")
  sleep 2
done

failed=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then failed=1; fi
done
printf '[%s] All workers exited (worker failure flag=%s). Aggregating...\n' "$(date)" "$failed"

if ! VERIFY_WORKER=1 VERIFY_AGGREGATE_ONLY=1 VERIFY_REQUIRE_COMPLETE=1 \
  Rscript --vanilla "$HARNESS" >>"$VERIFY_OUTPUT/logs/aggregate.log" 2>&1; then
  failed=1
fi

printf '[%s] Done. Inspect %s/logs/aggregate.log and run_status.csv\n' "$(date)" "$VERIFY_OUTPUT"
if (( failed != 0 )); then
  printf 'At least one run failed or the summary is incomplete. Rerun the same launcher to resume.\n' >&2
fi
exit "$failed"
