#!/usr/bin/env bash
# macOS: six disjoint outer-family workers, caffeinate, safe resume, aggregate once.
# Run from any directory: bash /path/to/Rproj/run_citrin_nested_ensemble_6workers_caffeinate.sh
# ENV overrides: CITRIN_DATA_DIR, CITRIN_ENSEMBLE_ENGINE, CITRIN_ENSEMBLE_OUTPUT,
#                CITRIN_ENSEMBLE_COMBAT, CITRIN_ENSEMBLE_STABILITY_B, etc.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
HARNESS="$SCRIPT_DIR/citrin_nested_ensemble_familyLOFO_v1_1.R"

# Mac: keep system awake while the launcher and all children execute.
# -i prevents idle sleep; -s adds system-sleep assertion while on AC power.
# Does NOT guarantee protection from closed-lid sleep, shutdown or power loss.
if [[ "$(uname -s)" == "Darwin" && "${CITRIN_ENSEMBLE_CAFFEINATED:-0}" != "1" ]]; then
  if ! command -v caffeinate >/dev/null 2>&1; then
    echo "ERROR: caffeinate is not on PATH" >&2
    exit 1
  fi
  export CITRIN_ENSEMBLE_CAFFEINATED=1
  exec caffeinate -is /bin/bash "$0" "$@"
fi

cd "$SCRIPT_DIR"
if [[ ! -f "$HARNESS" ]]; then echo "Missing: $HARNESS" >&2; exit 1; fi
if ! command -v Rscript >/dev/null 2>&1; then echo "Rscript not found in PATH" >&2; exit 1; fi

export CITRIN_DATA_DIR="${CITRIN_DATA_DIR:-$SCRIPT_DIR/../data}"
export CITRIN_FAMILY_FILE="${CITRIN_FAMILY_FILE:-$CITRIN_DATA_DIR/family_map.csv}"
if [[ -z "${CITRIN_ENSEMBLE_ENGINE:-}" ]]; then
  if [[ -f "$SCRIPT_DIR/citrin_44_nested_fixed_panel_v2_2.R" ]]; then
    export CITRIN_ENSEMBLE_ENGINE="$SCRIPT_DIR/citrin_44_nested_fixed_panel_v2_2.R"
  elif [[ -f "$SCRIPT_DIR/Claude/citrin_44_nested_fixed_panel_v2_2.R" ]]; then
    export CITRIN_ENSEMBLE_ENGINE="$SCRIPT_DIR/Claude/citrin_44_nested_fixed_panel_v2_2.R"
  else
    echo "v2.2 engine not found in Rproj directory or Claude/" >&2
    exit 1
  fi
fi
export CITRIN_ENSEMBLE_COMBAT="${CITRIN_ENSEMBLE_COMBAT:-no_combat}"
export CITRIN_ENSEMBLE_SEEDS="${CITRIN_ENSEMBLE_SEEDS:-2025,2026,2027,2028,2029,2030,2031,2032,2033,2034}"
export CITRIN_ENSEMBLE_STABILITY_B="${CITRIN_ENSEMBLE_STABILITY_B:-30}"
export CITRIN_ENSEMBLE_SUBSAMPLE="${CITRIN_ENSEMBLE_SUBSAMPLE:-sample}"
export CITRIN_ENSEMBLE_OUTER_LIMIT="${CITRIN_ENSEMBLE_OUTER_LIMIT:-0}"
export CITRIN_ENSEMBLE_WORKERS=6
export CITRIN_ENSEMBLE_FIT_FINAL="${CITRIN_ENSEMBLE_FIT_FINAL:-true}"
export CITRIN_ENSEMBLE_OUTPUT="${CITRIN_ENSEMBLE_OUTPUT:-$SCRIPT_DIR/results_rerun/citrin_nested_ensemble_v1_1_${CITRIN_ENSEMBLE_COMBAT}_innerLOFO}"
# Old v1's four-fold setting has no meaning in v1.1, even if exported before.
unset CITRIN_ENSEMBLE_INNER_FOLDS
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

if [[ ! -f "$CITRIN_DATA_DIR/raw_data.xlsx" || ! -f "$CITRIN_DATA_DIR/sample_info.csv" ||
      ! -f "$CITRIN_FAMILY_FILE" || ! -f "$CITRIN_ENSEMBLE_ENGINE" ]]; then
  echo "Missing input(s). Verify DATA_DIR, FAMILY_FILE, and ENGINE path." >&2
  echo "  DATA_DIR=$CITRIN_DATA_DIR" >&2
  echo "  FAMILY_FILE=$CITRIN_FAMILY_FILE" >&2
  echo "  ENGINE=$CITRIN_ENSEMBLE_ENGINE" >&2
  exit 1
fi

mkdir -p "$CITRIN_ENSEMBLE_OUTPUT/logs"
LOCK_DIR="$CITRIN_ENSEMBLE_OUTPUT/.launcher_lock"
THIS_HOST="$(hostname)"

acquire_launcher_lock() {
  local owner_host owner_pid old_dir
  if mkdir "$LOCK_DIR" 2>/dev/null; then
    printf 'host=%s\npid=%s\nutc=%s\n' "$THIS_HOST" "$$" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" \
      > "$LOCK_DIR/owner.txt"
    return 0
  fi
  if [[ ! -f "$LOCK_DIR/owner.txt" ]]; then
    echo "Launcher lock present but owner missing. Inspect manually: $LOCK_DIR" >&2
    return 1
  fi
  owner_host="$(sed -n '1s/^host=//p' "$LOCK_DIR/owner.txt")"
  owner_pid="$(sed -n '2s/^pid=//p' "$LOCK_DIR/owner.txt")"
  if [[ "$owner_host" != "$THIS_HOST" || ! "$owner_pid" =~ ^[0-9]+$ ]]; then
    echo "Launcher lock belongs to another/unknown host or PID; do NOT auto-delete: $LOCK_DIR" >&2
    return 1
  fi
  if kill -0 "$owner_pid" 2>/dev/null; then
    echo "Another launcher is active (PID $owner_pid); refusing duplicate run." >&2
    return 1
  fi
  old_dir="${LOCK_DIR}.abandoned.$$.${RANDOM}"
  if ! mv "$LOCK_DIR" "$old_dir"; then
    echo "Unable to rename dead launcher's lock. Inspect: $LOCK_DIR" >&2
    return 1
  fi
  # Do not delete if the owner changed while acquiring the abandoned lock.
  if [[ "$(sed -n '2s/^pid=//p' "$old_dir/owner.txt")" != "$owner_pid" ]]; then
    echo "Owner changed during lock recovery: $old_dir" >&2
    return 1
  fi
  rm -rf "$old_dir"
  echo "[RECOVER] Dead launcher lock (PID $owner_pid) removed." >&2
  mkdir "$LOCK_DIR"
  printf 'host=%s\npid=%s\nutc=%s\n' "$THIS_HOST" "$$" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" \
    > "$LOCK_DIR/owner.txt"
}

acquire_launcher_lock
children=()
cleanup() {
  local rc=$?
  trap - EXIT
  # bash 3.2 (macOS default) treats an empty array as unbound under set -u.
  for pid in ${children[@]+"${children[@]}"}; do
    if kill -0 "$pid" 2>/dev/null; then kill "$pid" 2>/dev/null || true; fi
  done
  for pid in ${children[@]+"${children[@]}"}; do wait "$pid" 2>/dev/null || true; done
  # Only the creating launcher may remove its lock.
  if [[ -f "$LOCK_DIR/owner.txt" ]] && \
      [[ "$(sed -n '2s/^pid=//p' "$LOCK_DIR/owner.txt")" == "$$" ]]; then
    rm -rf "$LOCK_DIR"
  fi
  exit "$rc"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
trap 'exit 129' HUP

LOGDIR="$CITRIN_ENSEMBLE_OUTPUT/logs/$(date '+%Y%m%d_%H%M%S')_pid$$"
mkdir -p "$LOGDIR"
echo "[START] $(date) v1.1 inner family-LOFO; workers=6; B=$CITRIN_ENSEMBLE_STABILITY_B"
echo "[OUTPUT] $CITRIN_ENSEMBLE_OUTPUT"
echo "[LOGS] $LOGDIR"
echo "[MODE] $CITRIN_ENSEMBLE_COMBAT"
echo "[NOTICE] Re-run with EXACT same settings and output folder to resume."

for j in 1 2 3 4 5 6; do
  CITRIN_ENSEMBLE_WORKER="$j" CITRIN_ENSEMBLE_AGGREGATE_ONLY=false \
    Rscript "$HARNESS" > "$LOGDIR/worker_${j}.log" 2>&1 &
  children+=("$!")
  sleep 0.3
done

failed=0
for idx in 0 1 2 3 4 5; do
  pid="${children[$idx]}"
  if wait "$pid"; then
    echo "[OK] Worker $((idx+1)) completed"
  else
    echo "[FAILED] Worker $((idx+1)); see $LOGDIR/worker_$((idx+1)).log" >&2
    failed=1
  fi
done
children=()
if [[ "$failed" -ne 0 ]]; then
  echo "[INCOMPLETE] Worker failure: no final aggregation. Correct the issue and resume." >&2
  exit 1
fi

CITRIN_ENSEMBLE_WORKER=1 CITRIN_ENSEMBLE_AGGREGATE_ONLY=true \
  Rscript "$HARNESS" > "$LOGDIR/aggregate.log" 2>&1 &
children=("$!")
if wait "${children[0]}"; then
  children=()
  echo "[SUCCESS] Complete. Inspect: $CITRIN_ENSEMBLE_OUTPUT/outer_performance_all_policies.csv"
  echo "[AGGREGATION LOG] $LOGDIR/aggregate.log"
else
  children=()
  echo "[AGGREGATION FAILED] Inspect $LOGDIR/aggregate.log" >&2
  exit 1
fi
