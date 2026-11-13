#!/bin/bash
# Launch the 100 label permutations with NW parallel workers, then collect the null distribution.
# Usage: NW=8 bash launch_perms.sh      (FIRST / LAST select a sub-range)
set -u
W="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"; REPO="$(cd "$W/../.." && pwd -P)"
NW=${NW:-8}; FIRST=${FIRST:-1}; LAST=${LAST:-100}
if [[ "$(uname -s)" == "Darwin" && "${PERM_CAFFEINATED:-0}" != "1" ]]; then export PERM_CAFFEINATED=1; exec caffeinate -is /bin/bash "$0" "$@"; fi
seq "$FIRST" "$LAST" | xargs -P "$NW" -I{} bash "$W/run_one_perm.sh" {}
(cd "$REPO" && Rscript 03_sensitivity/permutation/collect_perm_auc.R)
