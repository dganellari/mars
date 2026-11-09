#!/usr/bin/env bash
set -euo pipefail
[[ $# -eq 2 && -n ${SLURM_JOB_ID:-} ]] || exit 2
backend=$1
# Slurm copies this script into its spool directory; the submitter supplies the bundle.
bundle=$2
[[ $bundle == /* && -f $bundle/manifest.json ]] || exit 2
export PYTHONDONTWRITEBYTECODE=1
python3 "$bundle/verify.py"
work="${bundle%/*}/${bundle##*/}-results-$backend-$SLURM_JOB_ID"
args=("$backend" "$work" "$bundle/cornerstone")
[[ $backend != cpu ]] || args+=("$bundle/kokkos")
status=0
bash "$bundle/source/tests/amr/cluster.sh" "${args[@]}" || status=$?
python3 "$bundle/verify.py" || status=$?
if [[ $status -eq 0 ]]; then python3 "$bundle/audit.py" "$backend" "$work" || status=$?; fi
mkdir -p "$work"
printf '%s\n' "$status" > "$work/batch-exit.txt"
cp "$bundle/manifest.json" "$work/source-manifest.json"
if [[ $status -ne 0 && -f $work/build.log ]]; then tail -80 "$work/build.log"; fi
exit "$status"
