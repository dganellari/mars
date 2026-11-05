#!/usr/bin/env bash
set -euo pipefail
bundle=/capstor/scratch/cscs/gandanie/bp3-coarse-integration-retry-20261105
export PYTHONDONTWRITEBYTECODE=1
/usr/bin/python3 "$bundle/kit/check_bundle.py" "$bundle"
common=(--parsable --uenv-passthrough=ignore
  --uenv=/ritom/scratch/cscs/gandanie/.uenv-images/images/fe67719844db4b44e7a05cd88bc0dd1cfd5e3aeea9ba1a0734d98184c1e3e2d2/store.squashfs:/user-environment
  --view=default -A csstaff --exclusive --ntasks-per-node=4 --cpus-per-task=8 --export=ALL
  --exclude=nid005323,nid005617,nid005669,nid005679,nid005742,nid005820,nid005854,nid006377,nid005927)
phase=(prepare)
limit=02:00:00
if [[ $# -eq 1 ]]; then
  phase=(resume "$1")
  limit=00:45:00
elif [[ $# -gt 1 ]]; then
  printf 'Usage: submit.sh [existing native build directory]\n' >&2; exit 2
fi
preparation=$(sbatch "${common[@]}" --time="$limit" -N1 -n4 \
  --output=/capstor/scratch/cscs/gandanie/bp3_integrate_prepare_%j.log "$bundle/kit/job.sh" "${phase[@]}")
preparation=${preparation%%;*}
[[ $preparation =~ ^[0-9]+$ ]] || { printf 'Invalid job ID: %s\n' "$preparation" >&2; exit 2; }
printf 'Preparation job: %s\n' "$preparation"
trial=$(sbatch "${common[@]}" --time=01:00:00 -N16 -n64 --dependency="afterok:$preparation" --kill-on-invalid-dep=yes \
  --output=/capstor/scratch/cscs/gandanie/bp3_integrate_trial_%j.log \
  "$bundle/kit/job.sh" trial "/capstor/scratch/cscs/gandanie/bp3-coarse-integration/$preparation")
printf 'Dependent integration comparison: %s\n' "$trial"
