#!/usr/bin/env bash
set -euo pipefail
bundle=/capstor/scratch/cscs/gandanie/bp3-coarse-policy-kit
export PYTHONDONTWRITEBYTECODE=1
[[ $# -eq 0 ]] || { printf 'Usage: submit.sh\n' >&2; exit 2; }
/usr/bin/python3 "$bundle/kit/check_bundle.py" "$bundle"
common=(--parsable --uenv-passthrough=ignore
  --uenv=/ritom/scratch/cscs/gandanie/.uenv-images/images/fe67719844db4b44e7a05cd88bc0dd1cfd5e3aeea9ba1a0734d98184c1e3e2d2/store.squashfs:/user-environment
  --view=default -A csstaff --exclusive --ntasks-per-node=4 --cpus-per-task=8 --export=ALL
  --exclude=nid005323,nid005617,nid005669,nid005679,nid005742,nid005820,nid005854,nid006377,nid005927)
preparation=$(sbatch "${common[@]}" --time=03:00:00 -N1 -n4 \
  --output=/capstor/scratch/cscs/gandanie/bp3_policy_prepare_%j.log "$bundle/policy/job.sh" prepare)
preparation=${preparation%%;*}
[[ $preparation =~ ^[0-9]+$ ]] || { printf 'Invalid job ID: %s\n' "$preparation" >&2; exit 2; }
printf 'Preparation job: %s\n' "$preparation"
trial=$(sbatch "${common[@]}" --time=01:30:00 -N16 -n64 --dependency="afterok:$preparation" --kill-on-invalid-dep=yes \
  --output=/capstor/scratch/cscs/gandanie/bp3_policy_trial_%j.log \
  "$bundle/policy/job.sh" trial "/capstor/scratch/cscs/gandanie/bp3-coarse-policy/$preparation")
printf 'Dependent policy comparison: %s\n' "$trial"
