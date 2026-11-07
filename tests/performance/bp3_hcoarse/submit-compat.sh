#!/usr/bin/env bash
set -euo pipefail
bundle=/capstor/scratch/cscs/gandanie/bp3-hcoarse-compat-kit-20261106
export PYTHONDONTWRITEBYTECODE=1
[[ $# -eq 0 ]] || { printf 'Usage: submit-compat.sh\n' >&2; exit 2; }
/usr/bin/python3 "$bundle/kit/check_bundle.py" "$bundle"
common=(--parsable --uenv-passthrough=ignore
  --uenv=/ritom/scratch/cscs/gandanie/.uenv-images/images/fe67719844db4b44e7a05cd88bc0dd1cfd5e3aeea9ba1a0734d98184c1e3e2d2/store.squashfs:/user-environment
  --view=default -A csstaff --exclusive --ntasks-per-node=4 --cpus-per-task=8 --export=ALL
  --exclude=nid005323,nid005617,nid005669,nid005679,nid005742,nid005820,nid005854,nid006377,nid005927)
sbatch "${common[@]}" --time=02:00:00 -N2 -n8 \
  --output=/capstor/scratch/cscs/gandanie/bp3_hcoarse_compat_%j.log "$bundle/hcoarse/compat-job.sh"
