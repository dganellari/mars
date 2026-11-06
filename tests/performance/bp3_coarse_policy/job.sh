#!/usr/bin/env bash
set -euo pipefail
bundle=/capstor/scratch/cscs/gandanie/bp3-coarse-policy-kit
results=/capstor/scratch/cscs/gandanie/bp3-coarse-policy/${SLURM_JOB_ID:?}
export PYTHONDONTWRITEBYTECODE=1
case ${1:?prepare or trial} in
  prepare) exec /usr/bin/python3 "$bundle/policy/policy.py" prepare "$bundle" "$results" ;;
  trial) exec /usr/bin/python3 "$bundle/policy/policy.py" trial "$bundle" "$results" --prepared "${2:?preparation results}" ;;
  *) printf 'Unknown phase\n' >&2; exit 2 ;;
esac
