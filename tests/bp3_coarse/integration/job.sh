#!/usr/bin/env bash
set -euo pipefail
# Slurm executes a spool copy, so the source path must be explicit.
bundle=/capstor/scratch/cscs/gandanie/bp3-coarse-integration-retry-20261105
results=/capstor/scratch/cscs/gandanie/bp3-coarse-integration/${SLURM_JOB_ID:?}
export PYTHONDONTWRITEBYTECODE=1
case ${1:?prepare, resume or trial} in
  prepare) exec /usr/bin/python3 "$bundle/kit/run.py" prepare "$bundle" "$results" ;;
  resume) exec /usr/bin/python3 "$bundle/kit/run.py" resume "$bundle" "$results" --prepared "${2:?existing native build}" ;;
  trial) exec /usr/bin/python3 "$bundle/kit/run.py" trial "$bundle" "$results" --prepared "${2:?preparation results}" ;;
  *) printf 'Unknown phase\n' >&2; exit 2 ;;
esac
