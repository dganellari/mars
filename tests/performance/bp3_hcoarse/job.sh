#!/usr/bin/env bash
set -euo pipefail
bundle=/capstor/scratch/cscs/gandanie/bp3-hcoarse-resume-kit-20261106
results=/capstor/scratch/cscs/gandanie/bp3-hcoarse-integrated/${SLURM_JOB_ID:?}
export PYTHONDONTWRITEBYTECODE=1
case ${1:?prepare, resume or trial} in
  prepare) exec /usr/bin/python3 "$bundle/hcoarse/integrate.py" prepare "$bundle" "$results" ;;
  resume) exec /usr/bin/python3 "$bundle/hcoarse/integrate.py" resume "$bundle" "$results" --prepared "${2:?previous results}" ;;
  trial) exec /usr/bin/python3 "$bundle/hcoarse/integrate.py" trial "$bundle" "$results" --prepared "${2:?preparation results}" ;;
  *) printf 'Unknown phase\n' >&2; exit 2 ;;
esac
