#!/usr/bin/env bash
set -euo pipefail
bundle=/capstor/scratch/cscs/gandanie/bp3-paper-kit-20261107
results=/capstor/scratch/cscs/gandanie/bp3-paper/${SLURM_JOB_ID:?}
export PYTHONDONTWRITEBYTECODE=1
exec /usr/bin/python3 "$bundle/paper/driver.py" "${1:?phase}" "$bundle" "$results" "${@:2}"
