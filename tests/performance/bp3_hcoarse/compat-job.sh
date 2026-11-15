#!/usr/bin/env bash
set -euo pipefail
bundle=/capstor/scratch/cscs/gandanie/bp3-hcoarse-compat-kit-20261106
results=/capstor/scratch/cscs/gandanie/bp3-hcoarse-compat/${SLURM_JOB_ID:?}
export PYTHONDONTWRITEBYTECODE=1
exec /usr/bin/python3 "$bundle/hcoarse/compat.py" check "$bundle" "$results"
