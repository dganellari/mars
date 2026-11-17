#!/usr/bin/env bash
set -euo pipefail
kit=$(cd "$(dirname "$0")" && pwd)
export PYTHONDONTWRITEBYTECODE=1
python3 -B "$kit/verify.py" "$kit"
exec sbatch --account=csstaff --job-name=amr-nonconforming-diag \
    --nodes=2 --ntasks=8 --ntasks-per-node=4 --cpus-per-task=8 --mem=0 --exclusive \
    --time=02:00:00 --export=ALL --uenv-passthrough=use \
    --exclude=nid005323,nid005617,nid005669,nid005679,nid005742,nid005820,nid005854,nid005927,nid006377 \
    --output="${kit}-batch-%j.log" "$kit/batch.sh" "$kit"
