#!/usr/bin/env bash
set -euo pipefail
bundle=/capstor/scratch/cscs/gandanie/bp3-paper-kit-20261107-v2
export PYTHONDONTWRITEBYTECODE=1
[[ $# -eq 0 ]] || { printf 'Usage: submit.sh\n' >&2; exit 2; }
/usr/bin/python3 "$bundle/kit/check_bundle.py" "$bundle"
common=(--parsable --uenv-passthrough=ignore
  --uenv=/ritom/scratch/cscs/gandanie/.uenv-images/images/fe67719844db4b44e7a05cd88bc0dd1cfd5e3aeea9ba1a0734d98184c1e3e2d2/store.squashfs:/user-environment
  --view=default -A csstaff --exclusive --ntasks-per-node=4 --cpus-per-task=8 --export=ALL
  --exclude=nid005323,nid005617,nid005669,nid005679,nid005742,nid005820,nid005854,nid006377,nid005927)
prep=$(sbatch "${common[@]}" --time=04:00:00 -N2 -n8 \
  --output=/capstor/scratch/cscs/gandanie/bp3_paper_prepare_%j.log "$bundle/paper/job.sh" prepare)
prep=${prep%%;*}
[[ $prep =~ ^[0-9]+$ ]] || { printf 'Invalid preparation ID: %s\n' "$prep" >&2; exit 2; }
printf 'Preparation: %s\n' "$prep"
prepared=/capstor/scratch/cscs/gandanie/bp3-paper/$prep
standard=$(sbatch "${common[@]}" --time=04:00:00 -N16 -n64 --dependency="afterok:$prep" --kill-on-invalid-dep=yes \
  --output=/capstor/scratch/cscs/gandanie/bp3_paper_standard_%j.log "$bundle/paper/job.sh" standard "$prepared")
standard=${standard%%;*}
[[ $standard =~ ^[0-9]+$ ]] || { printf 'Invalid comparison ID: %s\n' "$standard" >&2; exit 2; }
printf 'Standard comparison: %s\n' "$standard"
# A failed large single-GPU MFEM case must remain visible, without hiding the 256-rank pilot.
large=$(sbatch "${common[@]}" --time=01:30:00 -N64 -n256 --dependency="afterok:$prep,afterany:$standard" --kill-on-invalid-dep=yes \
  --output=/capstor/scratch/cscs/gandanie/bp3_paper_large_%j.log "$bundle/paper/job.sh" large "$prepared")
printf '256-GPU pilot and strong scaling: %s\n' "$large"
