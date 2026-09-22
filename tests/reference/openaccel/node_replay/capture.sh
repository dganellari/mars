#!/usr/bin/env bash
set -euo pipefail
recipe=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
repo=$(cd -- "$recipe/../../../.." && pwd -P)
source=${1:-$(dirname -- "$repo")/OpenAccel}
baseline=${2:-$source/prgenv}
public_case=${3:-$source/mars-reference-inputs-20260920-v2}
[[ -f "$baseline/CMakeCache.txt" && -f "$public_case/manifest.json" ]] || {
    echo 'Usage: capture.sh OPENACCEL_SOURCE WORKING_BUILD PUBLIC_CHANNEL_BUNDLE' >&2
    exit 1
}
source=$(cd -- "$source" && pwd -P)
baseline=$(cd -- "$baseline" && pwd -P)
public_case=$(cd -- "$public_case" && pwd -P)
work=$(mktemp -d "$(dirname -- "$source")/OpenAccel-reference-nodes-$(date +%Y%m%d-%H%M%S)-XXXXXX")
trap 'status=$?; printf "%s\n" "$status" > "$work/exit.status"; printf "Saved reference: %s\n" "$work"' EXIT
python3 "$repo/scripts/prepare_openaccel_reference.py" --source "$source" \
    --output "$work/bundle" --public-case "$public_case" --include-nodes
python3 "$work/bundle/build_reference.py" --source "$source" --baseline-build "$baseline" \
    --destination "$work/source" --jobs 4 2>&1 | tee "$work/build.log"
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 \
    --cpus-per-task=1 --cpu-bind=cores --export=ALL --kill-on-bad-exit=1 \
    env OMP_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores \
    python3 "$work/bundle/run_openaccel_public.py" \
    --executable "$work/source/build/openaccel-3D.exe" --output "$work/run" --capture-interior
printf 'In the configured MARS CUDA build, run: bash %q %q\n' "$recipe/run.sh" "$work/run/exports/nodes"
