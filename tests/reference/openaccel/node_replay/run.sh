#!/usr/bin/env bash
set -euo pipefail
recipe=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
repo=$(cd -- "$recipe/../../../.." && pwd -P)
[[ $# -ge 1 && $# -le 2 ]] || { echo 'Usage: run.sh PUBLIC_EXPORTS_NODES [MARS_BUILD]' >&2; exit 1; }
capture=$1
build_root=${2:-$PWD}
build_root=$(cd -- "$build_root" && pwd -P)
replay="$build_root/examples/distributed/unstructured/mars_segregated_node_replay"
algebra="$build_root/examples/distributed/unstructured/mars_segregated_node_algebra_check"
[[ -x "$replay" && -x "$algebra" ]] || {
    echo "Build first: cmake --build '$build_root' --target mars_segregated_node_replay mars_segregated_node_algebra_check -j4" >&2
    exit 1
}
run_dir=$(mktemp -d "$build_root/node-replay-$(date +%Y%m%d-%H%M%S)-XXXXXX")
trap 'status=$?; printf "%s\n" "$status" > "$run_dir/exit.status"; printf "Saved replay: %s\n" "$run_dir"' EXIT
python3 "$repo/scripts/openaccel_node_check.py" "$capture" --pack "$run_dir/inputs.txt" 2>&1 | tee "$run_dir/prepare.log"
"$algebra" | tee "$run_dir/algebra.log"
git -C "$repo" rev-parse HEAD > "$run_dir/revision.txt"
sha256sum "$replay" "$algebra" "$run_dir/inputs.txt" \
    "$repo/backend/distributed/unstructured/fem/segregated/mars_segregated_node.hpp" \
    "$repo/examples/distributed/unstructured/mars_segregated_node_replay.cu" > "$run_dir/run-sha256.txt"
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 \
    --export=ALL --kill-on-bad-exit=1 \
    "$HOME/affinity/bind_numa.sh" "$replay" "$run_dir/inputs.txt" 2>&1 | tee "$run_dir/run.log"
