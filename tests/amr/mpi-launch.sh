#!/usr/bin/env bash
set -euo pipefail

# CTest supplies the rank count. Put multi-rank cases on two allocated nodes.
[[ $# -ge 3 && $1 == -n && -n ${SLURM_JOB_ID:-} ]] || exit 2
ranks=$2
shift 2
case $ranks in
    1) nodes=1; per_node=1 ;;
    2) nodes=2; per_node=1 ;;
    4) nodes=2; per_node=2 ;;
    *) echo "Unsupported validation rank count: $ranks" >&2; exit 2 ;;
esac
uenv=()
[[ -z ${MARS_AMR_UENV:-} ]] || uenv=(--uenv="$MARS_AMR_UENV")
exec srun "${uenv[@]}" --nodes="$nodes" --ntasks="$ranks" --ntasks-per-node="$per_node" \
    --cpus-per-task=8 --cpu-bind=cores --exclusive --export=ALL --kill-on-bad-exit=1 \
    --time=00:05:00 --uenv-passthrough=use --view=default "$@"
