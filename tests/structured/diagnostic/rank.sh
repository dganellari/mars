#!/usr/bin/env bash
set -euo pipefail
[[ $# == 4 ]] || exit 2
mode=$1
library=$2
prefix=$3
binary=$4
rank=${SLURM_PROCID:-${OMPI_COMM_WORLD_RANK:-${PMI_RANK:-unknown}}}
if [[ -r /proc/self/status ]]; then
    while IFS= read -r line; do
        case $line in Cpus_allowed_list:*|Mems_allowed_list:*) printf '%s\n' "$line" ;; esac
    done < /proc/self/status > "$prefix.rank$rank.affinity"
fi
unset KOKKOS_TOOLS_LIBS KOKKOS_PROFILE_LIBRARY
export MARS_DIAG_PREFIX=$prefix MARS_DIAG_FENCE=0 MARS_DIAG_MPI_TRACE=0 MARS_DIAG_CUDA=0
case $mode in
    control) exec "$binary" ;;
    mpi) export MARS_DIAG_MPI_TRACE=1 ;;
    mpi-quiet) ;;
    host) export KOKKOS_TOOLS_LIBS=$library MARS_DIAG_CUDA=1 ;;
    profile) export MARS_DIAG_FENCE=1 KOKKOS_TOOLS_LIBS=$library ;;
    *) exit 2 ;;
esac
case $(uname -s) in
    Darwin) export DYLD_INSERT_LIBRARIES="$library${DYLD_INSERT_LIBRARIES:+:$DYLD_INSERT_LIBRARIES}" DYLD_FORCE_FLAT_NAMESPACE=1 ;;
    Linux) export LD_PRELOAD="$library${LD_PRELOAD:+:$LD_PRELOAD}" ;;
    *) exit 2 ;;
esac
exec "$binary"
