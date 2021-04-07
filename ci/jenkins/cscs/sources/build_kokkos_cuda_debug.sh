#!/bin/bash

# $1 main builddir
# $2 main installdir

set -e
set -x
set -o pipefail

BUILDBASE=${1:-${SCRATCH}/build/kokkos_cuda}
INSTALLBASE=${2:-${PROJECT}/install/kokkos_cuda}

SCRIPTBASE=$(dirname $(realpath $0))
export BUILD_TYPE=Debug
"${SCRIPTBASE}/build_kokkos_cuda.sh" "$@"
