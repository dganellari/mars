#!/bin/bash

# $1 main builddir
# $2 main installdir

set -e
set -x
set -o pipefail


BUILDBASE=${1:-${SCRATCH}/build/kokkos_serial}
INSTALLBASE=${2:-${PROJECT}/install/kokkos_serial}

SCRIPTBASE=$(dirname $(realpath $0))

export BUILD_TYPE=Debug
"${SCRIPTBASE}/build_kokkos_serial.sh" "$@"
