#!/bin/bash

# $1 main builddir
# $2 main installdir

set -e
set -x
set -o pipefail


BUILDBASE=${1:-${SCRATCH}/build/cpu}
INSTALLBASE=${2:-${PROJECT}/install/cpu}

SCRIPTBASE=$(dirname $(realpath $0))

export BUILD_TYPE=Debug
"${SCRIPTBASE}/build_cpu.sh" "$@"
