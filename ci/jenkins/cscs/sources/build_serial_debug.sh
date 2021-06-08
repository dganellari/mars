#!/bin/bash

# $1 main builddir
# $2 main installdir

set -e
set -x
set -o pipefail


BUILDBASE=${1:-${SCRATCH}/build/serial}
INSTALLBASE=${2:-${PROJECT}/install/serial}

SCRIPTBASE=$(dirname $(realpath $0))

export BUILD_TYPE=Debug
"${SCRIPTBASE}/build_serial.sh" "$@"
