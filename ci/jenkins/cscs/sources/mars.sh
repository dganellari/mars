#!/bin/bash

set -x
set -e

# $1=MAINDIR
# $2=INSTALLDIR
# $4=BUILDDIR


MAINDIR=${1:-${SCRATCH}/build/mars}
INSTALLDIR=${2:-${MAINDIR}/install}
BUILDDIR=${4:-${MAINDIR}/build}


[[ -f ${INSTALLDIR}/build.log ]] && exit 0

MARCH=${MARCH:-"-march=native"}
BUILD_TYPE=${BUILD_TYPE:-"Release"}

BUILDSCRIPT_DIR=$(dirname $(realpath $0))
# When used the jenkins pipeline the src dir is the current directory if not push or pop dir are done before
SRCDIR=$(pwd)
echo "SRCDIR"
echo $SRCDIR

mkdir -p ${BUILDDIR}
mkdir -p ${INSTALLDIR}

MARS_USE_KOKKOS="OFF"

if [[ $BUILD_WITH_SERIAL_SUPPORT == 1 ]]; then
  MARS_USE_KOKKOS="ON"
fi

if [[ $BUILD_WITH_OMP_SUPPORT == 1 ]]; then
  OMP_FLAGS="-fopenmp"
  MARS_USE_KOKKOS="ON"
fi

CUDA_FLAGS=""
if [[ $BUILD_WITH_CUDA_SUPPORT == 1 ]]; then
#  if $CXX --version | grep -q GCC ; then
#    export CXX="${BUILDSCRIPT_DIR}/nvcc_wrapper"
#  fi
  PATCHES=""
  CUDA_FLAGS="-DMARS_USE_CUDA=ON"
  MARS_USE_KOKKOS="ON"
fi

PATCHES="$PATCHES"
for patch in $PATCHES ; do
  # if reverse apply succeeds, the patch has been applied already (we negate the check, i.e. we apply only if reverse apply does not succeed)
  if ! patch --dry-run -f -R -p1 < ${BUILDSCRIPT_DIR}/${patch} ; then
    patch -N -p1 < ${BUILDSCRIPT_DIR}/${patch}
  fi
done


pushd ${BUILDDIR}
cmake -DCMAKE_VERBOSE_MAKEFILE=ON \
      -DCMAKE_CXX_FLAGS_RELWITHDEBINFO="-O2 -g" \
      -DCMAKE_CXX_FLAGS="$OMP_FLAGS $MARCH" \
      -DCMAKE_INSTALL_PREFIX=${INSTALLDIR} \
      -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
      -DCMAKE_CXX_EXTENSIONS=OFF \
      -DCMAKE_CXX_STANDARD=14 \
      -DTRY_WITH_KOKKOS=${MARS_USE_KOKKOS} \
      ${CUDA_FLAGS} \
      ${SRCDIR}
make -j16
make -j16 install

setfacl -R -m m::rwx "${INSTALLDIR}"
