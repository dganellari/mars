#!/bin/bash

set -x
set -e

# $1=MAINDIR
# $2=INSTALLDIR
# $3=SRCDIR
# $4=BUILDDIR


MAINDIR=${1:-${SCRATCH}/build/mars}
INSTALLDIR=${2:-${MAINDIR}/install}
SRCDIR=${3:-${MAINDIR}/src}
BUILDDIR=${4:-${MAINDIR}/build}

[[ -f ${INSTALLDIR}/build.log ]] && exit 0

MARCH=${MARCH:-"-march=native"}
BUILD_TYPE=${BUILD_TYPE:-"Release"}

BUILDSCRIPT_DIR=$(dirname $(realpath $0))

echo "BUILDSCRIPT_DIR"
echo ${BUILDSCRIPT_DIR}
echo ${1}
echo ${3}
pwd -P
ls -lh
# TRILINOS_DIR=${TRILINOS_DIR:-${SCRATCH}/build/trilinos/install}

echo "Mars directories comp:"
echo ${SRCDIR}
echo ${MAINDIR}

mkdir -p ${SRCDIR}
mkdir -p ${BUILDDIR}
mkdir -p ${INSTALLDIR}

pushd ${MAINDIR}
if [[ -f ${SRCDIR}/CMakeLists.txt ]]; then
  # it seems we have cloned it already
  pushd ${SRCDIR}
#  git pull
else
  git clone --recurse-submodules ${GIT_DEPTH} -b master https://bitbucket.org/zulianp/mars.git ${SRCDIR}
  pushd ${SRCDIR}
fi

CUDA_FLAGS=""
if [[ $BUILD_WITH_CUDA_SUPPORT == 1 ]]; then
#  if $CXX --version | grep -q GCC ; then
#    export CXX="${BUILDSCRIPT_DIR}/nvcc_wrapper"
#  fi
  PATCHES=""
  CUDA_FLAGS="-DMARS_USE_CUDA=ON"
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
      -DCMAKE_CXX_FLAGS="-fopenmp $MARCH" \
      -DCMAKE_INSTALL_PREFIX=${INSTALLDIR} \
      -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
      -DCMAKE_CXX_EXTENSIONS=OFF \
      -DCMAKE_CXX_STANDARD=14 \
      -DTRY_WITH_KOKKOS=ON \
      ${CUDA_FLAGS} \
      ${SRCDIR}
make -j16
make -j16 install

setfacl -R -m m::rwx "${INSTALLDIR}"
