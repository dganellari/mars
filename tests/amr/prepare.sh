#!/usr/bin/env bash
set -euo pipefail

[[ $# -ge 3 && $# -le 4 ]] || { echo 'Usage: prepare.sh cpu|cuda|hip WORK CORNERSTONE_SOURCE [KOKKOS_SOURCE]' >&2; exit 2; }
backend=$1
work=$2
cornerstone=$3
source=$(cd "$(dirname "$0")/../.." && pwd)
mkdir -p "$work"
work=$(cd "$work" && pwd)
cornerstone=$(cd "$cornerstone" && pwd)
[[ -f $cornerstone/include/cstone/domain/domain.hpp ]] || exit 2
jobs=${MARS_AMR_BUILD_JOBS:-8}
cxx=${CXX:-mpicxx}
common=(-DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER="$cxx" -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
        -DCMAKE_INSTALL_PREFIX="$work/install" -DMARS_ENABLE_TESTS=OFF -DMARS_ENABLE_EXAMPLES=OFF
        -DMARS_ENABLE_FEM_EXAMPLES=OFF
        -DMARS_ENABLE_SERIAL_BACKEND=OFF -DMARS_ENABLE_HYPRE=OFF)
checks=(-DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER="$cxx" -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
        -DCMAKE_PREFIX_PATH="$work/install"
        -DMPIEXEC_EXECUTABLE="${MARS_AMR_MPI_LAUNCHER:-mpiexec}" -DMPIEXEC_NUMPROC_FLAG=-n)
if [[ ${MARS_AMR_EXODUS_CHECKS:-OFF} == ON ]]; then
    checks+=(-DMARS_AMR_EXODUS_CHECKS=ON -DMARS_CSTONE_INCLUDE_DIR="$cornerstone/include")
fi
case $backend in
    cpu)
        [[ $# -eq 4 ]] || { echo 'CPU validation also needs the pinned Kokkos source.' >&2; exit 2; }
        kokkos=$(cd "$4" && pwd)
        cmake -S "$kokkos" -B "$work/kokkos-build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=20 \
            -DCMAKE_CXX_COMPILER="$cxx" -DCMAKE_INSTALL_PREFIX="$work/kokkos-install" -DCMAKE_INSTALL_LIBDIR=lib \
            -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=OFF -DKokkos_ENABLE_CUDA=OFF \
            -DKokkos_ENABLE_HIP=OFF -DKokkos_ENABLE_TESTS=OFF
        cmake --build "$work/kokkos-build" -j "$jobs"
        cmake --install "$work/kokkos-build"
        common+=(-DMARS_ENABLE_KOKKOS=ON -DMARS_ENABLE_AMR_BACKEND=ON -DMARS_ENABLE_UNSTRUCTURED=OFF
                 -DMARS_ENABLE_DISTRIBUTED_BACKEND=OFF
                 -DMARS_ENABLE_MPI=OFF -DCMAKE_PREFIX_PATH="$work/kokkos-install"
                 -DKokkos_DIR="$work/kokkos-install/lib/cmake/Kokkos")
        checks+=(-DMARS_AMR_LEGACY_CHECKS=ON -DMARS_AMR_MPI_CHECKS=ON
                 -DKokkos_DIR="$work/kokkos-install/lib/cmake/Kokkos"
                 -DMARS_CSTONE_INCLUDE_DIR="$cornerstone/include")
        ;;
    cuda|hip)
        common+=(-DMARS_ENABLE_KOKKOS=OFF -DMARS_ENABLE_AMR_BACKEND=OFF -DMARS_ENABLE_MPI=ON
                 -DMARS_ENABLE_DISTRIBUTED_BACKEND=ON
                 -DMARS_ENABLE_UNSTRUCTURED=ON -DFETCHCONTENT_FULLY_DISCONNECTED=ON
                 -DCORNERSTONE_INCLUDE_DIRS="$cornerstone/include"
                 -Dcornerstone-octree_SOURCE_DIR="$cornerstone"
                 -Dcornerstone-octree_BINARY_DIR="$work/mars-build/cornerstone-build")
        if [[ $backend == cuda ]]; then
            common+=(-DMARS_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES="${MARS_AMR_GPU_ARCH:-90}")
        else
            common+=(-DMARS_ENABLE_HIP=ON -DHIP_GPU_ARCHITECTURES="${MARS_AMR_GPU_ARCH:-gfx942}")
        fi
        checks+=(-DMARS_AMR_GPU_CHECKS=ON)
        ;;
    *) echo "Unknown backend: $backend" >&2; exit 2 ;;
esac
cmake -S "$source" -B "$work/mars-build" "${common[@]}"
cmake --build "$work/mars-build" -j "$jobs"
cmake --install "$work/mars-build"
cmake -G "Unix Makefiles" -S "$source/tests/amr" -B "$work/checks" "${checks[@]}"
# Keep compiler peaks separate and identify the exact CUDA unit if one fails.
cmake --build "$work/checks" -j 1 -- -k
printf '%s\n' "$backend" > "$work/backend.txt"
echo "Prepared $backend validation in $work/checks; no runtime pass is implied."
