#!/bin/bash
# MFEM v4.8 (MPI + CUDA, no METIS) and mfem_bp3, with the compilers and the Hypre of an
# existing MARS build, so the baseline and MARS share the toolchain and BoomerAMG.
# Run from the MARS build directory, in the environment of its build:
#   bash ../tests/performance/mfem_bp3/build_mfem.sh [prefix] [cuda arch, default 90]
# prefix (default ./mfem-v4.8) holds MFEM's source, build and install; mfem_bp3 goes to the
# current directory.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
prefix=${1:-$PWD/mfem-v4.8}
arch=${2:-90}
cache=$PWD/CMakeCache.txt
[ -f "$cache" ] || { echo "no CMakeCache.txt here: run from the MARS build directory"; exit 1; }
value() { sed -n "s/^$1:[A-Z]*=//p" "$cache" | head -1; }

cxx=$(value CMAKE_CXX_COMPILER)
# Hypre as MARS found it: pkg-config, else find_path, else $HYPRE_DIR.
hypre=$(value HYPRE_PC_PREFIX)
if [ -z "$hypre" ]; then
    inc=$(value HYPRE_INCLUDE_DIR)
    hypre=${inc:+$(dirname "$inc")}
fi
hypre=${hypre:-${HYPRE_DIR:-}}
config=$hypre/include/HYPRE_config.h
[ -f "$config" ] || { echo "Hypre not found (prefix '$hypre')"; exit 1; }
echo "compiler $cxx"
echo "Hypre $hypre"
grep -E "define HYPRE_(RELEASE_VERSION|USING_CUDA|USING_UMPIRE|USING_GPU_AWARE_MPI)" "$config" || true
grep -q "define HYPRE_USING_CUDA" "$config" || { echo "this Hypre is not a CUDA build"; exit 1; }

mkdir -p "$prefix"
[ -d "$prefix/src" ] || git clone --depth 1 -b v4.8 https://github.com/mfem/mfem.git "$prefix/src"
cmake -S "$prefix/src" -B "$prefix/build" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$prefix/install" \
    -DCMAKE_CXX_COMPILER="$cxx" \
    -DMFEM_USE_MPI=YES -DMFEM_USE_METIS=NO -DMFEM_USE_CUDA=YES -DCUDA_ARCH=sm_$arch \
    -DHYPRE_DIR="$hypre" -DMFEM_ENABLE_EXAMPLES=NO -DMFEM_ENABLE_MINIAPPS=NO
cmake --build "$prefix/build" -j 32
cmake --install "$prefix/build"
make -f "$here/Makefile" MFEM_DIR="$prefix/install"
echo "MFEM BUILD OK"
