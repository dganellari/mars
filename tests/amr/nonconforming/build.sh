#!/usr/bin/env bash
set -euo pipefail
kit=$1
results=$2
cache=/capstor/scratch/cscs/gandanie/one-solver/hip1/mars/build-cuda/CMakeCache.txt
[[ -f $cache ]] || { echo "Missing original toolchain cache: $cache" >&2; exit 2; }
cc=$(sed -n 's/^CMAKE_C_COMPILER:FILEPATH=//p' "$cache")
cxx=$(sed -n 's/^CMAKE_CXX_COMPILER:FILEPATH=//p' "$cache")
nvcc=$(sed -n 's/^CMAKE_CUDA_COMPILER:[^=]*=//p' "$cache")
hypre=$(sed -n 's/^HYPRE_INCLUDE_DIR:PATH=//p' "$cache")
for path in "$cc" "$cxx" "$nvcc" "$hypre/HYPRE_config.h"; do
    [[ -f $path ]] || { echo "Toolchain unavailable: $path; use the original linalg/25.10:v1 environment" >&2; exit 2; }
done
python3 -B "$kit/verify.py" "$kit"
# Check shared audit math with the actual host compiler before the GPU build.
"$cxx" -std=c++20 -I "$kit/source" "$kit/source/tests/amr/tet_audit_host.cpp" \
    -o "$results/tet-audit-preflight"
"$results/tet-audit-preflight"
cmake -S "$kit/source" -B "$results/build" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER="$cc" -DCMAKE_CXX_COMPILER="$cxx" \
    -DCMAKE_CUDA_COMPILER="$nvcc" -DCMAKE_CUDA_ARCHITECTURES=90 -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
    -DMARS_ENABLE_CUDA=ON -DMARS_ENABLE_UNSTRUCTURED=ON -DMARS_ENABLE_FEM_EXAMPLES=ON \
    -DMARS_ENABLE_KOKKOS=OFF -DMARS_ENABLE_HYPRE=ON -DHYPRE_INCLUDE_DIR="$hypre" \
    -DMARS_ENABLE_MPI=ON -DMARS_ENABLE_RDMA=ON -DMARS_ENABLE_DISTRIBUTED_BACKEND=ON \
    -DMARS_ENABLE_TESTS=OFF -DCMAKE_DISABLE_FIND_PACKAGE_PkgConfig=ON \
    -DFETCHCONTENT_FULLY_DISCONNECTED=ON -DCORNERSTONE_INCLUDE_DIRS="$kit/cornerstone/include" \
    -Dcornerstone-octree_SOURCE_DIR="$kit/cornerstone" \
    -Dcornerstone-octree_BINARY_DIR="$results/build/cornerstone-build" \
    -DCMAKE_CUDA_COMPILER_LAUNCHER="python3;$kit/source/tests/amr/compile_cuda.py;$results/compile-receipts"
grep -q "^HYPRE_LIBRARY:FILEPATH=$(dirname "$hypre")/lib" "$results/build/CMakeCache.txt"
cmake --build "$results/build" --target mars_amr_galerkin -j 2
cmake -S "$kit/source/tests/amr" -B "$results/host-checks" -DCMAKE_BUILD_TYPE=Debug
cmake --build "$results/host-checks" --target tet_audit_host -j 2
ctest --test-dir "$results/host-checks" -R '^tet_audit_host$' --output-on-failure \
    --output-junit "$results/host-results.xml"
sha256sum "$results/build/examples/distributed/unstructured/mars_amr_galerkin" \
    "$hypre/HYPRE_config.h" "$cc" "$cxx" "$nvcc" > "$results/build-hashes.txt"
cp "$results/build/CMakeCache.txt" "$results/CMakeCache.txt"
cp "$results/build/compile_commands.json" "$results/compile_commands.json"
python3 -B "$kit/verify.py" "$kit"
