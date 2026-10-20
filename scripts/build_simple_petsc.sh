#!/usr/bin/env bash
# Build inside the MARS uenv so PETSc uses the same MPI and CUDA runtime.
set -euo pipefail
if (( $# != 1 )); then
    echo "Usage: bash scripts/build_simple_petsc.sh NEW_BUILD_DIRECTORY" >&2
    exit 2
fi
for tool in mpicc mpicxx nvcc python3 curl tar make sha256sum; do
    command -v "$tool" >/dev/null || { echo "Required tool missing in the MARS environment: $tool" >&2; exit 1; }
done
task_cc=$(command -v mpicc)
task_cxx=$(command -v mpicxx)
task_cuda=$(command -v nvcc)
mkdir "$1"
task_build=$(cd "$1" && pwd -P)
task_version=3.24.0
task_archive="petsc-$task_version.tar.gz"
curl --fail --location --retry 3 \
    "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/$task_archive" \
    --output "$task_build/$task_archive"
(
    cd "$task_build"
    printf '%s  %s\n' cc9063d80cae3ca87dd34586a92bac49613818a0689d9eac1bd91a799c5d0983 "$task_archive" | sha256sum --check
    tar -xzf "$task_archive"
)
task_source="$task_build/petsc-$task_version"
(
    cd "$task_source"
    PETSC_DIR="$task_source" python3 ./configure PETSC_ARCH=arch-mars-cuda \
        --with-cc="$task_cc" --with-cxx="$task_cxx" --with-fc=0 \
        --with-cuda=1 --with-cudac="$task_cuda" --with-cuda-arch=90 \
        --with-debugging=0 --with-shared-libraries=1 --with-fortran-bindings=0 \
        --with-scalar-type=real --with-precision=double \
        --with-x=0 --with-c2html=0 --with-make-np=4 --download-f2cblaslapack=1 \
        2>&1 | tee "$task_build/configure.log"
    make PETSC_DIR="$task_source" PETSC_ARCH=arch-mars-cuda all \
        2>&1 | tee "$task_build/build.log"
)
printf '\nPETSc built. CMake options:\n-DMARS_PETSC_ROOT=%s\n-DMARS_PETSC_ARCH=arch-mars-cuda\n' "$task_source"
