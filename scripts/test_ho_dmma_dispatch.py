#!/usr/bin/env python3
"""Compile the public launchers against recording stubs; no GPU execution."""
import os
from pathlib import Path
import subprocess
import tempfile

root = Path(__file__).resolve().parents[1]
source = (root / 'backend/distributed/unstructured/fem/mars_ho_laplacian_dmma.hpp').read_text()
source = source[source.index('// Per-quad-point metric (general/curved).'):source.index('// Matched-structure SCALAR baseline:')]
prologue = r'''
#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <type_traits>
using cudaError_t = int;
using cudaStream_t = int;
#ifdef __HIPCC__
constexpr int hoDmmaBlockSize = 64;
#else
constexpr int hoDmmaBlockSize = 128;
#endif
bool matrix_available = true; // FP64 matrix instructions: then the plane kernel, else the staged one
bool hoFp64MatrixAvailable() { return matrix_available; }
const double* expected_u;
double* expected_y;
const int* expected_dof;
const double *expected_d, *expected_g, *expected_w;
int calls = 0;
int record(int mode, const double* u, double* y, const int* dof, const double* d,
           const double* g, const double* w, size_t count, int stream) {
    assert(u == expected_u && y == expected_y && dof == expected_dof);
    assert(d == expected_d && g == expected_g && w == (mode ? expected_w : nullptr));
    assert(count == 3 && stream == 19);
    ++calls;
    return 37; // Check that the API forwards launch errors.
}
#ifdef __HIPCC__
template<int Mode>
int ho_laplacian_mfma_plane_launch(const double* u, double* y, const int* dof,
    const double* d, const double* g, const double* w, size_t count, int stream) {
    assert(matrix_available && "the MFMA kernel needs FP64 MFMA");
    return record(Mode,u,y,dof,d,g,w,count,stream);
}
#else
template<int Mode>
int ho_laplacian_dmma_plane_launch(const double* u, double* y, const int* dof,
    const double* d, const double* g, const double* w, size_t count, int stream) {
    assert(matrix_available && "the plane kernel needs SM80+");
    return record(Mode,u,y,dof,d,g,w,count,stream);
}
#endif
template<class T, int P, int Block, int Mode, int Tile>
int ho_laplacian_dmma_launch_impl(const T* u, T* y, const int* dof,
    const T* d, const T* g, const T* w, size_t count, int stream) {
    assert(!matrix_available && "with FP64 matrix instructions the plane kernel runs");
    static_assert(P == 7 && Block == hoDmmaBlockSize && Tile == 0);
    return record(Mode,u,y,dof,d,g,w,count,stream);
}
'''
main = r'''
int main() {
    double u[1],y[1],d[1],g[1],w[1]; int dof[1];
    expected_u=u;expected_y=y;expected_dof=dof;expected_d=d;expected_g=g;expected_w=w;
    for (bool available : {true, false}) {
        matrix_available = available;
        assert((ho_laplacian_dmma_launch<double,7>(u,y,dof,d,g,3,19)==37));
        assert((ho_laplacian_dmma_launch_affine<double,7>(u,y,dof,d,g,w,3,19)==37));
    }
    assert(calls==4);
}
'''
with tempfile.TemporaryDirectory(prefix='mars-dmma-dispatch-') as tmp:
    cpp = Path(tmp) / 'dispatch.cpp'
    exe = Path(tmp) / 'dispatch'
    cpp.write_text(prologue + source + main)
    for label, flags in [('CUDA', []), ('HIP', ['-D__HIPCC__'])]:
        subprocess.run([os.environ.get('CXX', 'c++'), '-std=c++20', '-Wall', '-Wextra',
                        *flags, str(cpp), '-o', str(exe)], check=True)
        subprocess.run([str(exe)], check=True)
        print(f'PASS {label}: public dispatch, arguments, non-default stream and error forwarding')
