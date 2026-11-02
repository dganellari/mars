#!/usr/bin/env python3
"""Host lane emulation of the CUDA DMMA plane kernel, not a GPU test."""
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
source = (ROOT / 'backend/distributed/unstructured/fem/mars_ho_laplacian_dmma_plane.hpp').read_text()
# Keep the kernel and its operand/address calculations; omit the CUDA headers. With MARS_DMMA_PLANE_DEVICE
# the header compiles the kernel body on the host and leaves out the launcher.
source = source[source.index('#if defined(MARS_DMMA_PLANE_DEVICE) ||'):]
prologue = r'''
#include <algorithm>
#include <atomic>
#include <barrier>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <thread>
#include <vector>
#define __host__
#define __device__
#define __global__
#define __forceinline__ inline
#define __launch_bounds__(...)
#define __shared__ static
#define MARS_DMMA_PLANE_DEVICE 1
struct Dim { unsigned x; };
thread_local Dim threadIdx{};
Dim blockIdx{};
std::barrier warp(32);
double operand_a[32],operand_b[32];
void __syncwarp() { warp.arrive_and_wait(); }
void __trap() { std::abort(); }
void atomicAdd(double* p,double v) { std::atomic_ref<double>(*p).fetch_add(v); }
// m8n8k4 f64 row.col, as verified on GH200 (marsir-compiler/marsir/backends/mlir_warp_reg.py): lane L holds
// A[L/4][L%4], B[L%4][L/4], C[L/4][2(L%4)] and C[L/4][2(L%4)+1].
void emulate_dmma(double a,double b,double& c0,double& c1) {
    const int lane=threadIdx.x%32, i=lane/4, q=lane%4;
    operand_a[lane]=a;operand_b[lane]=b;
    warp.arrive_and_wait();
    for(int k=0;k<4;++k) {
        c0+=operand_a[4*i+k]*operand_b[4*(2*q)+k];
        c1+=operand_a[4*i+k]*operand_b[4*(2*q+1)+k];
    }
    warp.arrive_and_wait();
}
'''
with tempfile.TemporaryDirectory(prefix='mars-dmma-plane-host-') as tmp:
    cpp = Path(tmp) / 'check.cpp'
    exe = Path(tmp) / 'check'
    cpp.write_text(prologue + source + (ROOT / 'tests/performance/dmma_plane_host.cpp').read_text())
    subprocess.run([os.environ.get('CXX', 'c++'), '-std=c++20', '-O1', '-pthread', '-I', str(ROOT / 'tests/performance'),
                    '-fsanitize=address,undefined', str(cpp), '-o', str(exe)], check=True)
    subprocess.run([str(exe)], check=True)
