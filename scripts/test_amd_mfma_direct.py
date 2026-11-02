#!/usr/bin/env python3
"""Host lane emulation of the actual MFMA kernel, not a GPU test."""
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
source = (ROOT / 'backend/distributed/unstructured/fem/mars_ho_laplacian_mfma_direct.hpp').read_text()
# Keep the kernel and its operand/address calculations; omit HIP headers and the launch syntax.
source = source[source.index('namespace mars::fem'):source.index('// Asynchronous additive apply')]
source += '\n#endif\n} // namespace mars::fem\n'
prologue = r'''
#include <algorithm>
#include <atomic>
#include <barrier>
#include <cmath>
#include <cstdio>
#include <thread>
#include <vector>
#define __host__
#define __device__
#define __global__
#define __forceinline__ inline
#define __launch_bounds__(...)
#define __shared__ static
#define __HIPCC__ 1
#define MARS_HAVE_FP64_MFMA 1
struct Dim { unsigned x; };
thread_local Dim threadIdx{};
Dim blockIdx{};
std::barrier wave(64);
double operand_a[64],operand_b[64];
void syncWarp() { wave.arrive_and_wait(); }
void atomicAdd(double* p,double v) { std::atomic_ref<double>(*p).fetch_add(v); }
double emulate_mfma(double a,double b,double acc,int,int,int) {
    int lane=threadIdx.x, row=lane/16, block=(lane%16)/4,col=lane%4;
    operand_a[lane]=a;operand_b[lane]=b;
    wave.arrive_and_wait();
    for(int k=0;k<4;++k) acc+=operand_a[16*k+4*block+row]*operand_b[16*k+4*block+col];
    wave.arrive_and_wait();
    return acc;
}
#define __builtin_amdgcn_mfma_f64_4x4x4f64 emulate_mfma
'''
with tempfile.TemporaryDirectory(prefix='mars-mfma-host-') as tmp:
    cpp = Path(tmp) / 'check.cpp'
    exe = Path(tmp) / 'check'
    cpp.write_text(prologue + source + (ROOT / 'tests/performance/mfma_direct_host.cpp').read_text())
    subprocess.run([os.environ.get('CXX', 'c++'), '-std=c++20', '-O1', '-pthread', '-I', str(ROOT / 'tests/performance'),
                    '-fsanitize=address,undefined', str(cpp), '-o', str(exe)], check=True)
    subprocess.run([str(exe)], check=True)
