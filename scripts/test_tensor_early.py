#!/usr/bin/env python3
"""Run exact kernel bodies sequentially under host sanitizers; GPU validation is separate."""
import os
import sys
from pathlib import Path
import subprocess
import tempfile

root = Path(__file__).resolve().parents[1]
fem = root / 'backend/distributed/unstructured/fem'
def clean(source):
    return '\n'.join(line for line in source.splitlines()
                     if not line.startswith(('#include', '#pragma once')))
base = (fem / 'mars_cvfem_hex_kernel.hpp').read_text().split('// Subdivide hex8')[0]
source = clean(base) + '\n}}\n'
candidate = 'early'
for name in ('mars_cvfem_hex_kernel_tensor.hpp', f'mars_cvfem_hex_kernel_tensor_{candidate}.hpp'):
    source += clean((fem / name).read_text()) + '\n'
prologue = """
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#define __host__
#define __device__
#define __global__
#define __constant__
#define __forceinline__ inline
#define MARS_TENSOR_EARLY_HOST_CHECK 1
struct Dim { unsigned x; };
Dim blockIdx{}, blockDim{}, threadIdx{};
inline void atomicAdd(double* p, double v) { *p += v; }
"""
with tempfile.TemporaryDirectory(prefix='mars-tensor-early-') as tmp:
    cpp = Path(tmp) / 'check.cpp'
    exe = Path(tmp) / 'check'
    cpp.write_text(prologue + source + (root / 'examples/distributed/unstructured/mars_cvfem_tensor_early_check.cu').read_text())
    subprocess.run([os.environ.get('CXX', 'c++'), '-std=c++20', '-O1', '-Wno-unknown-pragmas',
                    '-fsanitize=address,undefined', str(cpp), '-o', str(exe)], check=True)
    subprocess.run([str(exe)], check=True)

# The host kernel harness erases runtime calls; compile the actual stream selection separately.
gate = (root / 'examples/distributed/unstructured/mars_cvfem_tensor_early_check.cu').read_text()
selection = gate[gate.index('#if defined(__HIPCC__)', gate.index('Stream stream{};')):]
selection = selection[:selection.index('#endif') + len('#endif')]
with tempfile.TemporaryDirectory(prefix='mars-tensor-stream-') as tmp:
    cpp = Path(tmp) / 'stream.cpp'
    for platform, prefix, defines in [('HIP', 'hip', ['-D__HIPCC__']), ('CUDA', 'cuda', [])]:
        cpp.write_text('#define CHECK_GPU(call) (void)(call)\n'
                       'using Stream = void*;\n'
                       f'int {prefix}StreamCreateWithFlags(Stream*, unsigned);\n'
                       f'constexpr unsigned {prefix}StreamNonBlocking = 1;\n'
                       'void check() { Stream stream{};\n' + selection + '\n}\n')
        subprocess.run([os.environ.get('CXX', 'c++'), '-std=c++20', '-fsyntax-only',
                        *defines, str(cpp)], check=True)
        print(f'PASS {platform}: exact stream selection with only platform-native declarations')
