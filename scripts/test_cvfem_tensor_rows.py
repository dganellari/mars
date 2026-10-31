#!/usr/bin/env python3
"""Host emulation of the row map build and the tensor_rows kernels, not a GPU test: one std::thread per GPU
thread, a std::barrier per block. Checks tensor_rows and tensor_map entry by entry against tensor_faces."""
import os
from pathlib import Path
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
FEM = ROOT / 'backend/distributed/unstructured/fem'


def body(name):
    text = (FEM / name).read_text()
    text = re.sub(r'^\s*#\s*(include|pragma once).*$', '', text, flags=re.M)
    # kernel<<<grid, block[, shared, stream]>>>(args) -> hostLaunch(grid, block, call kernel, args); the
    # generic lambda lets the call deduce the kernel's template arguments, as a launch does
    def launch(m):
        depth, parts, cur = 0, [], ''
        for ch in m.group(2):
            depth += ch in '(<'
            depth -= ch in ')>'
            if ch == ',' and depth == 0:
                parts.append(cur)
                cur = ''
            else:
                cur += ch
        parts.append(cur)
        return f'hostLaunch({parts[0]}, {parts[1]}, [](auto... a) {{ {m.group(1)}(a...); }}, '
    return re.sub(r'([\w:]+(?:<[^<>;]*>)?)<<<(.+?)>>>\(', launch, text)


prologue = r'''
#include <algorithm>
#include <array>
#include <atomic>
#include <barrier>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <vector>
#define __host__
#define __device__
#define __global__
#define __constant__
#define __forceinline__ inline
#define __launch_bounds__(...)
#define __shared__ static
struct Dim { unsigned x = 0; };
struct alignas(16) uint4 { unsigned x, y, z, w; };
thread_local Dim threadIdx, blockIdx, blockDim, gridDim;
thread_local std::barrier<>* blockBarrier;
void __syncthreads() { blockBarrier->arrive_and_wait(); }
// stronger than a warp barrier, which is fine: every lane of a block calls it
void syncWarp() { blockBarrier->arrive_and_wait(); }
template<class K, class... A>
void hostLaunch(size_t grid, size_t block, K kernel, A... args)
{
    for (size_t b = 0; b < grid; ++b)
    {
        std::barrier<> barrier{std::ptrdiff_t(block)};
        std::vector<std::thread> threads;
        for (size_t t = 0; t < block; ++t)
            threads.emplace_back([&, t] {
                blockIdx.x = unsigned(b), threadIdx.x = unsigned(t), blockDim.x = unsigned(block);
                gridDim.x = unsigned(grid), blockBarrier = &barrier;
                kernel(args...);
            });
        for (auto& t : threads) t.join();
    }
}
double atomicAdd(double* p, double v) { return std::atomic_ref<double>(*p).fetch_add(v); }
int atomicAdd(int* p, int v) { return std::atomic_ref<int>(*p).fetch_add(v); }
int atomicMax(int* p, int v)
{
    std::atomic_ref<int> a(*p);
    int old = a.load();
    while (old < v && !a.compare_exchange_weak(old, v)) {}
    return old;
}
int atomicCAS(int* p, int expected, int desired)
{
    std::atomic_ref<int>(*p).compare_exchange_strong(expected, desired);
    return expected;
}
typedef int cudaError_t;
constexpr int cudaSuccess = 0;
enum cudaMemcpyKind { cudaMemcpyDeviceToHost };
inline int cudaMemcpy(void* d, const void* s, size_t n, cudaMemcpyKind) { std::memcpy(d, s, n); return 0; }
inline int cudaDeviceSynchronize() { return 0; }
#define cudaCheckCall(call) (void)(call)
#define cudaCheckError()
namespace cstone {
template<class T>
struct DeviceVector
{
    std::vector<T> v;
    DeviceVector() = default;
    explicit DeviceVector(size_t n) : v(n) {}
    DeviceVector(size_t n, T init) : v(n, init) {}
    T* data() { return v.data(); }
    const T* data() const { return v.data(); }
    void resize(size_t n) { v.resize(n); }
    size_t size() const { return v.size(); }
};
}
namespace thrust {
template<class T>
struct device_vector : std::vector<T> { using std::vector<T>::vector; };
template<class T>
T* raw_pointer_cast(T* p) { return p; }
struct Policy {};
constexpr Policy device{};
template<class It, class Out>
void exclusive_scan(Policy, It first, It last, Out out) { std::exclusive_scan(first, last, out, 0); }
}
'''

source = prologue + ''.join(body(f) for f in ['mars_cvfem_hex_kernel.hpp', 'mars_cvfem_hex_kernel_tensor.hpp',
                                               'mars_cvfem_hex_kernel_tensor_faces.hpp',
                                               'mars_cvfem_hex_kernel_tensor_rows.hpp', 'mars_cvfem_row_map.hpp'])
with tempfile.TemporaryDirectory(prefix='mars-tensor-rows-host-') as tmp:
    cpp = Path(tmp) / 'check.cpp'
    exe = Path(tmp) / 'check'
    cpp.write_text(source + (ROOT / 'tests/performance/tensor_rows_host.cpp').read_text())
    subprocess.run([os.environ.get('CXX', 'c++'), '-std=c++20', '-O1', '-pthread', '-fsanitize=address,undefined',
                    str(cpp), '-o', str(exe)], check=True)
    subprocess.run([str(exe)], check=True)
