#pragma once
// A Galerkin operator MARSIR generates from a `form = galerkin` spec
// (marsir-compiler/specs/bp1.op, bp3.op; PTX from marsir-mlir/test/bp_gpu_pipeline.py),
// run on element-local vectors. Kernel contract: {spec}_apply_p{p}(u, [G], y) with u, y
// E x n x n x n (n = p + 1 GLL nodes per axis) and G E x K x q x q x q (the K geometry
// components the spec reads), each an MLIR memref descriptor. The 1D matrices are globals
// in constant memory, {kernel}__m0 .. __m3: B q x n (the GLL basis at the q Gauss points),
// Bt its transpose, and when the spec reads a gradient D, the q x q Gauss-point derivative,
// and Dt its transpose; read at constant indices, they become fma operands instead of
// loads through L1. They are copied in when their pointers change. A `geometry = trilinear`
// spec (bp1_tri.op, bp3_tri.op) takes gq, wq (q: the Gauss points and weights) and X (E x
// 3 x 2 x 2 x 2: the corners by component, then corner bit along each reference axis)
// instead of G, and forms the geometry in registers. The launch shape is the published
// {kernel}__launch, or one warp per element. A matrix with a published {symbol}__fold =
// {rows, cols, sign} is read in the even-odd folded form (marsir-mlir --mir-even-odd,
// mars_marsir_fold.hpp), which is formed here on the device from the matrix given.

#include "backend/distributed/unstructured/marsir/mars_marsir_fold.hpp"
#include "backend/distributed/unstructured/marsir/mars_marsir_ptx_operator.hpp"

#include <string>
#include <vector>

namespace mars {
namespace marsir {

template <typename Real>
__global__ void fold_kernel(const Real* __restrict__ M, Real* __restrict__ out, long long r, long long c)
{
    const long long k = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (k < folded_size(r, c)) out[k] = folded_entry(M, r, c, k);
}

class GalerkinPtx {
public:
    // grad: the kernel takes D and Dt; K: geometry components (0: no G argument).
    GalerkinPtx(const std::string& ptx_path, const std::string& kernel, int p, int q, bool grad, int K)
        : p_(p), q_(q), grad_(grad), K_(K)
    {
        module_ = load_ptx(ptx_path);
        MARS_MARSIR_CU(cuModuleGetFunction(&fn_, module_, kernel.c_str()));
        CUdeviceptr global = 0;
        size_t bytes = 0;
        if (cuModuleGetGlobal(&global, &bytes, module_, (kernel + "__launch").c_str()) == CUDA_SUCCESS &&
            bytes == sizeof(long long[2])) {
            long long v[2];
            MARS_MARSIR_CU(cuMemcpyDtoH(v, global, sizeof(v)));
            threads_ = (unsigned)v[0];
            per_block_ = v[1];
        }
        for (int k = 0; k < (grad ? 4 : 2); ++k) {
            const std::string name = kernel + "__m" + std::to_string(k);
            MARS_MARSIR_CU(cuModuleGetGlobal(&matrix_[k], &matrix_bytes_[k], module_, name.c_str()));
            // Setup only: three integers, like __launch above.
            if (cuModuleGetGlobal(&global, &bytes, module_, (name + "__fold").c_str()) == CUDA_SUCCESS &&
                bytes == sizeof(long long[3])) {
                MARS_MARSIR_CU(cuMemcpyDtoH(fold_[k], global, sizeof(fold_[k])));
                if (folded_size(fold_[k][0], fold_[k][1]) * (long long)sizeof(double) != (long long)matrix_bytes_[k]) {
                    fprintf(stderr, "%s: %zu bytes, not the folded size of its __fold record\n", name.c_str(),
                            matrix_bytes_[k]);
                    std::abort();
                }
                MARS_MARSIR_CU(cuMemAlloc(&folded_[k], matrix_bytes_[k]));
            }
        }
    }
    ~GalerkinPtx()
    {
        for (CUdeviceptr f : folded_)
            if (f) cuMemFree(f);
        cuModuleUnload(module_);
    }
    GalerkinPtx(const GalerkinPtx&) = delete;
    GalerkinPtx& operator=(const GalerkinPtx&) = delete;

    // y = A u on E elements (u continuous, y unassembled), stored geometry G.
    void apply(const double* u, const double* B, const double* Bt, const double* D, const double* Dt,
               const double* G, double* y, long long E, cudaStream_t stream = 0) const
    {
        std::vector<Desc> geometry;
        if (K_ > 0) geometry.push_back(make(G, {E, K_, q_, q_, q_}));
        launch(u, B, Bt, D, Dt, geometry, y, E, stream);
    }

    // The same for a `geometry = trilinear` kernel (constructed with K = 0).
    void apply_trilinear(const double* u, const double* B, const double* Bt, const double* D, const double* Dt,
                         const double* gq, const double* wq, const double* X, double* y, long long E,
                         cudaStream_t stream = 0) const
    {
        launch(u, B, Bt, D, Dt, {make(gq, {q_}), make(wq, {q_}), make(X, {E, 3, 2, 2, 2})}, y, E, stream);
    }

private:
    // The driver API takes a pointer to every scalar of every descriptor, so each
    // descriptor is built in storage that outlives the launch call.
    struct Desc {
        const void* ptr;
        long long offset = 0, sizes[5], strides[5];
        int rank;
    };
    static Desc make(const void* ptr, std::vector<long long> sizes)
    {
        Desc d{ptr};
        d.rank = (int)sizes.size();
        long long stride = 1;
        for (int r = d.rank - 1; r >= 0; --r) {
            d.sizes[r] = sizes[r];
            d.strides[r] = stride;
            stride *= sizes[r];
        }
        return d;
    }

    void launch(const double* u, const double* B, const double* Bt, const double* D, const double* Dt,
                const std::vector<Desc>& geometry, double* y, long long E, cudaStream_t stream) const
    {
        if (E == 0) return;
        const long long n = p_ + 1;
        const double* source[4] = {B, Bt, D, Dt};
        for (int k = 0; k < (grad_ ? 4 : 2); ++k)
            if (source[k] != copied_[k]) {
                CUdeviceptr from = (CUdeviceptr)source[k];
                if (folded_[k]) {
                    const long long size = folded_size(fold_[k][0], fold_[k][1]);
                    fold_kernel<<<(unsigned)((size + 127) / 128), 128, 0, stream>>>(source[k], (double*)folded_[k],
                                                                                    fold_[k][0], fold_[k][1]);
                    if (cudaError_t e = cudaGetLastError(); e != cudaSuccess) {
                        fprintf(stderr, "fold_kernel: %s\n", cudaGetErrorString(e));
                        std::abort();
                    }
                    from = folded_[k];
                }
                MARS_MARSIR_CU(cuMemcpyDtoDAsync(matrix_[k], from, matrix_bytes_[k], (CUstream)stream));
                copied_[k] = source[k];
            }
        std::vector<Desc> descs{make(u, {E, n, n, n})};
        descs.insert(descs.end(), geometry.begin(), geometry.end());
        descs.push_back(make(y, {E, n, n, n}));
        std::vector<void*> params;
        for (Desc& d : descs) {
            params.push_back(&d.ptr);
            params.push_back(&d.ptr);
            params.push_back(&d.offset);
            for (int r = 0; r < d.rank; ++r) params.push_back(&d.sizes[r]);
            for (int r = 0; r < d.rank; ++r) params.push_back(&d.strides[r]);
        }
        MARS_MARSIR_CU(cuLaunchKernel(fn_, (unsigned)((E + per_block_ - 1) / per_block_), 1, 1, threads_, 1, 1, 0,
                                      (CUstream)stream, params.data(), nullptr));
    }

    CUmodule module_ = nullptr;
    CUfunction fn_ = nullptr;
    int p_, q_;
    bool grad_;
    int K_;
    unsigned threads_ = 32;
    long long per_block_ = 1;
    CUdeviceptr matrix_[4] = {};
    size_t matrix_bytes_[4] = {};
    mutable const double* copied_[4] = {};   // the matrices last copied into matrix_
    long long fold_[4][3] = {};               // {rows, cols, sign} of a folded matrix
    CUdeviceptr folded_[4] = {};              // its folded form, before the copy into matrix_
};

}  // namespace marsir
}  // namespace mars
