#pragma once
// A Galerkin operator MARSIR generates from a `form = galerkin` spec
// (marsir-compiler/specs/bp1.op, bp3.op; PTX from marsir-mlir/test/bp_gpu_pipeline.py),
// run on element-local vectors. Kernel contract: {spec}_apply_p{p}(u, B, Bt, [D, Dt],
// [G], y) with u, y E x n x n x n (n = p + 1 GLL nodes per axis), B q x n (the GLL basis
// at the q Gauss points), Bt its transpose, D the q x q Gauss-point derivative and Dt its
// transpose (when the spec reads a gradient), G E x K x q x q x q (the K geometry
// components the spec reads), each an MLIR memref descriptor. A `geometry = trilinear`
// spec (bp1_tri.op, bp3_tri.op) takes gq, wq (q: the Gauss points and weights) and X (E x
// 3 x 2 x 2 x 2: the corners by component, then corner bit along each reference axis)
// instead of G, and forms the geometry in registers. The launch shape is the published
// {kernel}__launch, or one warp per element.

#include "backend/distributed/unstructured/marsir/mars_marsir_ptx_operator.hpp"

#include <string>
#include <vector>

namespace mars {
namespace marsir {

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
    }
    ~GalerkinPtx() { cuModuleUnload(module_); }
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
        const long long n = p_ + 1, q = q_;
        std::vector<Desc> descs{make(u, {E, n, n, n}), make(B, {q, n}), make(Bt, {n, q})};
        if (grad_) {
            descs.push_back(make(D, {q, q}));
            descs.push_back(make(Dt, {q, q}));
        }
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
};

}  // namespace marsir
}  // namespace mars
