#pragma once
// The operator kernel MARSIR generates from high-level MLIR (a PTX file such as
// marsir-mlir/generated/hl_full_p7_sm90.ptx), run on element-local vectors. The PTX
// is JIT-loaded through the CUDA driver API into the context the runtime API has
// already made current, so it mixes freely with runtime-API allocations and kernels.
//
// Kernel contract (see marsir-mlir/test/hl_gpu_pipeline.py): laplacian_apply(U,
// Btil, Dtil, W, D, G, Y), one warp per element (grid = E, block = 32). Every memref
// argument is passed as an MLIR descriptor (allocated ptr, aligned ptr, offset,
// sizes, strides) and must be 16-byte aligned, which cudaMalloc guarantees.

#include <cuda.h>
#include <cuda_runtime.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

namespace mars {
namespace marsir {

#define MARS_MARSIR_CU(call)                                                     \
    do {                                                                         \
        CUresult r_ = (call);                                                    \
        if (r_ != CUDA_SUCCESS) {                                                \
            const char* s_ = nullptr;                                            \
            cuGetErrorString(r_, &s_);                                           \
            fprintf(stderr, "CUDA driver error %d (%s) at %s:%d\n", (int)r_,     \
                    s_ ? s_ : "?", __FILE__, __LINE__);                          \
            std::abort();                                                        \
        }                                                                        \
    } while (0)

class LaplacianPtx {
public:
    // p = 7: U, Y are E x 8 x 8 x 8; Btil, Dtil are 8 x 8 with a zero last row (the
    // face dimension padded to a full tensor-core tile); W, D are 8 x 8; G is
    // E x [3][7][3][8][8], component-major.
    explicit LaplacianPtx(const std::string& ptx_path, int p = 7) : p_(p)
    {
        std::ifstream in(ptx_path, std::ios::binary);
        if (!in) {
            fprintf(stderr, "cannot open MARSIR PTX %s\n", ptx_path.c_str());
            std::abort();
        }
        std::stringstream ss;
        ss << in.rdbuf();
        const std::string ptx = ss.str();
        cudaFree(nullptr);   // make sure the runtime's primary context exists
        MARS_MARSIR_CU(cuInit(0));
        MARS_MARSIR_CU(cuModuleLoadData(&module_, ptx.c_str()));
        MARS_MARSIR_CU(cuModuleGetFunction(&fn_, module_, "laplacian_apply"));
    }
    ~LaplacianPtx() { cuModuleUnload(module_); }
    LaplacianPtx(const LaplacianPtx&) = delete;
    LaplacianPtx& operator=(const LaplacianPtx&) = delete;

    // y = A u on E elements: u continuous, y unassembled.
    void apply(const double* d_u, const double* d_btil, const double* d_dtil, const double* d_w,
               const double* d_d, const double* d_g, double* d_y, long long E,
               cudaStream_t stream = 0) const
    {
        const long long n = p_ + 1, nn = n * n, n3 = nn * n;
        const long long g_elem = 3LL * p_ * 3 * nn;
        long long zero = 0;
        long long u_sz[4] = {E, n, n, n}, u_st[4] = {n3, nn, n, 1};
        long long m_sz[2] = {n, n}, m_st[2] = {n, 1};
        long long g_sz[6] = {E, 3, p_, 3, n, n};
        long long g_st[6] = {g_elem, (long long)p_ * 3 * nn, 3 * nn, nn, n, 1};
        // The driver API takes pointers to each argument value, so the device
        // pointers need addressable copies.
        const double *u = d_u, *bt = d_btil, *dt = d_dtil, *w = d_w, *d = d_d, *g = d_g;
        double* y = d_y;
        void* args[] = {
            &u,  &u,  &zero, &u_sz[0], &u_sz[1], &u_sz[2], &u_sz[3],
                             &u_st[0], &u_st[1], &u_st[2], &u_st[3],
            &bt, &bt, &zero, &m_sz[0], &m_sz[1], &m_st[0], &m_st[1],
            &dt, &dt, &zero, &m_sz[0], &m_sz[1], &m_st[0], &m_st[1],
            &w,  &w,  &zero, &m_sz[0], &m_sz[1], &m_st[0], &m_st[1],
            &d,  &d,  &zero, &m_sz[0], &m_sz[1], &m_st[0], &m_st[1],
            &g,  &g,  &zero, &g_sz[0], &g_sz[1], &g_sz[2], &g_sz[3], &g_sz[4], &g_sz[5],
                             &g_st[0], &g_st[1], &g_st[2], &g_st[3], &g_st[4], &g_st[5],
            &y,  &y,  &zero, &u_sz[0], &u_sz[1], &u_sz[2], &u_sz[3],
                             &u_st[0], &u_st[1], &u_st[2], &u_st[3],
        };
        MARS_MARSIR_CU(cuLaunchKernel(fn_, (unsigned)E, 1, 1, 32, 1, 1, 0, (CUstream)stream,
                                      args, nullptr));
    }

private:
    CUmodule module_ = nullptr;
    CUfunction fn_ = nullptr;
    int p_;
};

}  // namespace marsir
}  // namespace mars
