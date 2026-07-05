// Batched GPU gate + throughput measurement for the mir-generated DMMA PTX
// (sumfac_batched_sm90.ptx): Y[e] += D(8x8) * U[e](8x64), one warp per element,
// grid = E blocks. Self-contained (dlopen libcuda.so.1, no CUDA toolkit).
//
// ABI: 25 params = D static memref<8x8> (7 scalars) + U,Y dynamic
// memref<?x8x64> (9 scalars each: alloc, aligned, offset, 3 sizes, 3 strides).
//
// MEASUREMENT NOTE: this times ONE sum-factorization sweep, not the full HO
// Laplacian apply -- numbers are not directly comparable to mars_ho_dmma_bench
// (which does the whole operator). Reported: ms/apply, Gelem/s, and effective
// GB/s from the actual traffic (per element: read U 4KiB + read Y 4KiB +
// write Y 4KiB = 12 KiB; D is negligible/broadcast).
//
// Build:  g++ -O2 run_ptx_batched.cpp -o run_ptx_batched -ldl
// Run:    ./run_ptx_batched ../generated/sumfac_batched_sm90.ptx [numElements]

#include <dlfcn.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

typedef int CUresult;
typedef int CUdevice;
typedef unsigned long long CUdeviceptr;
typedef struct CUctx_st* CUcontext;
typedef struct CUmod_st* CUmodule;
typedef struct CUfunc_st* CUfunction;
typedef struct CUstream_st* CUstream;
typedef struct CUevent_st* CUevent;

static CUresult (*p_cuInit)(unsigned);
static CUresult (*p_cuDeviceGet)(CUdevice*, int);
static CUresult (*p_cuCtxCreate)(CUcontext*, unsigned, CUdevice);
static CUresult (*p_cuModuleLoadData)(CUmodule*, const void*);
static CUresult (*p_cuModuleGetFunction)(CUfunction*, CUmodule, const char*);
static CUresult (*p_cuMemAlloc)(CUdeviceptr*, size_t);
static CUresult (*p_cuMemcpyHtoD)(CUdeviceptr, const void*, size_t);
static CUresult (*p_cuMemcpyDtoH)(void*, CUdeviceptr, size_t);
static CUresult (*p_cuLaunchKernel)(CUfunction, unsigned, unsigned, unsigned,
                                    unsigned, unsigned, unsigned, unsigned,
                                    CUstream, void**, void**);
static CUresult (*p_cuCtxSynchronize)(void);
static CUresult (*p_cuEventCreate)(CUevent*, unsigned);
static CUresult (*p_cuEventRecord)(CUevent, CUstream);
static CUresult (*p_cuEventSynchronize)(CUevent);
static CUresult (*p_cuEventElapsedTime)(float*, CUevent, CUevent);
static CUresult (*p_cuGetErrorString)(CUresult, const char**);
static CUresult (*p_cuMemsetD8)(CUdeviceptr, unsigned char, size_t);

#define CK(x) do { CUresult r_ = (x); if (r_ != 0) { \
    const char* s_ = nullptr; if (p_cuGetErrorString) p_cuGetErrorString(r_, &s_); \
    fprintf(stderr, "CUDA error %d (%s) at %s:%d\n", r_, s_ ? s_ : "?", __FILE__, __LINE__); \
    exit(1); } } while (0)

static void* must_sym(void* lib, const char* name)
{
    void* p = dlsym(lib, name);
    if (!p) { fprintf(stderr, "missing driver symbol %s\n", name); exit(1); }
    return p;
}

int main(int argc, char** argv)
{
    void* lib = dlopen("libcuda.so.1", RTLD_NOW);
    if (!lib) lib = dlopen("libcuda.so", RTLD_NOW);
    if (!lib) { fprintf(stderr, "cannot dlopen libcuda.so.1: %s\n", dlerror()); return 1; }
    *(void**)&p_cuInit              = must_sym(lib, "cuInit");
    *(void**)&p_cuDeviceGet         = must_sym(lib, "cuDeviceGet");
    *(void**)&p_cuCtxCreate         = must_sym(lib, "cuCtxCreate_v2");
    *(void**)&p_cuModuleLoadData    = must_sym(lib, "cuModuleLoadData");
    *(void**)&p_cuModuleGetFunction = must_sym(lib, "cuModuleGetFunction");
    *(void**)&p_cuMemAlloc          = must_sym(lib, "cuMemAlloc_v2");
    *(void**)&p_cuMemcpyHtoD        = must_sym(lib, "cuMemcpyHtoD_v2");
    *(void**)&p_cuMemcpyDtoH        = must_sym(lib, "cuMemcpyDtoH_v2");
    *(void**)&p_cuLaunchKernel      = must_sym(lib, "cuLaunchKernel");
    *(void**)&p_cuCtxSynchronize    = must_sym(lib, "cuCtxSynchronize");
    *(void**)&p_cuEventCreate       = must_sym(lib, "cuEventCreate");
    *(void**)&p_cuEventRecord       = must_sym(lib, "cuEventRecord");
    *(void**)&p_cuEventSynchronize  = must_sym(lib, "cuEventSynchronize");
    *(void**)&p_cuEventElapsedTime  = must_sym(lib, "cuEventElapsedTime");
    *(void**)&p_cuGetErrorString    = must_sym(lib, "cuGetErrorString");
    *(void**)&p_cuMemsetD8          = must_sym(lib, "cuMemsetD8_v2");

    const char* ptxPath = argc > 1 ? argv[1] : "../generated/sumfac_batched_sm90.ptx";
    const long long E = argc > 2 ? atoll(argv[2]) : (1LL << 20);   // elements
    const char* kname = argc > 3 ? argv[3] : "sumfac_sweep_batched";
    const int wpb = argc > 4 ? atoi(argv[4]) : 1;   // warps (elements) per block
    if (E % wpb) { fprintf(stderr, "E must be divisible by warps-per-block\n"); return 1; }
    const unsigned grid = (unsigned)(E / wpb);

    FILE* f = fopen(ptxPath, "rb");
    if (!f) { perror(ptxPath); return 1; }
    fseek(f, 0, SEEK_END); long len = ftell(f); fseek(f, 0, SEEK_SET);
    std::vector<char> ptx(len + 1, 0);
    if (fread(ptx.data(), 1, len, f) != (size_t)len) { fclose(f); return 1; }
    fclose(f);

    CK(p_cuInit(0));
    CUdevice dev; CK(p_cuDeviceGet(&dev, 0));
    CUcontext ctx; CK(p_cuCtxCreate(&ctx, 0, dev));
    CUmodule mod; CK(p_cuModuleLoadData(&mod, ptx.data()));
    CUfunction fn; CK(p_cuModuleGetFunction(&fn, mod, kname));

    const int M = 8, K = 8, N = 64;
    const size_t elemDoubles = (size_t)K * N;                 // 512 doubles
    const size_t totalBytes = (size_t)E * elemDoubles * 8;

    std::vector<double> hD(M * K);
    std::vector<double> hU((size_t)E * elemDoubles);
    srand(42);
    for (auto& v : hD) v = 2.0 * rand() / RAND_MAX - 1.0;
    for (auto& v : hU) v = 2.0 * rand() / RAND_MAX - 1.0;

    CUdeviceptr dD, dU, dY;
    CK(p_cuMemAlloc(&dD, M * K * sizeof(double)));
    CK(p_cuMemAlloc(&dU, totalBytes));
    CK(p_cuMemAlloc(&dY, totalBytes));
    CK(p_cuMemcpyHtoD(dD, hD.data(), M * K * sizeof(double)));
    CK(p_cuMemcpyHtoD(dU, hU.data(), totalBytes));
    CK(p_cuMemsetD8(dY, 0, totalBytes));

    long long zero = 0;
    long long dSz[2] = {8, 8}, dSt[2] = {8, 1};
    long long bSz[3] = {E, 8, 64}, bSt[3] = {512, 64, 1};
    void* args[25] = {
        &dD, &dD, &zero, &dSz[0], &dSz[1], &dSt[0], &dSt[1],
        &dU, &dU, &zero, &bSz[0], &bSz[1], &bSz[2], &bSt[0], &bSt[1], &bSt[2],
        &dY, &dY, &zero, &bSz[0], &bSz[1], &bSz[2], &bSt[0], &bSt[1], &bSt[2],
    };
    // v4 (outlined "sweep_kernel"): gpu-kernel-outlining did not sink the loop
    // constants, so they arrive as params, in launch_func operand order:
    // [c1, c0][U desc 9][Y desc 9][D desc 7][vector<1x1xf64> zero pad][c0,c8,c4,c64]
    long long i1 = 1, i0 = 0, i8v = 8, i4v = 4, i64v = 64;
    double vpad = 0.0;
    // v5 ("sweepv5_kernel"): [c1,c0,c8][U 9][Y 9][D 7][vec pad][c4,c64]; launch grid=E/8, block=(32,8).
    void* argsV5[31] = {
        &i1, &i0, &i8v,
        &dU, &dU, &zero, &bSz[0], &bSz[1], &bSz[2], &bSt[0], &bSt[1], &bSt[2],
        &dY, &dY, &zero, &bSz[0], &bSz[1], &bSz[2], &bSt[0], &bSt[1], &bSt[2],
        &dD, &dD, &zero, &dSz[0], &dSz[1], &dSt[0], &dSt[1],
        &vpad,
        &i4v, &i64v,
    };
    void* argsV4[32] = {
        &i1, &i0,
        &dU, &dU, &zero, &bSz[0], &bSz[1], &bSz[2], &bSt[0], &bSt[1], &bSt[2],
        &dY, &dY, &zero, &bSz[0], &bSz[1], &bSz[2], &bSt[0], &bSt[1], &bSt[2],
        &dD, &dD, &zero, &dSz[0], &dSz[1], &dSt[0], &dSt[1],
        &vpad,
        &i0, &i8v, &i4v, &i64v,
    };
    void** kargs = args;
    if (strcmp(kname, "sweep_kernel") == 0) kargs = argsV4;
    if (strcmp(kname, "sweepv5_kernel") == 0) kargs = argsV5;

    // --- correctness: one apply into zeroed Y, spot-check 8 elements ---
    CK(p_cuLaunchKernel(fn, grid, 1, 1, 32, (unsigned)wpb, 1, 0, 0, kargs, nullptr));
    CK(p_cuCtxSynchronize());
    double err = 0;
    std::vector<double> out(elemDoubles), ref(elemDoubles);
    for (int probe = 0; probe < 8; ++probe) {
        long long e = (long long)((double)rand() / RAND_MAX * (E - 1));
        CK(p_cuMemcpyDtoH(out.data(), dY + (size_t)e * elemDoubles * 8,
                          elemDoubles * sizeof(double)));
        const double* u = &hU[(size_t)e * elemDoubles];
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < N; ++j) {
                double s = 0;
                for (int p = 0; p < K; ++p) s += hD[i * K + p] * u[p * N + j];
                ref[i * N + j] = s;
            }
        for (size_t i = 0; i < elemDoubles; ++i)
            err = fmax(err, fabs(out[i] - ref[i]));
    }
    printf("batched DMMA gate (E=%lld, %s, %d warps/block): spot max|Y - D*U| = %.3e  %s\n",
           E, kname, wpb, err, err < 1e-12 ? "PASS" : "FAIL");

    // --- throughput: time full-grid applies (Y keeps accumulating; fine for timing) ---
    CUevent t0, t1; CK(p_cuEventCreate(&t0, 0)); CK(p_cuEventCreate(&t1, 0));
    for (int w = 0; w < 3; ++w)
        CK(p_cuLaunchKernel(fn, grid, 1, 1, 32, (unsigned)wpb, 1, 0, 0, kargs, nullptr));
    CK(p_cuCtxSynchronize());
    const int iters = 50;
    CK(p_cuEventRecord(t0, 0));
    for (int i = 0; i < iters; ++i)
        CK(p_cuLaunchKernel(fn, grid, 1, 1, 32, (unsigned)wpb, 1, 0, 0, kargs, nullptr));
    CK(p_cuEventRecord(t1, 0)); CK(p_cuEventSynchronize(t1));
    float ms = 0; CK(p_cuEventElapsedTime(&ms, t0, t1));
    const double msPerApply = ms / iters;
    // Traffic per element: read U (4 KiB) + read Y + write Y (accumulate) = 12 KiB.
    const double gbPerApply = (double)E * 12.0 * 1024.0 / 1e9;
    printf("  %d applies: %.3f ms/apply | %.2f Gelem/s | ~%.0f GB/s effective (12KiB/elem)\n",
           iters, msPerApply, (double)E / (msPerApply * 1e-3) / 1e9,
           gbPerApply / (msPerApply * 1e-3));
    printf("  NOTE: one sweep only -- not comparable 1:1 to the full-Laplacian mars_ho_dmma_bench.\n");
    return err < 1e-12 ? 0 : 1;
}
