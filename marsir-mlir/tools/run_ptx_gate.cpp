// GPU execution gate for the mir-generated FP64 DMMA PTX (sumfac_sm90.ptx).
//
// SELF-CONTAINED: needs NO CUDA toolkit at build time. The CUDA *driver* API
// lives in libcuda.so.1 (installed with the GPU driver, always present on GPU
// nodes); cuda.h is only prototypes, so we declare the handful we use and
// dlopen the driver. IMPORTANT: dlsym must request the *_v2 ABI names that
// cuda.h normally maps via macros (the unsuffixed symbols are the legacy
// 32-bit-era ABI with different signatures).
//
// The kernel: Y += D(8x8) * U(8x64), the sum-factorization sweep, compiled by
// MLIR to FP64 tensor-core mma.sync.m8n8k4. ABI notes:
//  - MLIR lowers each memref<...> kernel arg to 7 scalars {alloc, aligned,
//    offset, size0, size1, stride0, stride1} -> 3 memrefs = 21 params.
//  - Warp-cooperative + accumulating: launch EXACTLY one warp (grid 1, block 32).
//
// Build:  g++ -O2 run_ptx_gate.cpp -o run_ptx_gate -ldl
// Run:    ./run_ptx_gate ../generated/sumfac_sm90.ptx

#include <dlfcn.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

// --- minimal CUDA driver API surface (see header comment re: _v2) ---
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

    const char* ptxPath = argc > 1 ? argv[1] : "../generated/sumfac_sm90.ptx";
    FILE* f = fopen(ptxPath, "rb");
    if (!f) { perror(ptxPath); return 1; }
    fseek(f, 0, SEEK_END); long len = ftell(f); fseek(f, 0, SEEK_SET);
    std::vector<char> ptx(len + 1, 0);
    if (fread(ptx.data(), 1, len, f) != (size_t)len) { fclose(f); return 1; }
    fclose(f);

    CK(p_cuInit(0));
    CUdevice dev; CK(p_cuDeviceGet(&dev, 0));
    CUcontext ctx; CK(p_cuCtxCreate(&ctx, 0, dev));
    CUmodule mod; CK(p_cuModuleLoadData(&mod, ptx.data()));   // driver JITs PTX->SASS
    CUfunction fn; CK(p_cuModuleGetFunction(&fn, mod, "sumfac_sweep"));

    const int M = 8, K = 8, N = 64;
    std::vector<double> hD(M * K), hU(K * N), hY(M * N, 0.0), ref(M * N, 0.0), out(M * N);
    srand(42);
    for (auto& v : hD) v = 2.0 * rand() / RAND_MAX - 1.0;
    for (auto& v : hU) v = 2.0 * rand() / RAND_MAX - 1.0;
    for (int i = 0; i < M; ++i)
        for (int p = 0; p < K; ++p) {
            const double d = hD[i * K + p];
            for (int j = 0; j < N; ++j) ref[i * N + j] += d * hU[p * N + j];
        }

    CUdeviceptr dD, dU, dY;
    CK(p_cuMemAlloc(&dD, M * K * sizeof(double)));
    CK(p_cuMemAlloc(&dU, K * N * sizeof(double)));
    CK(p_cuMemAlloc(&dY, M * N * sizeof(double)));
    CK(p_cuMemcpyHtoD(dD, hD.data(), M * K * sizeof(double)));
    CK(p_cuMemcpyHtoD(dU, hU.data(), K * N * sizeof(double)));
    CK(p_cuMemcpyHtoD(dY, hY.data(), M * N * sizeof(double)));

    long long zero = 0;
    long long dSz[2] = {8, 8},  dSt[2] = {8, 1};
    long long uSz[2] = {8, 64}, uSt[2] = {64, 1};   // Y: same shape/strides as U
    void* args[21] = {
        &dD, &dD, &zero, &dSz[0], &dSz[1], &dSt[0], &dSt[1],
        &dU, &dU, &zero, &uSz[0], &uSz[1], &uSt[0], &uSt[1],
        &dY, &dY, &zero, &uSz[0], &uSz[1], &uSt[0], &uSt[1],
    };

    CK(p_cuLaunchKernel(fn, 1, 1, 1, 32, 1, 1, 0, 0, args, nullptr));
    CK(p_cuCtxSynchronize());
    CK(p_cuMemcpyDtoH(out.data(), dY, M * N * sizeof(double)));

    double err = 0;
    for (int i = 0; i < M * N; ++i) err = fmax(err, fabs(out[i] - ref[i]));
    printf("mir-generated DMMA PTX gate: max|Y - D*U| = %.3e  %s\n",
           err, err < 1e-12 ? "PASS" : "FAIL");

    // Launch-latency sanity number only: a 1-warp toy kernel is NOT a benchmark.
    CUevent t0, t1; CK(p_cuEventCreate(&t0, 0)); CK(p_cuEventCreate(&t1, 0));
    const int iters = 10000;
    CK(p_cuEventRecord(t0, 0));
    for (int i = 0; i < iters; ++i)
        CK(p_cuLaunchKernel(fn, 1, 1, 1, 32, 1, 1, 0, 0, args, nullptr));
    CK(p_cuEventRecord(t1, 0)); CK(p_cuEventSynchronize(t1));
    float ms = 0; CK(p_cuEventElapsedTime(&ms, t0, t1));
    printf("  (%d launches, %.3f us/launch -- 1-warp toy, launch-latency bound)\n",
           iters, 1000.0f * ms / iters);
    return err < 1e-12 ? 0 : 1;
}
