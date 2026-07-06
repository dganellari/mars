// GPU gate for the emitter's rectangular B-sweep matmul (warp_bsweep_sm90.ptx):
// interp = Btil(8x8) @ U(8x64) -> 8x64, column-tiled to 8 output blocks x 2
// k-slabs = 16 mma. This is the sum-factorization B-sweep shape (Btil P->8
// padded @ U reshaped n x n^2); checks the column-tile indexing against
// hardware. One warp; driver-API only (dlopen libcuda, _v2 names).
// Build:  g++ -O2 run_bsweep.cpp -o run_bsweep -ldl
// Run:    ./run_bsweep ../generated/warp_bsweep_sm90.ptx

#include <dlfcn.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

typedef int CUresult;
typedef int CUdevice;
typedef unsigned long long CUdeviceptr;
typedef struct CUctx_st* CUcontext;
typedef struct CUmod_st* CUmodule;
typedef struct CUfunc_st* CUfunction;
typedef struct CUstream_st* CUstream;

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
    *(void**)&p_cuGetErrorString    = must_sym(lib, "cuGetErrorString");

    const char* ptxPath = argc > 1 ? argv[1] : "../generated/warp_bsweep_sm90.ptx";
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
    CUfunction fn; CK(p_cuModuleGetFunction(&fn, mod, "bsweep"));

    const int m = 8, k = 8, N = 64;
    std::vector<double> Bt(m * k), U(k * N), out(m * N), zero(m * N, 0.0), ref(m * N, 0.0);
    srand(42);
    for (auto& v : Bt) v = 2.0 * rand() / RAND_MAX - 1.0;
    for (auto& v : U)  v = 2.0 * rand() / RAND_MAX - 1.0;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < N; ++j) {
            double s = 0;
            for (int p = 0; p < k; ++p) s += Bt[i * k + p] * U[p * N + j];
            ref[i * N + j] = s;
        }

    CUdeviceptr dBt, dU, dI;
    CK(p_cuMemAlloc(&dBt, m * k * sizeof(double)));
    CK(p_cuMemAlloc(&dU, k * N * sizeof(double)));
    CK(p_cuMemAlloc(&dI, m * N * sizeof(double)));
    CK(p_cuMemcpyHtoD(dBt, Bt.data(), m * k * sizeof(double)));
    CK(p_cuMemcpyHtoD(dU, U.data(), k * N * sizeof(double)));
    CK(p_cuMemcpyHtoD(dI, zero.data(), m * N * sizeof(double)));

    long long zero_off = 0;
    long long btSz[2] = {8, 8},  btSt[2] = {8, 1};
    long long uSz[2] = {8, 64},  uSt[2] = {64, 1};
    void* args[21] = {
        &dBt, &dBt, &zero_off, &btSz[0], &btSz[1], &btSt[0], &btSt[1],
        &dU, &dU, &zero_off, &uSz[0], &uSz[1], &uSt[0], &uSt[1],
        &dI, &dI, &zero_off, &uSz[0], &uSz[1], &uSt[0], &uSt[1],
    };

    CK(p_cuLaunchKernel(fn, 1, 1, 1, 32, 1, 1, 0, 0, args, nullptr));
    CK(p_cuCtxSynchronize());
    CK(p_cuMemcpyDtoH(out.data(), dI, m * N * sizeof(double)));

    double err = 0;
    for (int i = 0; i < m * N; ++i) err = fmax(err, fabs(out[i] - ref[i]));
    printf("emitter B-sweep matmul gate: max|Interp - Btil@U| = %.3e  %s\n",
           err, err < 1e-12 ? "PASS" : "FAIL");
    return err < 1e-12 ? 0 : 1;
}
