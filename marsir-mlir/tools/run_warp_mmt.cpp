// GPU gate for the TRANSPOSED-B arm (warp_mmt_sm90.ptx): C := In @ W^T + C0,
// In/W/C all 8x8, two chained m8n8k4 slabs. W is stored [n,k] so its mma B
// fragment is read at [L/4, L%4]; this checks that transposed indexing against
// hardware, which IR inspection cannot. Same driver-API-only setup as
// run_warp_mma.cpp (dlopen libcuda, _v2 names).
// Build:  g++ -O2 run_warp_mmt.cpp -o run_warp_mmt -ldl
// Run:    ./run_warp_mmt ../generated/warp_mmt_sm90.ptx

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

    const char* ptxPath = argc > 1 ? argv[1] : "../generated/warp_mmt_sm90.ptx";
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
    CUfunction fn; CK(p_cuModuleGetFunction(&fn, mod, "warp_mmt"));

    const int M = 8, K = 8, N = 8;
    std::vector<double> hIn(M * K), hW(N * K), hC(M * N), ref(M * N), out(M * N);
    srand(42);
    for (auto& v : hIn) v = 2.0 * rand() / RAND_MAX - 1.0;
    for (auto& v : hW)  v = 2.0 * rand() / RAND_MAX - 1.0;
    for (auto& v : hC)  v = 2.0 * rand() / RAND_MAX - 1.0;
    ref = hC;
    // out[m,n] = C0[m,n] + sum_k In[m,k]*W[n,k]  (W stored [n,k])
    for (int m = 0; m < M; ++m)
        for (int n = 0; n < N; ++n) {
            double s = 0;
            for (int k = 0; k < K; ++k) s += hIn[m * K + k] * hW[n * K + k];
            ref[m * N + n] += s;
        }

    CUdeviceptr dIn, dW, dC;
    CK(p_cuMemAlloc(&dIn, M * K * sizeof(double)));
    CK(p_cuMemAlloc(&dW, N * K * sizeof(double)));
    CK(p_cuMemAlloc(&dC, M * N * sizeof(double)));
    CK(p_cuMemcpyHtoD(dIn, hIn.data(), M * K * sizeof(double)));
    CK(p_cuMemcpyHtoD(dW, hW.data(), N * K * sizeof(double)));
    CK(p_cuMemcpyHtoD(dC, hC.data(), M * N * sizeof(double)));

    long long zero = 0;
    long long iSz[2] = {8, 8}, iSt[2] = {8, 1};
    long long wSz[2] = {8, 8}, wSt[2] = {8, 1};
    long long cSz[2] = {8, 8}, cSt[2] = {8, 1};
    void* args[21] = {
        &dIn, &dIn, &zero, &iSz[0], &iSz[1], &iSt[0], &iSt[1],
        &dW, &dW, &zero, &wSz[0], &wSz[1], &wSt[0], &wSt[1],
        &dC, &dC, &zero, &cSz[0], &cSz[1], &cSt[0], &cSt[1],
    };

    CK(p_cuLaunchKernel(fn, 1, 1, 1, 32, 1, 1, 0, 0, args, nullptr));
    CK(p_cuCtxSynchronize());
    CK(p_cuMemcpyDtoH(out.data(), dC, M * N * sizeof(double)));

    double err = 0;
    for (int i = 0; i < M * N; ++i) err = fmax(err, fabs(out[i] - ref[i]));
    printf("warp->mmt (transposed-B) gate: max|C - (In@W^T + C0)| = %.3e  %s\n",
           err, err < 1e-12 ? "PASS" : "FAIL");
    return err < 1e-12 ? 0 : 1;
}
