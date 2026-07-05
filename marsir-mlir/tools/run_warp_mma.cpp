// GPU gate for the warp-distribution -> tensor-core bridge (warp_mma_sm90.ptx):
// C := A(8x4) * B(4x8) + C, one warp, one mma.sync.m8n8k4. The kernel was
// authored as a vector.contract INSIDE warp_execute_on_lane_0 and distributed
// by --mir-warp-distribute (WarpOpContractToMma) -- this checks the per-lane
// fragment indexing (A[L/4,L%4], B[L%4,L/4], C[L/4,2(L%4)]) against the
// hardware, which no amount of IR inspection can.
//
// Same driver-API-only setup as run_ptx_gate.cpp (dlopen libcuda, _v2 names).
// Build:  g++ -O2 run_warp_mma.cpp -o run_warp_mma -ldl
// Run:    ./run_warp_mma ../generated/warp_mma_sm90.ptx

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

    const char* ptxPath = argc > 1 ? argv[1] : "../generated/warp_mma_sm90.ptx";
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
    CUfunction fn; CK(p_cuModuleGetFunction(&fn, mod, "warp_mma"));

    const int M = 8, K = 4, N = 8;
    std::vector<double> hA(M * K), hB(K * N), hC(M * N), ref(M * N), out(M * N);
    srand(42);
    for (auto& v : hA) v = 2.0 * rand() / RAND_MAX - 1.0;
    for (auto& v : hB) v = 2.0 * rand() / RAND_MAX - 1.0;
    for (auto& v : hC) v = 2.0 * rand() / RAND_MAX - 1.0;
    ref = hC;
    for (int i = 0; i < M; ++i)
        for (int p = 0; p < K; ++p) {
            const double a = hA[i * K + p];
            for (int j = 0; j < N; ++j) ref[i * N + j] += a * hB[p * N + j];
        }

    CUdeviceptr dA, dB, dC;
    CK(p_cuMemAlloc(&dA, M * K * sizeof(double)));
    CK(p_cuMemAlloc(&dB, K * N * sizeof(double)));
    CK(p_cuMemAlloc(&dC, M * N * sizeof(double)));
    CK(p_cuMemcpyHtoD(dA, hA.data(), M * K * sizeof(double)));
    CK(p_cuMemcpyHtoD(dB, hB.data(), K * N * sizeof(double)));
    CK(p_cuMemcpyHtoD(dC, hC.data(), M * N * sizeof(double)));

    long long zero = 0;
    long long aSz[2] = {8, 4}, aSt[2] = {4, 1};
    long long bSz[2] = {4, 8}, bSt[2] = {8, 1};
    long long cSz[2] = {8, 8}, cSt[2] = {8, 1};
    void* args[21] = {
        &dA, &dA, &zero, &aSz[0], &aSz[1], &aSt[0], &aSt[1],
        &dB, &dB, &zero, &bSz[0], &bSz[1], &bSt[0], &bSt[1],
        &dC, &dC, &zero, &cSz[0], &cSz[1], &cSt[0], &cSt[1],
    };

    CK(p_cuLaunchKernel(fn, 1, 1, 1, 32, 1, 1, 0, 0, args, nullptr));
    CK(p_cuCtxSynchronize());
    CK(p_cuMemcpyDtoH(out.data(), dC, M * N * sizeof(double)));

    double err = 0;
    for (int i = 0; i < M * N; ++i) err = fmax(err, fabs(out[i] - ref[i]));
    printf("warp->mma bridge gate: max|C - (A*B + C0)| = %.3e  %s\n",
           err, err < 1e-12 ? "PASS" : "FAIL");
    return err < 1e-12 ? 0 : 1;
}
