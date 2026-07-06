// GPU gate for EMITTER-GENERATED warp-form code (warp_emit_selftest_sm90.ptx):
// the first hardware check of output from marsir-compiler's mlir_warp.py (not a
// hand-written .mlir). Kernel @two stages two chained contracts through shared
// memory:  S1 = M @ X  (standard arm) ; Out = S1 @ M^T  (transposed-B arm).
// S1 is workgroup/.shared (internal, not a param) -> 3 memref args = 21 params.
// One warp; driver-API only (dlopen libcuda, _v2 names).
// Build:  g++ -O2 run_warp_emit.cpp -o run_warp_emit -ldl
// Run:    ./run_warp_emit ../generated/warp_emit_selftest_sm90.ptx

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

    const char* ptxPath = argc > 1 ? argv[1] : "../generated/warp_emit_selftest_sm90.ptx";
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
    CUfunction fn; CK(p_cuModuleGetFunction(&fn, mod, "two"));

    const int n = 8, nn = 64;
    std::vector<double> M(nn), X(nn), out(nn), zero8(nn, 0.0);
    srand(42);
    for (auto& v : M) v = 2.0 * rand() / RAND_MAX - 1.0;
    for (auto& v : X) v = 2.0 * rand() / RAND_MAX - 1.0;

    // reference: S1 = M @ X ; Out = S1 @ M^T
    std::vector<double> S1(nn, 0.0), ref(nn, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double s = 0;
            for (int k = 0; k < n; ++k) s += M[i * n + k] * X[k * n + j];
            S1[i * n + j] = s;
        }
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double s = 0;
            for (int k = 0; k < n; ++k) s += S1[i * n + k] * M[j * n + k];  // M^T
            ref[i * n + j] = s;
        }

    CUdeviceptr dM, dX, dOut;
    CK(p_cuMemAlloc(&dM, nn * sizeof(double)));
    CK(p_cuMemAlloc(&dX, nn * sizeof(double)));
    CK(p_cuMemAlloc(&dOut, nn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dM, M.data(), nn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dX, X.data(), nn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dOut, zero8.data(), nn * sizeof(double)));

    long long zero = 0;
    long long sz[2] = {8, 8}, st[2] = {8, 1};
    auto pack = [&](void** a, CUdeviceptr* p) {
        a[0] = p; a[1] = p; a[2] = &zero;
        a[3] = &sz[0]; a[4] = &sz[1]; a[5] = &st[0]; a[6] = &st[1];
    };
    void* args[21];
    CUdeviceptr* ptrs[3] = {&dM, &dX, &dOut};
    for (int i = 0; i < 3; ++i) pack(&args[i * 7], ptrs[i]);

    CK(p_cuLaunchKernel(fn, 1, 1, 1, 32, 1, 1, 0, 0, args, nullptr));
    CK(p_cuCtxSynchronize());
    CK(p_cuMemcpyDtoH(out.data(), dOut, nn * sizeof(double)));

    double err = 0;
    for (int i = 0; i < nn; ++i) err = fmax(err, fabs(out[i] - ref[i]));
    printf("emitter-generated warp gate: max|Out - (M@X)@M^T| = %.3e  %s\n",
           err, err < 1e-12 ? "PASS" : "FAIL");
    return err < 1e-12 ? 0 : 1;
}
