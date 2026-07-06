// GPU gate for the single-direction contract chain (warp_one_dir_sm90.ptx):
// the real shape one direction of the HO operator takes in the warp-per-element
// target IR, staged across three warp regions (barriers between):
//   dt   = D @ If        (axis-0, standard arm, 2 k-slabs)
//   dts  = dt .* G       (flux, elementwise, fuses in-register)
//   tmp  = dts @ W^T     (axis-1, transposed-B arm, 2 slabs)
//   Out  = W @ tmp       (axis-0, standard arm, 2 slabs)
// One warp/block; scratch S1/S2 are plain global memrefs, gpu.barrier syncs the
// 32 lanes. This checks the cross-stage staging + both bridge arms against
// hardware (6 mma.sync total). Driver-API only (dlopen libcuda, _v2 names).
// Build:  g++ -O2 run_warp_one_dir.cpp -o run_warp_one_dir -ldl
// Run:    ./run_warp_one_dir ../generated/warp_one_dir_sm90.ptx

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

static void mm(const std::vector<double>& A, const std::vector<double>& B,
               std::vector<double>& C, int M, int K, int N)  // C = A@B
{
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j) {
            double s = 0;
            for (int k = 0; k < K; ++k) s += A[i * K + k] * B[k * N + j];
            C[i * N + j] = s;
        }
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

    const char* ptxPath = argc > 1 ? argv[1] : "../generated/warp_one_dir_sm90.ptx";
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
    CUfunction fn; CK(p_cuModuleGetFunction(&fn, mod, "one_dir"));

    const int n = 8, nn = 64;
    std::vector<double> If(nn), D(nn), G(nn), W(nn), zero8(nn, 0.0), out(nn);
    srand(42);
    for (auto& v : If) v = 2.0 * rand() / RAND_MAX - 1.0;
    for (auto& v : D)  v = 2.0 * rand() / RAND_MAX - 1.0;
    for (auto& v : G)  v = 2.0 * rand() / RAND_MAX - 1.0;
    for (auto& v : W)  v = 2.0 * rand() / RAND_MAX - 1.0;

    // reference
    std::vector<double> dt(nn), dts(nn), Wt(nn), tmp(nn), ref(nn);
    mm(D, If, dt, n, n, n);                       // dt = D @ If
    for (int i = 0; i < nn; ++i) dts[i] = dt[i] * G[i];
    for (int i = 0; i < n; ++i)                    // Wt = W^T
        for (int j = 0; j < n; ++j) Wt[i * n + j] = W[j * n + i];
    mm(dts, Wt, tmp, n, n, n);                     // tmp = dts @ W^T
    mm(W, tmp, ref, n, n, n);                      // Out = W @ tmp

    CUdeviceptr dIf, dD, dG, dW, dS1, dS2, dOut;
    for (auto* pp : {&dIf, &dD, &dG, &dW, &dS1, &dS2, &dOut})
        CK(p_cuMemAlloc(pp, nn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dIf, If.data(), nn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dD, D.data(), nn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dG, G.data(), nn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dW, W.data(), nn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dS1, zero8.data(), nn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dS2, zero8.data(), nn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dOut, zero8.data(), nn * sizeof(double)));

    long long zero = 0;
    long long sz[2] = {8, 8}, st[2] = {8, 1};
    auto pack = [&](void** a, CUdeviceptr* p) {
        a[0] = p; a[1] = p; a[2] = &zero;
        a[3] = &sz[0]; a[4] = &sz[1]; a[5] = &st[0]; a[6] = &st[1];
    };
    void* args[49];
    CUdeviceptr* ptrs[7] = {&dIf, &dD, &dG, &dW, &dS1, &dS2, &dOut};
    for (int i = 0; i < 7; ++i) pack(&args[i * 7], ptrs[i]);

    CK(p_cuLaunchKernel(fn, 1, 1, 1, 32, 1, 1, 0, 0, args, nullptr));
    CK(p_cuCtxSynchronize());
    CK(p_cuMemcpyDtoH(out.data(), dOut, nn * sizeof(double)));

    double err = 0;
    for (int i = 0; i < nn; ++i) err = fmax(err, fabs(out[i] - ref[i]));
    printf("warp one-direction chain gate: max|Out - W@((D@If).*G)@W^T| = %.3e  %s\n",
           err, err < 1e-12 ? "PASS" : "FAIL");
    return err < 1e-12 ? 0 : 1;
}
