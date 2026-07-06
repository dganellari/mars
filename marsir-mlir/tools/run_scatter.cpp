// GPU gate for the +/- plane scatter (warp_scatter_sm90.ptx), direction 0, two
// faces with the loop-carried plane overlap. Verifies the census's hardest
// blocker distributes AND is numerically correct:
//   intf0 = A[:,0:4]@B[0:4,:] ; intf1 = A[:,4:8]@B[4:8,:]
//   y[0] = -intf0 ; y[1] = intf0 - intf1 ; y[2] = intf1 ; y[3..7] = 0
// One warp; Sf/Y are global memrefs (real design = shared memory). Driver-API
// only (dlopen libcuda, _v2 names). Y is rank-3 -> 9-scalar memref ABI.
// Build:  g++ -O2 run_scatter.cpp -o run_scatter -ldl
// Run:    ./run_scatter ../generated/warp_scatter_sm90.ptx

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

    const char* ptxPath = argc > 1 ? argv[1] : "../generated/warp_scatter_sm90.ptx";
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
    CUfunction fn; CK(p_cuModuleGetFunction(&fn, mod, "scatter"));

    const int n = 8, nn = 64, nnn = 512;
    std::vector<double> A(nn), B(nn), Sf(nn, 0.0), Y(nnn, 0.0), ref(nnn, 0.0), out(nnn);
    srand(42);
    for (auto& v : A) v = 2.0 * rand() / RAND_MAX - 1.0;
    for (auto& v : B) v = 2.0 * rand() / RAND_MAX - 1.0;

    // intf0 = A[:,0:4]@B[0:4,:] ; intf1 = A[:,4:8]@B[4:8,:]
    std::vector<double> intf0(nn, 0.0), intf1(nn, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double s0 = 0, s1 = 0;
            for (int k = 0; k < 4; ++k)     s0 += A[i * n + k]     * B[k * n + j];
            for (int k = 4; k < 8; ++k)     s1 += A[i * n + k]     * B[k * n + j];
            intf0[i * n + j] = s0; intf1[i * n + j] = s1;
        }
    for (int i = 0; i < nn; ++i) {
        ref[0 * nn + i] = -intf0[i];
        ref[1 * nn + i] = intf0[i] - intf1[i];
        ref[2 * nn + i] = intf1[i];
    }

    CUdeviceptr dA, dB, dSf, dY;
    CK(p_cuMemAlloc(&dA, nn * sizeof(double)));
    CK(p_cuMemAlloc(&dB, nn * sizeof(double)));
    CK(p_cuMemAlloc(&dSf, nn * sizeof(double)));
    CK(p_cuMemAlloc(&dY, nnn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dA, A.data(), nn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dB, B.data(), nn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dSf, Sf.data(), nn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dY, Y.data(), nnn * sizeof(double)));

    long long zero = 0;
    long long sz2[2] = {8, 8}, st2[2] = {8, 1};
    long long ySz[3] = {8, 8, 8}, ySt[3] = {64, 8, 1};
    // A(7) B(7) Sf(7) Y(9) = 30 scalars
    void* args[30] = {
        &dA, &dA, &zero, &sz2[0], &sz2[1], &st2[0], &st2[1],
        &dB, &dB, &zero, &sz2[0], &sz2[1], &st2[0], &st2[1],
        &dSf, &dSf, &zero, &sz2[0], &sz2[1], &st2[0], &st2[1],
        &dY, &dY, &zero, &ySz[0], &ySz[1], &ySz[2], &ySt[0], &ySt[1], &ySt[2],
    };

    CK(p_cuLaunchKernel(fn, 1, 1, 1, 32, 1, 1, 0, 0, args, nullptr));
    CK(p_cuCtxSynchronize());
    CK(p_cuMemcpyDtoH(out.data(), dY, nnn * sizeof(double)));

    double err = 0;
    for (int i = 0; i < nnn; ++i) err = fmax(err, fabs(out[i] - ref[i]));
    printf("warp +/- plane scatter gate: max|Y - ref| = %.3e  %s\n",
           err, err < 1e-12 ? "PASS" : "FAIL");
    return err < 1e-12 ? 0 : 1;
}
