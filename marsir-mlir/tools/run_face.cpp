// GPU gate for the EMITTER-GENERATED single face of the Knaus apply
// (warp_face_sm90.ptx): the operator's repeating unit -- 4 contracts + the
// 3-term CVFEM flux, staged through shared memory (8 mma). Matches the Knaus
// per-face oracle:
//   dt2  = interp @ Dm^T ; dt1 = Dm @ interp
//   flux = g2*deriv + g0*dt2 + g1*dt1
//   tmp  = flux @ W^T ; intf = W @ tmp
// Inputs interp,deriv,Dm,W,g0,g1,g2 + output intf, all 8x8 -> 8 memrefs = 56
// params; dt1g/dt2g/flux/tmp are internal .shared. One warp; driver-API only.
// Build:  g++ -O2 run_face.cpp -o run_face -ldl
// Run:    ./run_face ../generated/warp_face_sm90.ptx

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

static void matmul(const double* A, const double* B, double* C, int n)  // C=A@B
{
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double s = 0;
            for (int k = 0; k < n; ++k) s += A[i * n + k] * B[k * n + j];
            C[i * n + j] = s;
        }
}
static void matmulT(const double* A, const double* B, double* C, int n)  // C=A@B^T
{
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double s = 0;
            for (int k = 0; k < n; ++k) s += A[i * n + k] * B[j * n + k];
            C[i * n + j] = s;
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

    const char* ptxPath = argc > 1 ? argv[1] : "../generated/warp_face_sm90.ptx";
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
    CUfunction fn; CK(p_cuModuleGetFunction(&fn, mod, "face"));

    const int n = 8, nn = 64;
    std::vector<double> interp(nn), deriv(nn), Dm(nn), W(nn),
                        g0(nn), g1(nn), g2(nn), out(nn), zero(nn, 0.0);
    srand(42);
    for (auto* v : {&interp, &deriv, &Dm, &W, &g0, &g1, &g2})
        for (auto& x : *v) x = 2.0 * rand() / RAND_MAX - 1.0;

    // oracle
    std::vector<double> dt2(nn), dt1(nn), flux(nn), tmp(nn), ref(nn);
    matmulT(interp.data(), Dm.data(), dt2.data(), n);   // dt2 = interp @ Dm^T
    matmul(Dm.data(), interp.data(), dt1.data(), n);    // dt1 = Dm @ interp
    for (int i = 0; i < nn; ++i)
        flux[i] = g2[i] * deriv[i] + g0[i] * dt2[i] + g1[i] * dt1[i];
    matmulT(flux.data(), W.data(), tmp.data(), n);      // tmp = flux @ W^T
    matmul(W.data(), tmp.data(), ref.data(), n);        // intf = W @ tmp

    CUdeviceptr d[8];
    std::vector<double>* host[7] = {&interp, &deriv, &Dm, &W, &g0, &g1, &g2};
    for (int i = 0; i < 8; ++i) CK(p_cuMemAlloc(&d[i], nn * sizeof(double)));
    for (int i = 0; i < 7; ++i) CK(p_cuMemcpyHtoD(d[i], host[i]->data(), nn * sizeof(double)));
    CK(p_cuMemcpyHtoD(d[7], zero.data(), nn * sizeof(double)));  // intf

    long long zoff = 0, sz[2] = {8, 8}, st[2] = {8, 1};
    void* args[56];
    for (int i = 0; i < 8; ++i) {
        void** a = &args[i * 7];
        a[0] = &d[i]; a[1] = &d[i]; a[2] = &zoff;
        a[3] = &sz[0]; a[4] = &sz[1]; a[5] = &st[0]; a[6] = &st[1];
    }

    CK(p_cuLaunchKernel(fn, 1, 1, 1, 32, 1, 1, 0, 0, args, nullptr));
    CK(p_cuCtxSynchronize());
    CK(p_cuMemcpyDtoH(out.data(), d[7], nn * sizeof(double)));

    double err = 0;
    for (int i = 0; i < nn; ++i) err = fmax(err, fabs(out[i] - ref[i]));
    printf("emitter face-chain gate: max|intf - W@((g2*deriv+g0*dt2+g1*dt1)@W^T)| = %.3e  %s\n",
           err, err < 1e-12 ? "PASS" : "FAIL");
    return err < 1e-12 ? 0 : 1;
}
