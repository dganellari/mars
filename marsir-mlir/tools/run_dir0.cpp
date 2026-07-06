// GPU gate for the EMITTER-GENERATED full single direction (dir 0) of the
// Knaus HO apply (warp_dir0_sm90.ptx): B-sweep + 7 faces (each: 4 contracts +
// 3-term flux) + the +/- plane scatter into y, all warp-per-element and staged
// in shared memory (88 mma). Matches the dir-0 Knaus oracle. Metric is passed
// as three per-face component arrays g0all/g1all/g2all (8xnxn) to isolate the
// composition from the interleaved-G stride-3 layout (added later).
// ABI: u2(8x64), Btil,Dtil,Dm,W(8x8), g0all,g1all,g2all,y3(8x8x8) = 71 params
// (rank-2 memref = 7 scalars, rank-3 = 9). One warp; driver-API only.
// Build:  g++ -O2 run_dir0.cpp -o run_dir0 -ldl
// Run:    ./run_dir0 ../generated/warp_dir0_sm90.ptx

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

static void mm(const double* A, const double* B, double* C, int n)   // C=A@B
{
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double s = 0;
            for (int k = 0; k < n; ++k) s += A[i * n + k] * B[k * n + j];
            C[i * n + j] = s;
        }
}
static void mmT(const double* A, const double* B, double* C, int n)  // C=A@B^T
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

    const char* ptxPath = argc > 1 ? argv[1] : "../generated/warp_dir0_sm90.ptx";
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
    CUfunction fn; CK(p_cuModuleGetFunction(&fn, mod, "dir0"));

    const int n = 8, P = 7, nn = 64, nnn = 512;
    std::vector<double> u2(nn * 8 / 8 * 8);       // 8x64
    u2.assign(8 * nn, 0.0);
    std::vector<double> Btil(nn), Dtil(nn), Dm(nn), W(nn),
                        g0(nnn), g1(nnn), g2(nnn), y(nnn, 0.0), ref(nnn, 0.0), out(nnn);
    srand(42);
    for (auto* v : {&u2, &Btil, &Dtil, &Dm, &W, &g0, &g1, &g2})
        for (auto& x : *v) x = 2.0 * rand() / RAND_MAX - 1.0;

    // oracle (dir 0): interp_all = Btil @ u2 ; deriv_all = Dtil @ u2  (8x64)
    std::vector<double> interp_all(8 * nn), deriv_all(8 * nn);
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < nn; ++j) {
            double si = 0, sd = 0;
            for (int k = 0; k < n; ++k) { si += Btil[i * n + k] * u2[k * nn + j];
                                          sd += Dtil[i * n + k] * u2[k * nn + j]; }
            interp_all[i * nn + j] = si; deriv_all[i * nn + j] = sd;
        }
    for (int l = 0; l < P; ++l) {
        const double* interp = &interp_all[l * nn];   // 8x8 face
        const double* deriv = &deriv_all[l * nn];
        double dt2[64], dt1[64], flux[64], tmp[64], intf[64];
        mmT(interp, Dm.data(), dt2, n);   // dt2 = interp @ Dm^T
        mm(Dm.data(), interp, dt1, n);    // dt1 = Dm @ interp
        for (int i = 0; i < nn; ++i)
            flux[i] = g2[l * nn + i] * deriv[i] + g0[l * nn + i] * dt2[i]
                    + g1[l * nn + i] * dt1[i];
        mmT(flux, W.data(), tmp, n);      // tmp = flux @ W^T
        mm(W.data(), tmp, intf, n);       // intf = W @ tmp
        for (int i = 0; i < nn; ++i) { ref[l * nn + i] -= intf[i];
                                       ref[(l + 1) * nn + i] += intf[i]; }
    }

    // upload
    CUdeviceptr du2, dBt, dDt, dDm, dW, dg0, dg1, dg2, dy;
    CK(p_cuMemAlloc(&du2, 8 * nn * sizeof(double)));
    CK(p_cuMemAlloc(&dBt, nn * sizeof(double)));
    CK(p_cuMemAlloc(&dDt, nn * sizeof(double)));
    CK(p_cuMemAlloc(&dDm, nn * sizeof(double)));
    CK(p_cuMemAlloc(&dW, nn * sizeof(double)));
    CK(p_cuMemAlloc(&dg0, nnn * sizeof(double)));
    CK(p_cuMemAlloc(&dg1, nnn * sizeof(double)));
    CK(p_cuMemAlloc(&dg2, nnn * sizeof(double)));
    CK(p_cuMemAlloc(&dy, nnn * sizeof(double)));
    CK(p_cuMemcpyHtoD(du2, u2.data(), 8 * nn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dBt, Btil.data(), nn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dDt, Dtil.data(), nn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dDm, Dm.data(), nn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dW, W.data(), nn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dg0, g0.data(), nnn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dg1, g1.data(), nnn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dg2, g2.data(), nnn * sizeof(double)));
    CK(p_cuMemcpyHtoD(dy, y.data(), nnn * sizeof(double)));

    long long z = 0;
    long long u2Sz[2] = {8, 64}, u2St[2] = {64, 1};
    long long sz2[2] = {8, 8}, st2[2] = {8, 1};
    long long sz3[3] = {8, 8, 8}, st3[3] = {64, 8, 1};
    std::vector<void*> args;
    auto add2 = [&](CUdeviceptr* p, long long* sz, long long* st) {
        args.push_back(p); args.push_back(p); args.push_back(&z);
        args.push_back(&sz[0]); args.push_back(&sz[1]);
        args.push_back(&st[0]); args.push_back(&st[1]);
    };
    auto add3 = [&](CUdeviceptr* p) {
        args.push_back(p); args.push_back(p); args.push_back(&z);
        args.push_back(&sz3[0]); args.push_back(&sz3[1]); args.push_back(&sz3[2]);
        args.push_back(&st3[0]); args.push_back(&st3[1]); args.push_back(&st3[2]);
    };
    add2(&du2, u2Sz, u2St);
    add2(&dBt, sz2, st2); add2(&dDt, sz2, st2);
    add2(&dDm, sz2, st2); add2(&dW, sz2, st2);
    add3(&dg0); add3(&dg1); add3(&dg2); add3(&dy);

    CK(p_cuLaunchKernel(fn, 1, 1, 1, 32, 1, 1, 0, 0, args.data(), nullptr));
    CK(p_cuCtxSynchronize());
    CK(p_cuMemcpyDtoH(out.data(), dy, nnn * sizeof(double)));

    double err = 0;
    for (int i = 0; i < nnn; ++i) err = fmax(err, fabs(out[i] - ref[i]));
    printf("emitter dir-0 operator gate: max|y - Knaus_dir0(u)| = %.3e  %s\n",
           err, err < 1e-12 ? "PASS" : "FAIL");
    return err < 1e-12 ? 0 : 1;
}
