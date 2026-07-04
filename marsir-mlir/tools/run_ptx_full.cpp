// GPU gate + timing for the FUSED FULL-OPERATOR kernel (full_apply_pN_sm90.ptx):
// the whole Knaus Alg-2 apply, one thread per element, 128 threads/block, ONE
// kernel launch. Self-contained (dlopen libcuda.so.1).
//
// Launch-operand order (from gpu.launch_func):
//   [i64 1][i64 0][i64 tpb][U desc 11][G desc 15][Y desc 11][f64 0.0]
//   [i64 n][i64 n*n][i64 P][Btil desc 7][Dtil desc 7][D desc 7][W desc 7]
// U/Y: memref<?xnxnxnxf64>; G: memref<?x3xPxnxnx3xf64>; ops: PxN / nxn static.
//
// Correctness: spot-check vs an in-process CPU port of applyHoCvfemElement.
// NOTE: thread-per-element SCALAR kernel (no tensor cores yet) with large
// per-thread locals at p=7 -- a fusion baseline, not the tuned endpoint.
//
// Build:  g++ -O2 run_ptx_full.cpp -o run_ptx_full -ldl
// Run:    ./run_ptx_full ../generated/full_apply_p3_sm90.ptx 1048576 3

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

// CPU port of the Knaus Alg-2 apply (applyHoCvfemElement) for one element.
static void oracle(int p, const double* u, const double* Bt, const double* Dt,
                   const double* Dm, const double* W, const double* G,
                   double* y)
{
    const int n = p + 1, nn = n * n;
    auto idx = [&](int d, int nrm, int t1, int t2) {
        if (d == 0) return nrm * nn + t1 * n + t2;
        if (d == 1) return t1 * nn + nrm * n + t2;
        return t1 * nn + t2 * n + nrm;
    };
    for (int i = 0; i < n * nn; ++i) y[i] = 0.0;
    std::vector<double> interp(nn), deriv(nn), flux(nn), tmp(nn), intf(nn);
    for (int d = 0; d < 3; ++d)
        for (int l = 0; l < p; ++l) {
            for (int s = 0; s < n; ++s)
                for (int r = 0; r < n; ++r) {
                    double bi = 0, di = 0;
                    for (int q = 0; q < n; ++q) {
                        bi += Bt[l * n + q] * u[idx(d, q, s, r)];
                        di += Dt[l * n + q] * u[idx(d, q, s, r)];
                    }
                    interp[s * n + r] = bi; deriv[s * n + r] = di;
                }
            for (int s = 0; s < n; ++s)
                for (int r = 0; r < n; ++r) {
                    double dt2 = 0, dt1 = 0;
                    for (int q = 0; q < n; ++q) dt2 += Dm[r * n + q] * interp[s * n + q];
                    for (int q = 0; q < n; ++q) dt1 += Dm[s * n + q] * interp[q * n + r];
                    const double* g = G + (((d * p + l) * n + s) * n + r) * 3;
                    flux[s * n + r] = g[2] * deriv[s * n + r] + g[0] * dt2 + g[1] * dt1;
                }
            for (int s = 0; s < n; ++s)
                for (int r = 0; r < n; ++r) { double v = 0;
                    for (int q = 0; q < n; ++q) v += W[r * n + q] * flux[s * n + q];
                    tmp[s * n + r] = v; }
            for (int s = 0; s < n; ++s)
                for (int r = 0; r < n; ++r) { double v = 0;
                    for (int q = 0; q < n; ++q) v += W[s * n + q] * tmp[q * n + r];
                    intf[s * n + r] = v; }
            for (int s = 0; s < n; ++s)
                for (int r = 0; r < n; ++r) {
                    y[idx(d, l, s, r)]     -= intf[s * n + r];
                    y[idx(d, l + 1, s, r)] += intf[s * n + r];
                }
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
    *(void**)&p_cuEventCreate       = must_sym(lib, "cuEventCreate");
    *(void**)&p_cuEventRecord       = must_sym(lib, "cuEventRecord");
    *(void**)&p_cuEventSynchronize  = must_sym(lib, "cuEventSynchronize");
    *(void**)&p_cuEventElapsedTime  = must_sym(lib, "cuEventElapsedTime");
    *(void**)&p_cuGetErrorString    = must_sym(lib, "cuGetErrorString");
    *(void**)&p_cuMemsetD8          = must_sym(lib, "cuMemsetD8_v2");

    const char* ptxPath = argc > 1 ? argv[1] : "../generated/full_apply_p3_sm90.ptx";
    const long long E = argc > 2 ? atoll(argv[2]) : (1LL << 20);
    const int p = argc > 3 ? atoi(argv[3]) : 3;
    const int tpb = argc > 4 ? atoi(argv[4]) : 128;
    if (E % tpb) { fprintf(stderr, "E must be divisible by tpb\n"); return 1; }
    const int n = p + 1, nn = n * n, n3 = nn * n;
    const long long gElem = 3LL * p * nn * 3;

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
    CUfunction fn; CK(p_cuModuleGetFunction(&fn, mod, "laplacian_apply_batched_kernel"));

    std::vector<double> hU((size_t)E * n3), hBt((size_t)p * n), hDt((size_t)p * n),
        hDm((size_t)nn), hW((size_t)nn), hG((size_t)E * gElem);
    srand(42);
    auto rnd = [] { return 2.0 * rand() / RAND_MAX - 1.0; };
    for (auto& x : hU) x = rnd();
    for (auto& x : hBt) x = rnd();
    for (auto& x : hDt) x = rnd();
    for (auto& x : hDm) x = rnd();
    for (auto& x : hW) x = rnd();
    for (auto& x : hG) x = rnd();

    CUdeviceptr dU, dY, dG, dBt, dDt, dDm, dW;
    CK(p_cuMemAlloc(&dU, hU.size() * 8));
    CK(p_cuMemAlloc(&dY, hU.size() * 8));
    CK(p_cuMemAlloc(&dG, hG.size() * 8));
    CK(p_cuMemAlloc(&dBt, hBt.size() * 8));
    CK(p_cuMemAlloc(&dDt, hDt.size() * 8));
    CK(p_cuMemAlloc(&dDm, hDm.size() * 8));
    CK(p_cuMemAlloc(&dW, hW.size() * 8));
    CK(p_cuMemcpyHtoD(dU, hU.data(), hU.size() * 8));
    CK(p_cuMemsetD8(dY, 0, hU.size() * 8));
    CK(p_cuMemcpyHtoD(dG, hG.data(), hG.size() * 8));
    CK(p_cuMemcpyHtoD(dBt, hBt.data(), hBt.size() * 8));
    CK(p_cuMemcpyHtoD(dDt, hDt.data(), hDt.size() * 8));
    CK(p_cuMemcpyHtoD(dDm, hDm.data(), hDm.size() * 8));
    CK(p_cuMemcpyHtoD(dW, hW.data(), hW.size() * 8));

    long long i1 = 1, i0 = 0, itpb = tpb, in = n, inn = nn, iP = p, zero = 0;
    double fzero = 0.0;
    long long uSz[4] = {E, n, n, n}, uSt[4] = {(long long)n3, (long long)nn, n, 1};
    long long gSz[6] = {E, 3, p, n, n, 3};
    long long gSt[6] = {gElem, (long long)p * nn * 3, (long long)nn * 3,
                        (long long)n * 3, 3, 1};
    long long oSz[2] = {p, n}, oSt[2] = {n, 1};       // Btil/Dtil (P x n)
    long long sSz[2] = {n, n}, sSt[2] = {n, 1};       // D/W (n x n)
    void* args[] = {
        &i1, &i0, &itpb,
        &dU, &dU, &zero, &uSz[0], &uSz[1], &uSz[2], &uSz[3],
                         &uSt[0], &uSt[1], &uSt[2], &uSt[3],
        &dG, &dG, &zero, &gSz[0], &gSz[1], &gSz[2], &gSz[3], &gSz[4], &gSz[5],
                         &gSt[0], &gSt[1], &gSt[2], &gSt[3], &gSt[4], &gSt[5],
        &dY, &dY, &zero, &uSz[0], &uSz[1], &uSz[2], &uSz[3],
                         &uSt[0], &uSt[1], &uSt[2], &uSt[3],
        &fzero, &in, &inn, &iP,
        &dBt, &dBt, &zero, &oSz[0], &oSz[1], &oSt[0], &oSt[1],
        &dDt, &dDt, &zero, &oSz[0], &oSz[1], &oSt[0], &oSt[1],
        &dDm, &dDm, &zero, &sSz[0], &sSz[1], &sSt[0], &sSt[1],
        &dW,  &dW,  &zero, &sSz[0], &sSz[1], &sSt[0], &sSt[1],
    };

    const unsigned grid = (unsigned)(E / tpb);
    CK(p_cuLaunchKernel(fn, grid, 1, 1, (unsigned)tpb, 1, 1, 0, 0, args, nullptr));
    CK(p_cuCtxSynchronize());

    double err = 0;
    std::vector<double> out(n3), ref(n3);
    for (int probe = 0; probe < 8; ++probe) {
        long long e = (long long)((double)rand() / RAND_MAX * (E - 1));
        CK(p_cuMemcpyDtoH(out.data(), dY + (size_t)e * n3 * 8, n3 * 8));
        oracle(p, &hU[(size_t)e * n3], hBt.data(), hDt.data(), hDm.data(),
               hW.data(), &hG[(size_t)e * gElem], ref.data());
        for (int i = 0; i < n3; ++i) err = fmax(err, fabs(out[i] - ref[i]));
    }
    printf("FULL operator gate (p=%d, E=%lld, fused 1-kernel): spot max|err| = %.3e  %s\n",
           p, E, err, err < 1e-11 ? "PASS" : "FAIL");

    CUevent t0, t1; CK(p_cuEventCreate(&t0, 0)); CK(p_cuEventCreate(&t1, 0));
    const int iters = 20;
    for (int w = 0; w < 3; ++w)
        CK(p_cuLaunchKernel(fn, grid, 1, 1, (unsigned)tpb, 1, 1, 0, 0, args, nullptr));
    CK(p_cuCtxSynchronize());
    CK(p_cuEventRecord(t0, 0));
    for (int i = 0; i < iters; ++i)
        CK(p_cuLaunchKernel(fn, grid, 1, 1, (unsigned)tpb, 1, 1, 0, 0, args, nullptr));
    CK(p_cuEventRecord(t1, 0)); CK(p_cuEventSynchronize(t1));
    float ms = 0; CK(p_cuEventElapsedTime(&ms, t0, t1));
    const double msApply = ms / iters;
    // Useful traffic: u + y-write + per-point metric G (dominant, like the hand
    // kernel's d_G stream).
    const double bytesElem = (2.0 * n3 + gElem) * 8.0;
    printf("  %d applies: %.3f ms/apply | %.1f ns/elem | ~%.0f GB/s useful (%.1f KiB/elem)\n",
           iters, msApply, msApply * 1e6 / (double)E,
           (double)E * bytesElem / (msApply * 1e-3) / 1e9, bytesElem / 1024.0);
    printf("  NOTE: scalar thread-per-element fusion baseline (no tensor cores yet).\n");
    return err < 1e-11 ? 0 : 1;
}
