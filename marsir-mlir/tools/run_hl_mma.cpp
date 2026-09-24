// GPU gate + timing for the kernel MARSIR generates from HIGH-LEVEL mir by
// PASSES ONLY (marsir-mlir/make_hl_mma.sh). One WARP per element, grid = E,
// blockDim = 32 -- the register-resident form, so this is the counterpart of
// run_batched.cpp but for the pass-generated kernel rather than the Python-
// emitted one.
//
// Checks 8 random elements against the same Knaus oracle run_ptx_full.cpp uses,
// then times the kernel. Element data is a deterministic hash of (element, slot)
// so any element can be regenerated on the host without a full mirror.
//
// Build:  g++ -O2 run_hl_mma.cpp -o run_hl_mma -ldl
// Run:    ./run_hl_mma ../generated/hl_full_sm90.ptx [E] [p]
#include <dlfcn.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "knaus_oracle.h"

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

    const char* ptxPath = argc > 1 ? argv[1] : "../generated/hl_full_sm90.ptx";
    const long long E = argc > 2 ? atoll(argv[2]) : (1LL << 20);
    const int p = argc > 3 ? atoi(argv[3]) : 7;
    const unsigned WARP = 32;   // one warp per element; grid = E
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
    CUfunction fn; CK(p_cuModuleGetFunction(&fn, mod, "laplacian_apply"));

    std::vector<double> hBt((size_t)p * n), hDt((size_t)p * n),
        hDm((size_t)nn), hW((size_t)nn);
    srand(42);
    auto rnd = [] { return 2.0 * rand() / RAND_MAX - 1.0; };
    for (auto& x : hBt) x = rnd();
    for (auto& x : hDt) x = rnd();
    for (auto& x : hDm) x = rnd();
    for (auto& x : hW) x = rnd();

    std::vector<double> elemU((size_t)n3), elemG((size_t)gElem);
    auto fillElem = [&](long long e) {
        for (size_t i = 0; i < (size_t)n3; ++i)    elemU[i] = hashVal(e, i, 0u);
        for (size_t i = 0; i < (size_t)gElem; ++i) elemG[i] = hashVal(e, i, 0x9e37u);
    };

    CUdeviceptr dU, dY, dG, dBt, dDt, dDm, dW;
    const size_t uBytes = (size_t)E * n3 * 8, gBytes = (size_t)E * gElem * 8;
    CK(p_cuMemAlloc(&dU, uBytes));
    CK(p_cuMemAlloc(&dY, uBytes));
    CK(p_cuMemAlloc(&dG, gBytes));
    CK(p_cuMemAlloc(&dBt, hBt.size() * 8));
    CK(p_cuMemAlloc(&dDt, hDt.size() * 8));
    CK(p_cuMemAlloc(&dDm, hDm.size() * 8));
    CK(p_cuMemAlloc(&dW, hW.size() * 8));
    {   // Chunked staging: bounded host memory, and O(E/chunk) uploads rather than
        // the O(E) driver calls a per-element loop would cost.
        const size_t perElem = (size_t)(n3 + gElem) * 8;
        size_t chunk = (64u << 20) / (perElem ? perElem : 1);
        if (chunk == 0) chunk = 1;
        if (chunk > (size_t)E) chunk = (size_t)E;
        std::vector<double> sU(chunk * n3), sG(chunk * gElem);
        for (long long e0 = 0; e0 < E; e0 += (long long)chunk) {
            const size_t cnt = (size_t)((E - e0) < (long long)chunk ? (E - e0)
                                                                    : (long long)chunk);
            for (size_t c = 0; c < cnt; ++c) {
                fillElem(e0 + (long long)c);
                std::copy(elemU.begin(), elemU.end(), sU.begin() + c * n3);
                std::copy(elemG.begin(), elemG.end(), sG.begin() + c * gElem);
            }
            CK(p_cuMemcpyHtoD(dU + (size_t)e0 * n3 * 8, sU.data(), cnt * n3 * 8));
            CK(p_cuMemcpyHtoD(dG + (size_t)e0 * gElem * 8, sG.data(), cnt * gElem * 8));
        }
    }
    CK(p_cuMemsetD8(dY, 0, uBytes));
    CK(p_cuMemcpyHtoD(dBt, hBt.data(), hBt.size() * 8));
    CK(p_cuMemcpyHtoD(dDt, hDt.data(), hDt.size() * 8));
    CK(p_cuMemcpyHtoD(dDm, hDm.data(), hDm.size() * 8));
    CK(p_cuMemcpyHtoD(dW, hW.data(), hW.size() * 8));

    long long zero = 0;
    long long uSz[4] = {E, n, n, n}, uSt[4] = {(long long)n3, (long long)nn, n, 1};
    long long gSz[6] = {E, 3, p, 3, n, n};   // [dir][face][component][row][col]
    long long gSt[6] = {gElem, (long long)p * 3 * nn, 3LL * nn, (long long)nn, n, 1};
    long long oSz[2] = {p, n}, oSt[2] = {n, 1};       // Btil/Dtil (P x n)
    long long sSz[2] = {n, n}, sSt[2] = {n, 1};       // D/W (n x n)
    void* args[] = {
        // laplacian_apply(U, Btil, Dtil, Dm, W, G, Y) -- MLIR memref descriptors:
        // allocated ptr, aligned ptr, offset, sizes[rank], strides[rank].
        &dU,  &dU,  &zero, &uSz[0], &uSz[1], &uSz[2], &uSz[3],
                           &uSt[0], &uSt[1], &uSt[2], &uSt[3],
        &dBt, &dBt, &zero, &oSz[0], &oSz[1], &oSt[0], &oSt[1],
        &dDt, &dDt, &zero, &oSz[0], &oSz[1], &oSt[0], &oSt[1],
        // ORDER: the operator is emitted as (u, Btil, Dtil, W, D, G) -- W BEFORE D.
        // Both are memref<8x8>, so swapping them is type-invisible and silently
        // yields a plausible-looking wrong answer.
        &dW,  &dW,  &zero, &sSz[0], &sSz[1], &sSt[0], &sSt[1],
        &dDm, &dDm, &zero, &sSz[0], &sSz[1], &sSt[0], &sSt[1],
        &dG,  &dG,  &zero, &gSz[0], &gSz[1], &gSz[2], &gSz[3], &gSz[4], &gSz[5],
                           &gSt[0], &gSt[1], &gSt[2], &gSt[3], &gSt[4], &gSt[5],
        &dY,  &dY,  &zero, &uSz[0], &uSz[1], &uSz[2], &uSz[3],
                           &uSt[0], &uSt[1], &uSt[2], &uSt[3],
    };

    const unsigned grid = (unsigned)E;
    CK(p_cuLaunchKernel(fn, grid, 1, 1, WARP, 1, 1, 0, 0, args, nullptr));
    CK(p_cuCtxSynchronize());

    double err = 0;
    std::vector<double> out(n3), ref(n3);
    // Probes: element 0, the last element, then random ones. Printing |out| next
    // to |ref| per probe separates the failure modes: |out| = 0 everywhere means
    // Y was never written; correct at e = 0 but zero elsewhere means every warp
    // wrote element 0; nonzero-but-wrong means the math itself is off.
    const int nProbe = 10;
    for (int probe = 0; probe < nProbe; ++probe) {
        long long e = probe == 0 ? 0 : probe == 1 ? E - 1
                    : (long long)((double)rand() / RAND_MAX * (E - 1));
        CK(p_cuMemcpyDtoH(out.data(), dY + (size_t)e * n3 * 8, n3 * 8));
        fillElem(e);
        oracle(p, elemU.data(), hBt.data(), hDt.data(), hDm.data(),
               hW.data(), elemG.data(), ref.data());
        double eo = 0, mo = 0, mr = 0;
        for (int i = 0; i < n3; ++i) {
            eo = fmax(eo, fabs(out[i] - ref[i]));
            mo = fmax(mo, fabs(out[i]));
            mr = fmax(mr, fabs(ref[i]));
        }
        printf("  probe e=%-8lld max|out|=%.3e  max|ref|=%.3e  max|out-ref|=%.3e  y[0..2]=% .4e % .4e % .4e  ref=% .4e % .4e % .4e\n",
               e, mo, mr, eo, out[0], out[1], out[2], ref[0], ref[1], ref[2]);
        err = fmax(err, eo);
    }
    printf("HIGH-LEVEL mir -> mma gate (p=%d, E=%lld, warp/elem): spot max|err| = %.3e  %s\n",
           p, E, err, err < 1e-11 ? "PASS" : "FAIL");

    CUevent t0, t1; CK(p_cuEventCreate(&t0, 0)); CK(p_cuEventCreate(&t1, 0));
    const int iters = 20;
    for (int w = 0; w < 3; ++w)
        CK(p_cuLaunchKernel(fn, grid, 1, 1, WARP, 1, 1, 0, 0, args, nullptr));
    CK(p_cuCtxSynchronize());
    CK(p_cuEventRecord(t0, 0));
    for (int i = 0; i < iters; ++i)
        CK(p_cuLaunchKernel(fn, grid, 1, 1, WARP, 1, 1, 0, 0, args, nullptr));
    CK(p_cuEventRecord(t1, 0)); CK(p_cuEventSynchronize(t1));
    float ms = 0; CK(p_cuEventElapsedTime(&ms, t0, t1));
    const double msApply = ms / iters;
    // Useful traffic: u + y-write + per-point metric G (dominant, like the hand
    // kernel's d_G stream).
    const double bytesElem = (2.0 * n3 + gElem) * 8.0;
    printf("  %d applies: %.3f ms/apply | %.1f ns/elem | ~%.0f GB/s useful (%.1f KiB/elem)\n",
           iters, msApply, msApply * 1e6 / (double)E,
           (double)E * bytesElem / (msApply * 1e-3) / 1e9, bytesElem / 1024.0);
    return err < 1e-11 ? 0 : 1;
}
