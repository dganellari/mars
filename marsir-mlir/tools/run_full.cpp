// GPU gate for the EMITTER-GENERATED FULL 3-direction HO operator
// (warp_full_sm90.ptx): all three directions (each B-sweep + 7 faces + plane
// scatter) accumulating into y, warp-per-element, staged in shared memory (264
// mma, 10.5 KB shared/block). Correctness form: the input is presented in the
// three moveaxis layouts u0/u1/u2 (transpose done host-side for now). Metric is
// 9 per-(dir,component) arrays. Matches the full Knaus oracle.
// ABI: u0,u1,u2(8x64), Btil,Dtil,Dm,W(8x8), g{d}{c}all(8x8x8) x9, y3(8x8x8)
//      = 3*7 + 4*7 + 9*9 + 9 = 139 params. One warp; driver-API only.
// Build:  g++ -O2 run_full.cpp -o run_full -ldl
// Run:    ./run_full ../generated/warp_full_sm90.ptx

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
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) {
        double s = 0; for (int k = 0; k < n; ++k) s += A[i*n+k]*B[k*n+j];
        C[i*n+j] = s;
    }
}
static void mmT(const double* A, const double* B, double* C, int n)  // C=A@B^T
{
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) {
        double s = 0; for (int k = 0; k < n; ++k) s += A[i*n+k]*B[j*n+k];
        C[i*n+j] = s;
    }
}

int main(int argc, char** argv)
{
    void* lib = dlopen("libcuda.so.1", RTLD_NOW);
    if (!lib) lib = dlopen("libcuda.so", RTLD_NOW);
    if (!lib) { fprintf(stderr, "cannot dlopen libcuda.so.1: %s\n", dlerror()); return 1; }
    *(void**)&p_cuInit = must_sym(lib, "cuInit");
    *(void**)&p_cuDeviceGet = must_sym(lib, "cuDeviceGet");
    *(void**)&p_cuCtxCreate = must_sym(lib, "cuCtxCreate_v2");
    *(void**)&p_cuModuleLoadData = must_sym(lib, "cuModuleLoadData");
    *(void**)&p_cuModuleGetFunction = must_sym(lib, "cuModuleGetFunction");
    *(void**)&p_cuMemAlloc = must_sym(lib, "cuMemAlloc_v2");
    *(void**)&p_cuMemcpyHtoD = must_sym(lib, "cuMemcpyHtoD_v2");
    *(void**)&p_cuMemcpyDtoH = must_sym(lib, "cuMemcpyDtoH_v2");
    *(void**)&p_cuLaunchKernel = must_sym(lib, "cuLaunchKernel");
    *(void**)&p_cuCtxSynchronize = must_sym(lib, "cuCtxSynchronize");
    *(void**)&p_cuGetErrorString = must_sym(lib, "cuGetErrorString");

    const char* ptxPath = argc > 1 ? argv[1] : "../generated/warp_full_sm90.ptx";
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
    CUfunction fn; CK(p_cuModuleGetFunction(&fn, mod, "full"));

    const int n = 8, P = 7, nn = 64, nnn = 512;
    std::vector<double> u(nnn), Btil(nn), Dtil(nn), Dm(nn), W(nn), y(nnn, 0.0),
                        ref(nnn, 0.0), out(nnn);
    std::vector<std::vector<double>> g(9, std::vector<double>(nnn));
    srand(42);
    for (auto* v : {&u, &Btil, &Dtil, &Dm, &W})
        for (auto& x : *v) x = 2.0 * rand() / RAND_MAX - 1.0;
    for (auto& gg : g) for (auto& x : gg) x = 2.0 * rand() / RAND_MAX - 1.0;

    // moveaxis input layouts (u[i,j,k] = u[i*64+j*8+k])
    std::vector<double> u0(8*nn), u1(8*nn), u2(8*nn);
    for (int q = 0; q < n; ++q) for (int s = 0; s < n; ++s) for (int r = 0; r < n; ++r) {
        u0[q*nn + s*n + r] = u[q*nn + s*n + r];  // u[q,s,r]
        u1[q*nn + s*n + r] = u[s*nn + q*n + r];  // u[s,q,r]
        u2[q*nn + s*n + r] = u[s*nn + r*n + q];  // u[s,r,q]
    }
    const double* uf[3] = {u0.data(), u1.data(), u2.data()};

    // oracle
    for (int d = 0; d < 3; ++d) {
        std::vector<double> ia(8*nn), da(8*nn);  // B-sweeps (rectangular 8x8 @ 8x64)
        for (int i = 0; i < n; ++i) for (int j = 0; j < nn; ++j) {
            double si = 0, sd = 0;
            for (int k = 0; k < n; ++k) { si += Btil[i*n+k]*uf[d][k*nn+j];
                                          sd += Dtil[i*n+k]*uf[d][k*nn+j]; }
            ia[i*nn+j] = si; da[i*nn+j] = sd;
        }
        const double* g0 = g[d*3+0].data(), *g1 = g[d*3+1].data(), *g2 = g[d*3+2].data();
        for (int l = 0; l < P; ++l) {
            const double* interp = &ia[l*nn]; const double* deriv = &da[l*nn];
            double dt2[64], dt1[64], flux[64], tmp[64], intf[64];
            mmT(interp, Dm.data(), dt2, n);
            mm(Dm.data(), interp, dt1, n);
            for (int i = 0; i < nn; ++i)
                flux[i] = g2[l*nn+i]*deriv[i] + g0[l*nn+i]*dt2[i] + g1[l*nn+i]*dt1[i];
            mmT(flux, W.data(), tmp, n);
            mm(W.data(), tmp, intf, n);
            for (int s = 0; s < n; ++s) for (int r = 0; r < n; ++r) {
                int a, b;  // linear positions of Y[l][s,r] and Y[l+1][s,r]
                if (d == 0) { a = l*nn + s*n + r;       b = (l+1)*nn + s*n + r; }
                else if (d == 1) { a = s*nn + l*n + r;  b = s*nn + (l+1)*n + r; }
                else { a = s*nn + r*n + l;              b = s*nn + r*n + (l+1); }
                ref[a] -= intf[s*n+r]; ref[b] += intf[s*n+r];
            }
        }
    }

    // device buffers, in kernel arg order
    auto alloc = [&](const double* h, int cnt) {
        CUdeviceptr d; CK(p_cuMemAlloc(&d, cnt*sizeof(double)));
        CK(p_cuMemcpyHtoD(d, h, cnt*sizeof(double))); return d;
    };
    std::vector<CUdeviceptr> dev_ptrs;
    dev_ptrs.push_back(alloc(u0.data(), 8*nn));
    dev_ptrs.push_back(alloc(u1.data(), 8*nn));
    dev_ptrs.push_back(alloc(u2.data(), 8*nn));
    dev_ptrs.push_back(alloc(Btil.data(), nn));
    dev_ptrs.push_back(alloc(Dtil.data(), nn));
    dev_ptrs.push_back(alloc(Dm.data(), nn));
    dev_ptrs.push_back(alloc(W.data(), nn));
    for (int i = 0; i < 9; ++i) dev_ptrs.push_back(alloc(g[i].data(), nnn));
    CUdeviceptr dy = alloc(y.data(), nnn);
    dev_ptrs.push_back(dy);

    long long z = 0;
    long long u2Sz[2] = {8, 64}, u2St[2] = {64, 1};
    long long sz2[2] = {8, 8}, st2[2] = {8, 1};
    long long sz3[3] = {8, 8, 8}, st3[3] = {64, 8, 1};
    std::vector<void*> args;
    auto add2 = [&](CUdeviceptr* p, long long* sz, long long* st) {
        args.push_back(p); args.push_back(p); args.push_back(&z);
        args.push_back(&sz[0]); args.push_back(&sz[1]); args.push_back(&st[0]); args.push_back(&st[1]);
    };
    auto add3 = [&](CUdeviceptr* p) {
        args.push_back(p); args.push_back(p); args.push_back(&z);
        args.push_back(&sz3[0]); args.push_back(&sz3[1]); args.push_back(&sz3[2]);
        args.push_back(&st3[0]); args.push_back(&st3[1]); args.push_back(&st3[2]);
    };
    add2(&dev_ptrs[0], u2Sz, u2St); add2(&dev_ptrs[1], u2Sz, u2St); add2(&dev_ptrs[2], u2Sz, u2St);
    for (int i = 3; i < 7; ++i) add2(&dev_ptrs[i], sz2, st2);
    for (int i = 7; i < 16; ++i) add3(&dev_ptrs[i]);
    add3(&dev_ptrs[16]);

    CK(p_cuLaunchKernel(fn, 1, 1, 1, 32, 1, 1, 0, 0, args.data(), nullptr));
    CK(p_cuCtxSynchronize());
    CK(p_cuMemcpyDtoH(out.data(), dy, nnn*sizeof(double)));

    double err = 0;
    for (int i = 0; i < nnn; ++i) err = fmax(err, fabs(out[i] - ref[i]));
    printf("emitter FULL operator gate: max|y - Knaus(u)| = %.3e  %s\n",
           err, err < 1e-12 ? "PASS" : "FAIL");
    return err < 1e-12 ? 0 : 1;
}
