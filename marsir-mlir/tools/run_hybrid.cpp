// Throughput measurement for the EMITTER-GENERATED batched HO operator with a
// CHEAP per-element metric (warp_hybrid_sm90.ptx): g{d}{c} are Exnxn
// (one nxn/elem, shared across faces) -> ~8x less metric traffic, arithmetic
// intensity above the FP64 ridge (compute-bound regime).: one warp per element, grid = E blocks. Spot-checks
// element 0 against the full Knaus oracle, then times the kernel over E
// elements and reports ms, elements/s, DOF/s, and an effective GB/s.
//
// This is the FP64 tensor-core number for the warp-per-element structure -- to
// be read against the prior hand-DMMA result (1.06x vs a matched scalar twin,
// LSU-instruction-bound). A matched scalar twin from the same emitter is the
// clean apples-to-apples follow-up; this establishes where the TC kernel lands.
//
// ABI: U(?x8x8x8)=11, Btil,Dtil,Dm,W(8x8)=7 each, 9x metric(?x8x8x8)=11 each,
//      Y(?x8x8x8)=11 -> 149 params. Launch grid=(E,1,1) block=(32,1,1).
// Build:  g++ -O2 run_batched.cpp -o run_batched -ldl
// Run:    ./run_batched ../generated/warp_batched_sm90.ptx [E]

#include <dlfcn.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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
static CUresult (*p_cuMemsetD8)(CUdeviceptr, unsigned char, size_t);
static CUresult (*p_cuMemcpyHtoD)(CUdeviceptr, const void*, size_t);
static CUresult (*p_cuMemcpyDtoD)(CUdeviceptr, CUdeviceptr, size_t);
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

#define CK(x) do { CUresult r_ = (x); if (r_ != 0) { \
    const char* s_ = nullptr; if (p_cuGetErrorString) p_cuGetErrorString(r_, &s_); \
    fprintf(stderr, "CUDA error %d (%s) at %s:%d\n", r_, s_ ? s_ : "?", __FILE__, __LINE__); \
    exit(1); } } while (0)

static void* sy(void* lib, const char* n)
{ void* p = dlsym(lib, n); if (!p) { fprintf(stderr, "missing %s\n", n); exit(1); } return p; }

static void mm(const double* A, const double* B, double* C, int n)
{ for (int i=0;i<n;++i) for (int j=0;j<n;++j){double s=0;for(int k=0;k<n;++k)s+=A[i*n+k]*B[k*n+j];C[i*n+j]=s;} }
static void mmT(const double* A, const double* B, double* C, int n)
{ for (int i=0;i<n;++i) for (int j=0;j<n;++j){double s=0;for(int k=0;k<n;++k)s+=A[i*n+k]*B[j*n+k];C[i*n+j]=s;} }

int main(int argc, char** argv)
{
    void* lib = dlopen("libcuda.so.1", RTLD_NOW);
    if (!lib) lib = dlopen("libcuda.so", RTLD_NOW);
    if (!lib) { fprintf(stderr, "cannot dlopen libcuda: %s\n", dlerror()); return 1; }
    *(void**)&p_cuInit = sy(lib,"cuInit");
    *(void**)&p_cuDeviceGet = sy(lib,"cuDeviceGet");
    *(void**)&p_cuCtxCreate = sy(lib,"cuCtxCreate_v2");
    *(void**)&p_cuModuleLoadData = sy(lib,"cuModuleLoadData");
    *(void**)&p_cuModuleGetFunction = sy(lib,"cuModuleGetFunction");
    *(void**)&p_cuMemAlloc = sy(lib,"cuMemAlloc_v2");
    *(void**)&p_cuMemsetD8 = sy(lib,"cuMemsetD8_v2");
    *(void**)&p_cuMemcpyHtoD = sy(lib,"cuMemcpyHtoD_v2");
    *(void**)&p_cuMemcpyDtoD = sy(lib,"cuMemcpyDtoD_v2");
    *(void**)&p_cuMemcpyDtoH = sy(lib,"cuMemcpyDtoH_v2");
    *(void**)&p_cuLaunchKernel = sy(lib,"cuLaunchKernel");
    *(void**)&p_cuCtxSynchronize = sy(lib,"cuCtxSynchronize");
    *(void**)&p_cuEventCreate = sy(lib,"cuEventCreate");
    *(void**)&p_cuEventRecord = sy(lib,"cuEventRecord");
    *(void**)&p_cuEventSynchronize = sy(lib,"cuEventSynchronize");
    *(void**)&p_cuEventElapsedTime = sy(lib,"cuEventElapsedTime");
    *(void**)&p_cuGetErrorString = sy(lib,"cuGetErrorString");

    const char* ptxPath = argc > 1 ? argv[1] : "../generated/warp_hybrid_sm90.ptx";
    long E = argc > 2 ? atol(argv[2]) : 65536;
    FILE* f = fopen(ptxPath, "rb");
    if (!f) { perror(ptxPath); return 1; }
    fseek(f,0,SEEK_END); long len=ftell(f); fseek(f,0,SEEK_SET);
    std::vector<char> ptx(len+1,0);
    if (fread(ptx.data(),1,len,f)!=(size_t)len){fclose(f);return 1;} fclose(f);

    CK(p_cuInit(0));
    CUdevice dev; CK(p_cuDeviceGet(&dev,0));
    CUcontext ctx; CK(p_cuCtxCreate(&ctx,0,dev));
    CUmodule mod; CK(p_cuModuleLoadData(&mod, ptx.data()));
    CUfunction fn; CK(p_cuModuleGetFunction(&fn, mod, "full_hybrid"));

    const int n=8, P=7, nn=64, nnn=512;
    // element-0 data + oracle
    std::vector<double> u(nnn), Btil(nn), Dtil(nn), Dm(nn), W(nn), ref(nnn,0.0);
    std::vector<std::vector<double>> g(9, std::vector<double>(nn));  // per-element nxn
    srand(42);
    for (auto* v : {&u,&Btil,&Dtil,&Dm,&W}) for (auto& x:*v) x=2.0*rand()/RAND_MAX-1.0;
    for (auto& gg:g) for (auto& x:gg) x=2.0*rand()/RAND_MAX-1.0;
    std::vector<double> u0(8*nn),u1(8*nn),u2(8*nn);
    for(int q=0;q<n;++q)for(int s=0;s<n;++s)for(int r=0;r<n;++r){
        u0[q*nn+s*n+r]=u[q*nn+s*n+r]; u1[q*nn+s*n+r]=u[s*nn+q*n+r]; u2[q*nn+s*n+r]=u[s*nn+r*n+q]; }
    const double* uf[3]={u0.data(),u1.data(),u2.data()};
    for(int d=0;d<3;++d){
        std::vector<double> ia(8*nn),da(8*nn);
        for(int i=0;i<n;++i)for(int j=0;j<nn;++j){double si=0,sd=0;
            for(int k=0;k<n;++k){si+=Btil[i*n+k]*uf[d][k*nn+j];sd+=Dtil[i*n+k]*uf[d][k*nn+j];}
            ia[i*nn+j]=si;da[i*nn+j]=sd;}
        const double *g0=g[d*3+0].data(),*g1=g[d*3+1].data(),*g2=g[d*3+2].data();
        for(int l=0;l<P;++l){const double* ip=&ia[l*nn];const double* dv=&da[l*nn];
            double dt2[64],dt1[64],flux[64],tmp[64],intf[64];
            mmT(ip,Dm.data(),dt2,n); mm(Dm.data(),ip,dt1,n);
            for(int i=0;i<nn;++i)flux[i]=g2[i]*dv[i]+g0[i]*dt2[i]+g1[i]*dt1[i];  // per-element metric, all faces
            mmT(flux,W.data(),tmp,n); mm(W.data(),tmp,intf,n);
            for(int s=0;s<n;++s)for(int r=0;r<n;++r){int a,b;
                if(d==0){a=l*nn+s*n+r;b=(l+1)*nn+s*n+r;}
                else if(d==1){a=s*nn+l*n+r;b=s*nn+(l+1)*n+r;}
                else{a=s*nn+r*n+l;b=s*nn+r*n+(l+1);}
                ref[a]-=intf[s*n+r];ref[b]+=intf[s*n+r];}}
    }

    // device: batched U/metric/Y (E elements), broadcast element-0 data.
    size_t elemB = nnn*sizeof(double);          // u/y element = 8x8x8
    size_t metB  = nn*sizeof(double);            // metric element = nxn
    auto allocBatchedSz = [&](const double* elem0, size_t eb)->CUdeviceptr {
        CUdeviceptr d; CK(p_cuMemAlloc(&d, (size_t)E*eb));
        CK(p_cuMemcpyHtoD(d, elem0, eb));
        for (long e=1;e<E;++e) CK(p_cuMemcpyDtoD(d+(size_t)e*eb, d, eb));
        return d;
    };
    CUdeviceptr dU = allocBatchedSz(u.data(), elemB);
    CUdeviceptr dG[9]; for(int i=0;i<9;++i) dG[i]=allocBatchedSz(g[i].data(), metB);
    CUdeviceptr dY; CK(p_cuMemAlloc(&dY,(size_t)E*elemB)); CK(p_cuMemsetD8(dY,0,(size_t)E*elemB));
    auto allocOp = [&](const double* h)->CUdeviceptr {
        CUdeviceptr d; CK(p_cuMemAlloc(&d,nn*sizeof(double))); CK(p_cuMemcpyHtoD(d,h,nn*sizeof(double))); return d; };
    CUdeviceptr dBt=allocOp(Btil.data()),dDt=allocOp(Dtil.data()),dDm=allocOp(Dm.data()),dW=allocOp(W.data());

    long long z=0;
    long long b4Sz[4]={E,8,8,8}, b4St[4]={512,64,8,1};
    long long m3Sz[3]={E,8,8}, m3St[3]={64,8,1};
    long long sz2[2]={8,8}, st2[2]={8,1};
    std::vector<void*> args;
    auto add4=[&](CUdeviceptr* p){ args.push_back(p);args.push_back(p);args.push_back(&z);
        for(int i=0;i<4;++i)args.push_back(&b4Sz[i]); for(int i=0;i<4;++i)args.push_back(&b4St[i]); };
    auto add3m=[&](CUdeviceptr* p){ args.push_back(p);args.push_back(p);args.push_back(&z);
        for(int i=0;i<3;++i)args.push_back(&m3Sz[i]); for(int i=0;i<3;++i)args.push_back(&m3St[i]); };
    auto add2=[&](CUdeviceptr* p){ args.push_back(p);args.push_back(p);args.push_back(&z);
        args.push_back(&sz2[0]);args.push_back(&sz2[1]);args.push_back(&st2[0]);args.push_back(&st2[1]); };
    add4(&dU); add2(&dBt); add2(&dDt); add2(&dDm); add2(&dW);
    for(int i=0;i<9;++i) add3m(&dG[i]); add4(&dY);

    // correctness spot-check (element 0)
    CK(p_cuLaunchKernel(fn,(unsigned)E,1,1,32,1,1,0,0,args.data(),nullptr));
    CK(p_cuCtxSynchronize());
    std::vector<double> out0(nnn);
    CK(p_cuMemcpyDtoH(out0.data(), dY, elemB));
    double err=0; for(int i=0;i<nnn;++i) err=fmax(err,fabs(out0[i]-ref[i]));
    printf("HYBRID (register-resident face chain) operator: E=%ld  spot-check(elem0) max|y-Knaus| = %.3e  %s\n",
           E, err, err<1e-12?"PASS":"FAIL");

    // timing
    CUevent t0,t1; CK(p_cuEventCreate(&t0,0)); CK(p_cuEventCreate(&t1,0));
    const int iters=20;
    for(int w=0;w<3;++w) CK(p_cuLaunchKernel(fn,(unsigned)E,1,1,32,1,1,0,0,args.data(),nullptr));
    CK(p_cuCtxSynchronize());
    CK(p_cuEventRecord(t0,0));
    for(int it=0;it<iters;++it) CK(p_cuLaunchKernel(fn,(unsigned)E,1,1,32,1,1,0,0,args.data(),nullptr));
    CK(p_cuEventRecord(t1,0)); CK(p_cuEventSynchronize(t1));
    float ms=0; CK(p_cuEventElapsedTime(&ms,t0,t1)); ms/=iters;
    double dofs = (double)E*nnn;                         // (p+1)^3 DOF/elem
    double bytes = (double)E*(2*elemB + 9*metB);         // u + y + 9 per-elem metric
    printf("  %.3f ms/apply  |  %.2f Melem/s  |  %.2f GDOF/s  |  ~%.0f GB/s (u+metric+y)\n",
           ms, E/ms/1e3, dofs/ms/1e6, bytes/ms/1e6);
    return err<1e-12?0:1;
}
