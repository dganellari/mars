// GPU gate for the PASS-GENERATED register-resident chain
// (warp_chain_pass_sm90.ptx). The INPUT to mir-opt was high-level 8x8
// vector.contract ops -- no fragments, no shuffles, no warp regions. The
// --mir-chain-contracts PASS invented the per-lane mma + shuffle relayout.
// This validates that a real MLIR lowering (not a Python emitter) produces
// correct register-resident tensor-core code:   Out = (A @ B1) @ B2
// Build: g++ -O2 run_chain_pass.cpp -o run_chain_pass -ldl
#include <dlfcn.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
typedef int CUresult; typedef int CUdevice; typedef unsigned long long CUdeviceptr;
typedef struct CUctx_st* CUcontext; typedef struct CUmod_st* CUmodule;
typedef struct CUfunc_st* CUfunction; typedef struct CUstream_st* CUstream;
static CUresult(*p_cuInit)(unsigned);static CUresult(*p_cuDeviceGet)(CUdevice*,int);
static CUresult(*p_cuCtxCreate)(CUcontext*,unsigned,CUdevice);static CUresult(*p_cuModuleLoadData)(CUmodule*,const void*);
static CUresult(*p_cuModuleGetFunction)(CUfunction*,CUmodule,const char*);static CUresult(*p_cuMemAlloc)(CUdeviceptr*,size_t);
static CUresult(*p_cuMemcpyHtoD)(CUdeviceptr,const void*,size_t);static CUresult(*p_cuMemcpyDtoH)(void*,CUdeviceptr,size_t);
static CUresult(*p_cuLaunchKernel)(CUfunction,unsigned,unsigned,unsigned,unsigned,unsigned,unsigned,unsigned,CUstream,void**,void**);
static CUresult(*p_cuCtxSynchronize)(void);static CUresult(*p_cuGetErrorString)(CUresult,const char**);
#define CK(x) do{CUresult r_=(x);if(r_){const char*s_=0;if(p_cuGetErrorString)p_cuGetErrorString(r_,&s_);fprintf(stderr,"CUDA err %d (%s) %s:%d\n",r_,s_?s_:"?",__FILE__,__LINE__);exit(1);}}while(0)
static void* sy(void*l,const char*n){void*p=dlsym(l,n);if(!p){fprintf(stderr,"missing %s\n",n);exit(1);}return p;}
int main(int argc,char**argv){
  void*lib=dlopen("libcuda.so.1",RTLD_NOW); if(!lib)lib=dlopen("libcuda.so",RTLD_NOW);
  if(!lib){fprintf(stderr,"no libcuda\n");return 1;}
  *(void**)&p_cuInit=sy(lib,"cuInit");*(void**)&p_cuDeviceGet=sy(lib,"cuDeviceGet");*(void**)&p_cuCtxCreate=sy(lib,"cuCtxCreate_v2");
  *(void**)&p_cuModuleLoadData=sy(lib,"cuModuleLoadData");*(void**)&p_cuModuleGetFunction=sy(lib,"cuModuleGetFunction");
  *(void**)&p_cuMemAlloc=sy(lib,"cuMemAlloc_v2");*(void**)&p_cuMemcpyHtoD=sy(lib,"cuMemcpyHtoD_v2");*(void**)&p_cuMemcpyDtoH=sy(lib,"cuMemcpyDtoH_v2");
  *(void**)&p_cuLaunchKernel=sy(lib,"cuLaunchKernel");*(void**)&p_cuCtxSynchronize=sy(lib,"cuCtxSynchronize");*(void**)&p_cuGetErrorString=sy(lib,"cuGetErrorString");
  const char*ptx=argc>1?argv[1]:"../generated/warp_chain_pass_sm90.ptx";
  FILE*f=fopen(ptx,"rb"); if(!f){perror(ptx);return 1;} fseek(f,0,SEEK_END);long L=ftell(f);fseek(f,0,SEEK_SET);
  std::vector<char> buf(L+1,0); if(fread(buf.data(),1,L,f)!=(size_t)L){fclose(f);return 1;} fclose(f);
  CK(p_cuInit(0)); CUdevice d; CK(p_cuDeviceGet(&d,0)); CUcontext c; CK(p_cuCtxCreate(&c,0,d));
  CUmodule m; CK(p_cuModuleLoadData(&m,buf.data())); CUfunction fn; CK(p_cuModuleGetFunction(&fn,m,"chain2"));
  const int n=8,nn=64; double A[64],B1[64],B2[64],T[64],ref[64],out[64]; srand(42);
  for(auto&x:A)x=2.0*rand()/RAND_MAX-1; for(auto&x:B1)x=2.0*rand()/RAND_MAX-1; for(auto&x:B2)x=2.0*rand()/RAND_MAX-1;
  for(int i=0;i<n;++i)for(int j=0;j<n;++j){double s=0;for(int k=0;k<n;++k)s+=A[i*n+k]*B1[k*n+j];T[i*n+j]=s;}
  for(int i=0;i<n;++i)for(int j=0;j<n;++j){double s=0;for(int k=0;k<n;++k)s+=T[i*n+k]*B2[k*n+j];ref[i*n+j]=s;}
  CUdeviceptr dA,dB1,dB2,dO;
  CK(p_cuMemAlloc(&dA,nn*8));CK(p_cuMemAlloc(&dB1,nn*8));CK(p_cuMemAlloc(&dB2,nn*8));CK(p_cuMemAlloc(&dO,nn*8));
  CK(p_cuMemcpyHtoD(dA,A,nn*8));CK(p_cuMemcpyHtoD(dB1,B1,nn*8));CK(p_cuMemcpyHtoD(dB2,B2,nn*8));
  double zero[64]={0}; CK(p_cuMemcpyHtoD(dO,zero,nn*8));
  long long z=0,sz[2]={8,8},st[2]={8,1}; void*args[28];
  CUdeviceptr* ps[4]={&dA,&dB1,&dB2,&dO};
  for(int i=0;i<4;++i){void**a=&args[i*7];a[0]=ps[i];a[1]=ps[i];a[2]=&z;a[3]=&sz[0];a[4]=&sz[1];a[5]=&st[0];a[6]=&st[1];}
  CK(p_cuLaunchKernel(fn,1,1,1,32,1,1,0,0,args,0)); CK(p_cuCtxSynchronize()); CK(p_cuMemcpyDtoH(out,dO,nn*8));
  double e=0; for(int i=0;i<nn;++i)e=fmax(e,fabs(out[i]-ref[i]));
  printf("PASS-generated register-resident chain gate: max|out-(A@B1)@B2| = %.3e  %s\n",e,e<1e-12?"PASS":"FAIL");
  return e<1e-12?0:1;
}
