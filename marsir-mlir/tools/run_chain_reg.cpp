// GPU gate for the register-resident mma chain (warp_chain_reg_sm90.ptx): two
// m8n8k4 mma where mma2's A operand is a gpu.shuffle RELAYOUT of mma1's C
// fragment -- no shared memory, no barrier (the foundation of register-resident
// tensor-core chaining). Validates the shuffle formula numerically:
//   D1 = C + A1@B1 ; A2 = D1[:,0:4] ; out = D1 + A2@B2
// A1,B1,B2 random; C random. Driver-API only. grid 1, block (32,1,1).
// Build: g++ -O2 run_chain_reg.cpp -o run_chain_reg -ldl
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
  const char*ptx=argc>1?argv[1]:"../generated/warp_chain_reg_sm90.ptx";
  FILE*f=fopen(ptx,"rb"); if(!f){perror(ptx);return 1;} fseek(f,0,SEEK_END);long L=ftell(f);fseek(f,0,SEEK_SET);
  std::vector<char> buf(L+1,0); if(fread(buf.data(),1,L,f)!=(size_t)L){fclose(f);return 1;} fclose(f);
  CK(p_cuInit(0)); CUdevice d; CK(p_cuDeviceGet(&d,0)); CUcontext c; CK(p_cuCtxCreate(&c,0,d));
  CUmodule m; CK(p_cuModuleLoadData(&m,buf.data())); CUfunction fn; CK(p_cuModuleGetFunction(&fn,m,"chain_mma"));
  double A1[32],B1[32],B2[32],C[64],out[64],ref[64]; srand(42);
  for(auto&x:A1)x=2.0*rand()/RAND_MAX-1; for(auto&x:B1)x=2.0*rand()/RAND_MAX-1;
  for(auto&x:B2)x=2.0*rand()/RAND_MAX-1; for(auto&x:C)x=2.0*rand()/RAND_MAX-1;
  double D1[64]; for(int i=0;i<8;++i)for(int j=0;j<8;++j){double s=C[i*8+j];for(int k=0;k<4;++k)s+=A1[i*4+k]*B1[k*8+j];D1[i*8+j]=s;}
  for(int i=0;i<8;++i)for(int j=0;j<8;++j){double s=D1[i*8+j];for(int k=0;k<4;++k)s+=D1[i*8+k]*B2[k*8+j];ref[i*8+j]=s;} // A2=D1[:,0:4]
  CUdeviceptr dA1,dB1,dB2,dC;
  CK(p_cuMemAlloc(&dA1,32*8));CK(p_cuMemAlloc(&dB1,32*8));CK(p_cuMemAlloc(&dB2,32*8));CK(p_cuMemAlloc(&dC,64*8));
  CK(p_cuMemcpyHtoD(dA1,A1,32*8));CK(p_cuMemcpyHtoD(dB1,B1,32*8));CK(p_cuMemcpyHtoD(dB2,B2,32*8));CK(p_cuMemcpyHtoD(dC,C,64*8));
  long long z=0, a4Sz[2]={8,4},a4St[2]={4,1}, b4Sz[2]={4,8},b4St[2]={8,1}, cSz[2]={8,8},cSt[2]={8,1};
  void*args[28]={&dA1,&dA1,&z,&a4Sz[0],&a4Sz[1],&a4St[0],&a4St[1],
                 &dB1,&dB1,&z,&b4Sz[0],&b4Sz[1],&b4St[0],&b4St[1],
                 &dB2,&dB2,&z,&b4Sz[0],&b4Sz[1],&b4St[0],&b4St[1],
                 &dC,&dC,&z,&cSz[0],&cSz[1],&cSt[0],&cSt[1]};
  CK(p_cuLaunchKernel(fn,1,1,1,32,1,1,0,0,args,0)); CK(p_cuCtxSynchronize()); CK(p_cuMemcpyDtoH(out,dC,64*8));
  double e=0; for(int i=0;i<64;++i)e=fmax(e,fabs(out[i]-ref[i]));
  printf("register-resident mma chain (shuffle relayout) gate: max|out-ref| = %.3e  %s\n",e,e<1e-12?"PASS":"FAIL");
  return e<1e-12?0:1;
}
