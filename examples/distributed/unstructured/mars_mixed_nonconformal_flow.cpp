#include "mars_mixed_nonconformal_runtime.hpp"
#include "mars_mixed_nonconformal_linear.hpp"
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);MPI_Comm_rank(MPI_COMM_WORLD,&rankId);MPI_Comm_size(MPI_COMM_WORLD,&rankCount);int result=1;
    try {
#if defined(__CUDACC__) || defined(__HIPCC__)
        MPI_Comm local;MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rankId,MPI_INFO_NULL,&local);int local_rank,devices;
        MPI_Comm_rank(local,&local_rank);check(cudaGetDeviceCount(&devices)==cudaSuccess && devices>0,"mixed GPU selection");
        check(cudaSetDevice(devices==1?0:local_rank%devices)==cudaSuccess,"mixed selected GPU");MPI_Comm_free(&local);
#endif
        std::string path,output;int budget=2000;double tolerance=1e-9;bool setup_only=false,high_resolution=true;
        for(int i=1;i<argc;++i) {std::string arg=argv[i];
            if(arg=="--partitioned-mesh" && i+1<argc) path=argv[++i];else if(arg=="--output" && i+1<argc) output=argv[++i];
            else if(arg=="--iterations" && i+1<argc) budget=std::stoi(argv[++i]);else if(arg=="--tolerance" && i+1<argc) tolerance=std::stod(argv[++i]);
            else if(arg=="--advection" && i+1<argc) {
                const std::string scheme=argv[++i];
                if(scheme!="upwind" && scheme!="high-resolution") throw std::runtime_error("advection must be upwind or high-resolution");
                high_resolution=scheme=="high-resolution";
            }
            else if(arg=="--setup-only") setup_only=true;else throw std::runtime_error("unknown or incomplete mixed flow option");}
        if(path.empty() || budget<1 || budget>100000 || !std::isfinite(tolerance) || tolerance<=0) throw std::runtime_error("invalid mixed flow options");
        mixedflow::hypre_check(HYPRE_Initialize());
#if defined(__CUDACC__) || defined(__HIPCC__)
        mixedflow::hypre_check(HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE));mixedflow::hypre_check(HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE));mixedflow::hypre_check(HYPRE_SetGpuAwareMPI(1));
#endif
        {auto mesh=mixedflow::read(path);mixedflow::Channel channel(mesh);channel.view_.controls.high_resolution=high_resolution;
         if(rankId==0) std::printf("Advection: %s\n",high_resolution?"high-resolution":"upwind");
         bool ok=setup_only || channel.solve(budget,tolerance,output);if(rankId==0) std::printf("%s: mixed %s\n",ok?"PASS":"FAIL",setup_only?"setup only":"flow convergence (field comparison separate)");result=ok?0:2;}
        mixedflow::hypre_check(HYPRE_Finalize());
    } catch(const std::exception& error) {std::fprintf(stderr,"rank %d mixed flow: %s\n",rankId,error.what());MPI_Abort(MPI_COMM_WORLD,1);}
    MPI_Finalize();return result;
}
