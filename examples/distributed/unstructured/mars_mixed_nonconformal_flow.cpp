#include "mars_mixed_nonconformal_runtime.hpp"
#include "mars_mixed_nonconformal_linear.hpp"
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
#ifdef MARS_MIXED_PETSC
    mixedflow::petsc_check(PetscInitializeNoArguments());
#endif
    MPI_Comm_rank(MPI_COMM_WORLD,&rankId);MPI_Comm_size(MPI_COMM_WORLD,&rankCount);int result=1;
    try {
#if defined(__CUDACC__) || defined(__HIPCC__)
        MPI_Comm local;MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rankId,MPI_INFO_NULL,&local);int local_rank,devices;
        MPI_Comm_rank(local,&local_rank);check(cudaGetDeviceCount(&devices)==cudaSuccess && devices>0,"mixed GPU selection");
        const int selected=devices==1?0:local_rank%devices;
        check(cudaSetDevice(selected)==cudaSuccess,"mixed selected GPU");MPI_Comm_free(&local);
        char bus[32]{};check(cudaDeviceGetPCIBusId(bus,sizeof(bus),selected)==cudaSuccess,"mixed GPU identity");
        std::printf("GPU assembly: rank %d device %d PCI %s\n",rankId,selected,bus);
#endif
        std::string path,output,interface_name="compatibility";int budget=2000;double tolerance=1e-9;bool setup_only=false,high_resolution=true;
        mixedflow::FlowProfile profile=mixedflow::FlowProfile::mars;
        mixedflow::InterfaceMethod interface_method=mixedflow::InterfaceMethod::compatibility;
        for(int i=1;i<argc;++i) {std::string arg=argv[i];
            if(arg=="--partitioned-mesh" && i+1<argc) path=argv[++i];else if(arg=="--output" && i+1<argc) output=argv[++i];
            else if(arg=="--iterations" && i+1<argc) budget=std::stoi(argv[++i]);else if(arg=="--tolerance" && i+1<argc) tolerance=std::stod(argv[++i]);
            else if(arg=="--profile" && i+1<argc) {
                const std::string name=argv[++i];
                if(name=="mars") profile=mixedflow::FlowProfile::mars;
                else if(name=="openaccel") profile=mixedflow::FlowProfile::openaccel;
                else throw std::runtime_error("profile must be mars or openaccel");
            }
            else if(arg=="--interface" && i+1<argc) {
                interface_name=argv[++i];
                if(interface_name=="compatibility") interface_method=mixedflow::InterfaceMethod::compatibility;
                else if(interface_name=="conservative") interface_method=mixedflow::InterfaceMethod::conservative;
                else if(interface_name=="legacy") interface_method=mixedflow::InterfaceMethod::legacy;
                else throw std::runtime_error("interface must be compatibility, conservative or legacy");
            }
            else if(arg=="--advection" && i+1<argc) {
                const std::string scheme=argv[++i];
                if(scheme!="upwind" && scheme!="high-resolution") throw std::runtime_error("advection must be upwind or high-resolution");
                high_resolution=scheme=="high-resolution";
            }
            else if(arg=="--setup-only") setup_only=true;else throw std::runtime_error("unknown or incomplete mixed flow option");}
        if(path.empty() || budget<1 || budget>100000 || !std::isfinite(tolerance) || tolerance<=0) throw std::runtime_error("invalid mixed flow options");
#ifndef MARS_MIXED_PETSC
        if(profile==mixedflow::FlowProfile::openaccel) throw std::runtime_error("OpenAccel profile requires MARS_ENABLE_PETSC");
#endif
        if(profile==mixedflow::FlowProfile::openaccel && interface_method!=mixedflow::InterfaceMethod::compatibility)
            throw std::runtime_error("OpenAccel profile requires --interface compatibility");
#if (defined(__CUDACC__) || defined(__HIPCC__)) && !defined(HYPRE_USING_GPU)
        if(profile!=mixedflow::FlowProfile::openaccel) throw std::runtime_error("GPU MARS profile requires GPU-enabled Hypre");
#endif
        mixedflow::hypre_check(HYPRE_Initialize());
#if defined(__CUDACC__) || defined(__HIPCC__)
        mixedflow::hypre_check(HYPRE_SetMemoryLocation(profile==mixedflow::FlowProfile::openaccel?HYPRE_MEMORY_HOST:HYPRE_MEMORY_DEVICE));mixedflow::hypre_check(HYPRE_SetExecutionPolicy(profile==mixedflow::FlowProfile::openaccel?HYPRE_EXEC_HOST:HYPRE_EXEC_DEVICE));mixedflow::hypre_check(HYPRE_SetGpuAwareMPI(1));
#else
        mixedflow::hypre_check(HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST));
        mixedflow::hypre_check(HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST));
#endif
        {auto mesh=mixedflow::read(path);mixedflow::Channel channel(mesh,interface_method,profile);channel.view_.controls.high_resolution=high_resolution;
         if(rankId==0) std::printf("Profile: %s\n",profile==mixedflow::FlowProfile::openaccel?"openaccel":"mars");
         if(rankId==0) std::printf("Interface: %s\n",interface_name.c_str());
         if(rankId==0) std::printf("Advection: %s\n",high_resolution?"high-resolution":"upwind");
         bool ok=setup_only || channel.solve(budget,tolerance,output);if(rankId==0) std::printf("%s: mixed %s\n",ok?"PASS":"FAIL",setup_only?"setup only":profile==mixedflow::FlowProfile::openaccel?"reference equation convergence (mass balance and field comparison separate)":"flow convergence (field comparison separate)");result=ok?0:2;}
        mixedflow::hypre_check(HYPRE_Finalize());
    } catch(const std::exception& error) {std::fprintf(stderr,"rank %d mixed flow: %s\n",rankId,error.what());MPI_Abort(MPI_COMM_WORLD,1);}
#ifdef MARS_MIXED_PETSC
    mixedflow::petsc_check(PetscFinalize());
#endif
    MPI_Finalize();return result;
}
