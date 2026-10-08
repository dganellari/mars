#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_ls.h>
#ifdef __CUDACC__
#include <cuda_runtime.h>
#endif
#include "../../../../backend/distributed/unstructured/solvers/mars_hypre_pressure_settings.hpp"
#include "../../../../backend/distributed/unstructured/fem/segregated/mars_segregated_pressure_capture.hpp"
#include <iostream>
#include <dlfcn.h>

namespace frozen=mars::segregated::frozen;
namespace settings=mars::fem::pressure_settings;
using settings::checked;

template<class T> struct Input {
    T* p=nullptr;
    explicit Input(std::vector<T>& h) {
#ifdef __CUDACC__
        frozen::require(cudaMalloc(reinterpret_cast<void**>(&p),h.size()*sizeof(T))==cudaSuccess);
        frozen::require(cudaMemcpy(p,h.data(),h.size()*sizeof(T),cudaMemcpyHostToDevice)==cudaSuccess);
#else
        p=h.data();
#endif
    }
    ~Input() {
#ifdef __CUDACC__
        cudaFree(p);
#endif
    }
    void output(std::vector<T>& h) {
#ifdef __CUDACC__
        frozen::require(cudaMemcpy(h.data(),p,h.size()*sizeof(T),cudaMemcpyDeviceToHost)==cudaSuccess);
#endif
    }
};

int main(int argc,char** argv) {
    checked(hypre_MPI_Init(&argc,&argv));
    int rank=0,ranks=0; checked(hypre_MPI_Comm_rank(hypre_MPI_COMM_WORLD,&rank)); checked(hypre_MPI_Comm_size(hypre_MPI_COMM_WORLD,&ranks));
    try {
        frozen::require(argc==4);
#if defined(HYPRE_USING_GPU) && !defined(__CUDACC__)
        throw std::runtime_error("GPU Hypre replay requires a CUDA compilation");
#endif
#ifdef __CUDACC__
        int devices=0,local=0; hypre_MPI_Comm node;
        checked(MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&node));
        checked(MPI_Comm_rank(node,&local)); checked(MPI_Comm_free(&node));
        frozen::require(cudaGetDeviceCount(&devices)==cudaSuccess && devices>0);
        frozen::require(cudaSetDevice(local%devices)==cudaSuccess);
#endif
        checked(HYPRE_Init());
#ifdef __CUDACC__
        checked(HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE)); checked(HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE));
        checked(HYPRE_SetSpMVUseVendor(0));
        const auto memory=HYPRE_MEMORY_DEVICE;
#else
        checked(HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST)); checked(HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST));
        const auto memory=HYPRE_MEMORY_HOST;
#endif
        {
        const std::filesystem::path capture(argv[1]), configuration(argv[2]), output(argv[3]);
        const auto part=frozen::Part::read(capture/frozen::part_name(rank));
        frozen::require(part.rank==std::uint64_t(rank) && part.ranks==std::uint64_t(ranks));
        std::ifstream cfg(std::filesystem::is_directory(configuration)?configuration/frozen::part_name(rank,".settings"):configuration);
        frozen::require(bool(cfg)); const auto controls=settings::read(cfg);
        frozen::require(controls.count("method") && (controls.at("method")==0 || controls.at("method")==1));
        frozen::require(controls.count("rtol") && controls.count("atol") && controls.at("rtol")>0
            && controls.at("rtol")<1 && controls.at("atol")>=0);
        // A reference configuration must retain the captured acceptance target.
        frozen::require(part.maximum && controls.at("rtol")==part.relative && controls.at("atol")==part.absolute);
#ifdef __CUDACC__
        frozen::require(controls.count("relaxtype") && controls.at("relaxtype")==18
            && controls.count("coarserelax") && controls.at("coarserelax")==18);
#endif
        if (!rank) {
            frozen::require(std::filesystem::create_directory(output));
            std::filesystem::permissions(output,std::filesystem::perms::owner_all);
        }
        checked(hypre_MPI_Barrier(hypre_MPI_COMM_WORLD));
        {
            Dl_info hypre_library{},mpi_library{};
            frozen::require(dladdr(reinterpret_cast<void*>(&HYPRE_BoomerAMGCreate),&hypre_library)!=0);
#ifndef HYPRE_SEQUENTIAL
            frozen::require(dladdr(reinterpret_cast<void*>(&MPI_Init),&mpi_library)!=0);
#endif
            frozen::Writer identity(output/frozen::part_name(rank,".libraries"));
            const std::string text=std::string(hypre_library.dli_fname)+"\n"+(mpi_library.dli_fname?mpi_library.dli_fname:"")+"\n";
            identity.bytes(text.data(),text.size()); identity.finish();
        }
        const auto rows=part.rows();
        std::vector<HYPRE_BigInt> ids(rows),cols(part.nnz);
        std::vector<HYPRE_Int> lengths(rows);
        auto values=part.values,b=part.rhs; std::vector<double> x(rows,0);
        for (std::size_t row=0;row<rows;++row) { ids[row]=HYPRE_BigInt(part.first+row); lengths[row]=part.offsets[row+1]-part.offsets[row]; }
        for (std::size_t k=0;k<cols.size();++k) cols[k]=HYPRE_BigInt(part.map[part.columns[k]]);
        Input<HYPRE_BigInt> d_ids(ids),d_cols(cols); Input<HYPRE_Int> d_lengths(lengths);
        Input<double> d_values(values),d_b(b),d_x(x);
        HYPRE_IJMatrix matrix; HYPRE_IJVector rhs,solution;
        checked(HYPRE_IJMatrixCreate(hypre_MPI_COMM_WORLD,ids.front(),ids.back(),ids.front(),ids.back(),&matrix));
        checked(HYPRE_IJMatrixSetObjectType(matrix,HYPRE_PARCSR));
        checked(HYPRE_IJMatrixInitialize_v2(matrix,memory));
        checked(HYPRE_IJMatrixSetValues(matrix,HYPRE_Int(rows),d_lengths.p,d_ids.p,d_cols.p,d_values.p));
        checked(HYPRE_IJMatrixAssemble(matrix));
        for (auto* vector:{&rhs,&solution}) {
            checked(HYPRE_IJVectorCreate(hypre_MPI_COMM_WORLD,ids.front(),ids.back(),vector));
            checked(HYPRE_IJVectorSetObjectType(*vector,HYPRE_PARCSR)); checked(HYPRE_IJVectorInitialize_v2(*vector,memory));
        }
        checked(HYPRE_IJVectorSetValues(rhs,HYPRE_Int(rows),d_ids.p,d_b.p)); checked(HYPRE_IJVectorAssemble(rhs));
        checked(HYPRE_IJVectorSetValues(solution,HYPRE_Int(rows),d_ids.p,d_x.p)); checked(HYPRE_IJVectorAssemble(solution));
        HYPRE_ParCSRMatrix a; HYPRE_ParVector pb,px;
        checked(HYPRE_IJMatrixGetObject(matrix,reinterpret_cast<void**>(&a)));
        checked(HYPRE_IJVectorGetObject(rhs,reinterpret_cast<void**>(&pb)));
        checked(HYPRE_IJVectorGetObject(solution,reinterpret_cast<void**>(&px)));
        const bool flex=controls.at("method")!=0;
        HYPRE_Solver solver,amg;
        checked((flex?HYPRE_ParCSRFlexGMRESCreate:HYPRE_ParCSRGMRESCreate)(hypre_MPI_COMM_WORLD,&solver));
        checked(HYPRE_BoomerAMGCreate(&amg)); settings::apply(controls,solver,amg,flex);
        checked((flex?HYPRE_ParCSRFlexGMRESSetPrecond:HYPRE_ParCSRGMRESSetPrecond)(solver,HYPRE_BoomerAMGSolve,HYPRE_BoomerAMGSetup,amg));
        checked((flex?HYPRE_ParCSRFlexGMRESSetup:HYPRE_ParCSRGMRESSetup)(solver,a,pb,px));
        std::ostringstream effective; settings::write(effective,settings::snapshot(solver,amg,flex));
        const int error=(flex?HYPRE_ParCSRFlexGMRESSolve:HYPRE_ParCSRGMRESSolve)(solver,a,pb,px);
        const int global_error=HYPRE_GetError(); HYPRE_ClearAllErrors();
        HYPRE_Int iterations=0,converged=0; HYPRE_Real relative=0;
        checked((flex?HYPRE_FlexGMRESGetNumIterations:HYPRE_GMRESGetNumIterations)(solver,&iterations));
        checked((flex?HYPRE_FlexGMRESGetConverged:HYPRE_GMRESGetConverged)(solver,&converged));
        checked((flex?HYPRE_FlexGMRESGetFinalRelativeResidualNorm:HYPRE_GMRESGetFinalRelativeResidualNorm)(solver,&relative));
        checked(HYPRE_IJVectorGetValues(solution,HYPRE_Int(rows),d_ids.p,d_x.p)); d_x.output(x);
        frozen::Writer result(output/frozen::part_name(rank,".solution")); result.array(x); result.finish();
        effective<<std::setprecision(17)<<"result_solve_error "<<error<<"\nresult_global_error "<<global_error
            <<"\nresult_fatal_error "<<bool((error|global_error)&~HYPRE_ERROR_CONV)
            <<"\nresult_iterations "<<iterations<<"\nresult_converged "<<converged<<"\nresult_reported "<<relative<<'\n';
        frozen::Writer info(output/frozen::part_name(rank,".report")); const auto text=effective.str(); info.bytes(text.data(),text.size()); info.finish();
        checked((flex?HYPRE_ParCSRFlexGMRESDestroy:HYPRE_ParCSRGMRESDestroy)(solver)); checked(HYPRE_BoomerAMGDestroy(amg));
        checked(HYPRE_IJVectorDestroy(rhs)); checked(HYPRE_IJVectorDestroy(solution)); checked(HYPRE_IJMatrixDestroy(matrix));
        checked(hypre_MPI_Barrier(hypre_MPI_COMM_WORLD));
        if (!rank) { frozen::Writer done(output/"complete"); const std::string text="mars-pressure-replay-v1\n"+std::to_string(ranks)+"\n"; done.bytes(text.data(),text.size()); done.finish(); }
        }
        checked(HYPRE_Finalize()); checked(hypre_MPI_Finalize());
        // Completion is evidence capture, never a convergence verdict.
        return 0;
    } catch (...) {
        if (!rank) std::cerr<<"ERROR: private pressure replay failed\n";
        hypre_MPI_Abort(hypre_MPI_COMM_WORLD,1); return 1;
    }
}
