#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_ls.h>
#ifdef __CUDACC__
#include <cuda_runtime.h>
#endif
#include "../../../../backend/distributed/unstructured/solvers/mars_hypre_pressure_profile.hpp"
#include "../../../../backend/distributed/unstructured/fem/segregated/mars_segregated_pressure_capture.hpp"
#include <iostream>
#include <dlfcn.h>

namespace frozen=mars::segregated::frozen;
namespace settings=mars::fem::pressure_settings;
using settings::checked;
#include "pressure_replay_recovery.hpp"

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
        checked(hypre_ForceSyncComputeStream());
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
        frozen::require(bool(cfg)); auto controls=settings::read(cfg);
        const double recovery_requested=controls.count("recovery_rounds")?controls.at("recovery_rounds"):0;
        controls.erase("recovery_rounds");
        frozen::require(recovery_requested==0 || recovery_requested==3);
        const int recovery_rounds=int(recovery_requested);
        const int maximum_recovery=recovery::maximum(recovery_rounds,hypre_MPI_COMM_WORLD);
        frozen::require(recovery::maximum(recovery_rounds!=maximum_recovery,hypre_MPI_COMM_WORLD)==0);
        frozen::require(controls.count("method") && (controls.at("method")==0 || controls.at("method")==1));
        frozen::require(controls.count("rtol") && controls.count("atol") && controls.at("rtol")>0
            && controls.at("rtol")<1 && controls.at("atol")>=0);
        // A reference configuration must retain the captured acceptance target.
        frozen::require(part.maximum && controls.at("rtol")==part.relative && controls.at("atol")==part.absolute);
        const bool gpu_profile=controls.count("relax_down") || controls.count("relax_up");
        if (gpu_profile) settings::validate_gpu_profile(controls);
        frozen::require(!recovery_rounds || (gpu_profile && controls.at("maxiter")<=INT32_MAX/4));
        if(recovery_rounds) {
            // IJ must preserve the captured coefficients without duplicate-column sums.
            for(std::size_t row=0;row<part.rows();++row) {
                std::vector<std::int64_t> keys;
                for(int k=part.offsets[row];k<part.offsets[row+1];++k) keys.push_back(part.map[part.columns[k]]);
                std::sort(keys.begin(),keys.end());
                frozen::require(std::adjacent_find(keys.begin(),keys.end())==keys.end());
                frozen::require(part.rhs[row]==0 || part.rhs[row]*part.rhs[row]>=std::numeric_limits<double>::min());
            }
        }
#ifdef __CUDACC__
        frozen::require(gpu_profile || (controls.count("relaxtype") && controls.at("relaxtype")==18
            && controls.count("coarserelax") && controls.at("coarserelax")==18));
        // Recovery requires a GPU-aware MPI runtime, including Hypre's own exchanges.
        if(recovery_rounds) checked(HYPRE_SetGpuAwareMPI(1));
#endif
        if (!rank) {
            frozen::require(std::filesystem::create_directory(output));
            std::filesystem::permissions(output,std::filesystem::perms::owner_all);
            if(recovery_rounds) frozen::require(std::filesystem::create_directory(output/"initial"));
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
        checked(HYPRE_BoomerAMGCreate(&amg));
        if (gpu_profile) {
            settings::Values krylov;
            for (const auto* key:{"method","kdim","miniter","maxiter","rtol","atol"}) krylov[key]=controls.at(key);
            settings::apply(krylov,solver,amg,flex);
            settings::apply_gpu_amg_profile(controls,amg);
        } else settings::apply(controls,solver,amg,flex);
        checked((flex?HYPRE_ParCSRFlexGMRESSetPrecond:HYPRE_ParCSRGMRESSetPrecond)(solver,HYPRE_BoomerAMGSolve,HYPRE_BoomerAMGSetup,amg));
        checked((flex?HYPRE_ParCSRFlexGMRESSetup:HYPRE_ParCSRGMRESSetup)(solver,a,pb,px));
        std::ostringstream effective; settings::write(effective,settings::snapshot(solver,amg,flex));
        const auto initial=recovery::solve(solver,flex,a,pb,px);
        effective<<std::setprecision(17)<<"result_solve_error "<<initial.error<<"\nresult_global_error "<<initial.global
            <<"\nresult_fatal_error "<<bool((initial.error|initial.global)&~HYPRE_ERROR_CONV)
            <<"\nresult_iterations "<<initial.iterations<<"\nresult_converged "<<initial.converged<<"\nresult_reported "<<initial.relative<<'\n';
        if(recovery_rounds) {
            checked(HYPRE_IJVectorGetValues(solution,HYPRE_Int(rows),d_ids.p,d_x.p)); d_x.output(x);
            frozen::Writer before(output/"initial"/frozen::part_name(rank,".solution")); before.array(x); before.finish();
            std::ostringstream trace;
            const auto recovered=recovery::run(reinterpret_cast<hypre_ParCSRMatrix*>(a),reinterpret_cast<hypre_ParVector*>(pb),
                reinterpret_cast<hypre_ParVector*>(px),solver,flex,controls,memory,recovery_rounds,initial,trace);
            const auto restored=settings::snapshot(solver,amg,flex);
            frozen::require(restored.at("rtol")==controls.at("rtol") && restored.at("atol")==controls.at("atol")
                && restored.at("miniter")==controls.at("miniter") && restored.at("maxiter")==controls.at("maxiter"));
            effective<<"recovery_requested_rounds "<<recovery_rounds<<"\nrecovery_rounds "<<recovered.rounds
                <<"\nrecovery_iterations "<<recovered.iterations<<"\nrecovery_stop "<<recovered.stop
                <<"\nrecovery_controls_restored 1\n";
            frozen::Writer history(output/frozen::part_name(rank,".recovery"));
            const auto text=trace.str(); history.bytes(text.data(),text.size()); history.finish();
        }
        checked(HYPRE_IJVectorGetValues(solution,HYPRE_Int(rows),d_ids.p,d_x.p)); d_x.output(x);
        frozen::Writer result(output/frozen::part_name(rank,".solution")); result.array(x); result.finish();
        frozen::Writer info(output/frozen::part_name(rank,".report")); const auto text=effective.str(); info.bytes(text.data(),text.size()); info.finish();
        checked((flex?HYPRE_ParCSRFlexGMRESDestroy:HYPRE_ParCSRGMRESDestroy)(solver)); checked(HYPRE_BoomerAMGDestroy(amg));
        checked(HYPRE_IJVectorDestroy(rhs)); checked(HYPRE_IJVectorDestroy(solution)); checked(HYPRE_IJMatrixDestroy(matrix));
        checked(hypre_MPI_Barrier(hypre_MPI_COMM_WORLD));
        if (!rank) {
            const std::string text="mars-pressure-replay-v1\n"+std::to_string(ranks)+"\n";
            if(recovery_rounds) { frozen::Writer done(output/"initial"/"complete"); done.bytes(text.data(),text.size()); done.finish(); }
            frozen::Writer done(output/"complete"); done.bytes(text.data(),text.size()); done.finish();
        }
        }
        checked(HYPRE_Finalize()); checked(hypre_MPI_Finalize());
        // Completion is evidence capture, never a convergence verdict.
        return 0;
    } catch (...) {
        if (!rank) std::cerr<<"ERROR: private pressure replay failed\n";
        hypre_MPI_Abort(hypre_MPI_COMM_WORLD,1); return 1;
    }
}
